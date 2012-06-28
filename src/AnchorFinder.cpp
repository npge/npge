/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <set>
#include <map>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include "po.hpp"

#include "AnchorFinder.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "BloomFilter.hpp"
#include "Exception.hpp"
#include "make_hash.hpp"

namespace bloomrepeats {

AnchorFinder::AnchorFinder():
    min_fragments_(2),
    anchor_size_(ANCHOR_SIZE),
    only_ori_(0),
    workers_(1) {
    set_palindromes_elimination(true);
}

void AnchorFinder::add_options(po::options_description& desc) const {
    desc.add_options()
    ("anchors-min-fragments",
     po::value<size_t>()->default_value(min_fragments()),
     "min number of fragments in a block to accept this block")
    ("anchor-size", po::value<size_t>()->default_value(anchor_size()),
     "anchor size")
    ("no-palindromes", po::bool_switch(), "eliminate palindromes (default)")
    ("palindromes", po::bool_switch(), "do not eliminate palindromes")
    ("only-ori", po::value<int>()->default_value(only_ori()),
     "consider only specified ori; 0 = consider both ori")
    ("workers", po::value<int>()->default_value(workers()),
     "number of threads used to find anchors. "
     "Using >= 2 workers may (very unlikely) cause races, "
     "since bloom filter is not protected by a mutex. "
     "Such a races may cause some anchors not to be found. "
     "The smallest piece of work, passed to a worker, is one sequence. "
     "So it is useless to set workers > sequences.")
   ;
}

void AnchorFinder::apply_options(po::variables_map& vm) {
    set_min_fragments(vm["anchors-min-fragments"].as<size_t>());
    if (vm["anchor-size"].as<size_t>() == 0) {
        throw Exception("'anchor-size' set to 0");
    }
    set_anchor_size(vm["anchor-size"].as<size_t>());
    if (vm["no-palindromes"].as<bool>() && vm["palindromes"].as<bool>()) {
        throw Exception("both 'no-palindromes' and 'palindromes' specified");
    }
    if (vm["no-palindromes"].as<bool>()) {
        set_palindromes_elimination(true);
    } else if (vm["palindromes"].as<bool>()) {
        set_palindromes_elimination(false);
    }
    if (std::abs(vm["only-ori"].as<int>()) > 1) {
        throw Exception("'only-ori' must be -1, 0 or 1");
    }
    set_only_ori(vm["only-ori"].as<int>());
    if (std::abs(vm["workers"].as<int>()) < 1) {
        throw Exception("'workers' number must be >= 1");
    }
    set_workers(vm["workers"].as<int>());
}

void AnchorFinder::add_sequence(SequencePtr s) {
    seqs_.push_back(s);
    if (block_set_) {
        block_set_->add_sequence(s);
    }
}

typedef std::set<size_t> Possible;

static void test_and_add(SequencePtr s, BloomFilter& filter, size_t anchor_size,
                         Possible& p, int ori_to_add, int only_ori,
                         boost::mutex* mutex) {
    bool prev[3] = {false, false, false};
    size_t prev_hash[3] = {0, 0, 0};
    Fragment f(s);
    s->make_first_fragment(f, anchor_size, only_ori);
    while (only_ori ? s->next_fragment_keeping_ori(f) : s->next_fragment(f)) {
        bool add = only_ori || f.ori() == ori_to_add;
        size_t hash;
        if (prev_hash[f.ori() + 1] == 0) {
            hash = f.hash();
        } else {
            char remove_char = f.raw_at(f.ori() == 1 ? -1 : anchor_size);
            char add_char = f.at(f.ori() == 1 ? -1 : 0);
            hash = reuse_hash(prev_hash[f.ori() + 1], anchor_size,
                              remove_char, add_char, f.ori() == 1);
        }
        prev_hash[f.ori() + 1] = hash;
        if (add && filter.test_and_add(hash) || !add && filter.test(hash)) {
            if (!prev[f.ori() + 1]) {
                prev[f.ori() + 1] = true;
                if (mutex) {
                    mutex->lock();
                }
                p.insert(hash);
                if (mutex) {
                    mutex->unlock();
                }
            }
        } else {
            prev[f.ori() + 1] = false;
        }
    }
}

typedef std::map<std::string, BlockPtr> StrToBlock;

static void find_blocks(SequencePtr s, size_t anchor_size, const Possible& p,
                        StrToBlock& str_to_block, int only_ori,
                        boost::mutex* mutex) {
    size_t prev_hash[3] = {0, 0, 0};
    Fragment f(s);
    s->make_first_fragment(f, anchor_size, only_ori);
    while (only_ori ? s->next_fragment_keeping_ori(f) : s->next_fragment(f)) {
        size_t hash;
        if (prev_hash[f.ori() + 1] == 0) {
            hash = f.hash();
        } else {
            char remove_char = f.raw_at(f.ori() == 1 ? -1 : anchor_size);
            char add_char = f.at(f.ori() == 1 ? -1 : 0);
            hash = reuse_hash(prev_hash[f.ori() + 1], anchor_size,
                              remove_char, add_char, f.ori() == 1);
        }
        prev_hash[f.ori() + 1] = hash;
        if (p.find(hash) != p.end()) {
            std::string key = f.str();
            BlockPtr block;
            if (mutex) {
                mutex->lock();
            }
            if (str_to_block.find(key) != str_to_block.end()) {
                block = str_to_block[key];
            } else {
                std::string complement_key = key;
                complement(complement_key);
                if (str_to_block.find(complement_key) != str_to_block.end() &&
                        !only_ori) {
                    if (mutex) {
                        mutex->unlock();
                    }
                    continue;
                } else {
                    block = str_to_block[key] = Block::create_new();
                }
            }
            block->insert(Fragment::create_new(f));
            if (mutex) {
                mutex->unlock();
            }
        }
    }
}

typedef boost::function<void()> Task;
typedef std::vector<Task> Tasks;

void process_some_seqs(Tasks& tasks, boost::mutex* mutex) {
    while (true) {
        Task task;
        if (mutex) {
            mutex->lock();
        }
        if (!tasks.empty()) {
            task = tasks.back();
            tasks.pop_back();
        }
        if (mutex) {
            mutex->unlock();
        }
        if (!task.empty()) {
            task();
        } else {
            break;
        }
    }
}

typedef boost::thread Thread;
typedef boost::shared_ptr<Thread> ThreadPtr;
typedef std::vector<ThreadPtr> Threads;

static void do_tasks(Tasks& tasks, int workers, boost::mutex* mutex) {
    Threads threads;
    for (int i = 1; i < workers; i++) {
        threads.push_back(boost::make_shared<boost::thread>(
                              boost::bind(process_some_seqs, boost::ref(tasks),
                                          mutex)));
    }
    process_some_seqs(tasks, mutex);
    BOOST_FOREACH (ThreadPtr thread, threads) {
        thread->join();
    }
}

void AnchorFinder::run() {
    boost::mutex* mutex = workers_ == 1 ? 0 : new boost::mutex();
    if (!anchor_handler_) {
        return;
    }
    size_t length_sum = 0;
    BOOST_FOREACH (SequencePtr s, seqs_) {
        length_sum += s->size();
    }
    if (std::log(length_sum) / std::log(4) > anchor_size_) {
        length_sum = std::pow(4, anchor_size_);
    }
    float error_prob = 1.0 / length_sum;
    std::set<size_t> possible_anchors;
    {
        BloomFilter filter(length_sum, error_prob);
        Tasks tasks;
        BOOST_FOREACH (SequencePtr s, seqs_) {
            tasks.push_back(
                boost::bind(test_and_add, s, boost::ref(filter),
                            anchor_size_, boost::ref(possible_anchors),
                            add_ori_, only_ori_, mutex));
        }
        do_tasks(tasks, workers_, mutex);
    }
    StrToBlock str_to_block;
    Tasks tasks;
    BOOST_FOREACH (SequencePtr s, seqs_) {
        tasks.push_back(
            boost::bind(find_blocks, s, anchor_size_,
                        boost::ref(possible_anchors),
                        boost::ref(str_to_block), only_ori_, mutex));
    }
    do_tasks(tasks, workers_, mutex);
    BOOST_FOREACH (const StrToBlock::value_type& key_and_block, str_to_block) {
        BlockPtr block = key_and_block.second;
        if (block->size() >= min_fragments_) {
            anchor_handler_(block);
        }
    }
    delete mutex;
}

void AnchorFinder::set_block_set(BlockSetPtr block_set) {
    anchor_handler_ = boost::bind(&BlockSet::insert, block_set.get(), _1);
    BOOST_FOREACH (SequencePtr seq, seqs_) {
        block_set->add_sequence(seq);
    }
    block_set_ = block_set;
}

bool AnchorFinder::palindromes_elimination() const {
    return add_ori_ == -Sequence::FIRST_ORI;
}

void AnchorFinder::set_palindromes_elimination(bool eliminate) {
    add_ori_ = eliminate ? -Sequence::FIRST_ORI : Sequence::FIRST_ORI;
}

}

