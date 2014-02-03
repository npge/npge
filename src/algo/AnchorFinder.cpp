/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <set>
#include <map>
#include "boost-xtime.hpp"
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
#include "simple_task.hpp"

namespace bloomrepeats {

AnchorFinder::AnchorFinder():
    anchor_size_(ANCHOR_SIZE),
    only_ori_(0) {
    set_palindromes_elimination(true);
}

void AnchorFinder::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("anchor-size", po::value<size_t>()->default_value(anchor_size()),
     "anchor size")
    ("no-palindromes",
     po::value<bool>()->default_value(palindromes_elimination()),
     "eliminate palindromes")
    ("only-ori", po::value<int>()->default_value(only_ori()),
     "consider only specified ori; 0 = consider both ori")
   ;
}

void AnchorFinder::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("anchor-size")) {
        if (vm["anchor-size"].as<size_t>() == 0) {
            throw Exception("'anchor-size' set to 0");
        }
        set_anchor_size(vm["anchor-size"].as<size_t>());
    }
    if (vm.count("no-palindromes")) {
        set_palindromes_elimination(vm["no-palindromes"].as<bool>());
    }
    if (std::abs(vm["only-ori"].as<int>()) > 1) {
        throw Exception("'only-ori' must be -1, 0 or 1");
    }
    set_only_ori(vm["only-ori"].as<int>());
}

static int ns_in_fragment(const Fragment& f) {
    int result = 0;
    for (int i = 0; i < f.length(); i++) {
        if (f.raw_at(i) == 'N') {
            result += 1;
        }
    }
    return result;
}

typedef std::set<size_t> Possible;

static void test_and_add(SequencePtr s, BloomFilter& filter, size_t anchor_size,
                         Possible& p, int ori_to_add, int only_ori,
                         boost::mutex* mutex) {
    bool prev[3] = {false, false, false};
    size_t prev_hash[3] = {0, 0, 0};
    Fragment f(s);
    s->make_first_fragment(f, anchor_size, only_ori);
    int Ns = 0; // number of 'N' inside the fragment (* 2 , ori)
    while (only_ori ? s->next_fragment_keeping_ori(f) : s->next_fragment(f)) {
        bool add = only_ori || f.ori() == ori_to_add;
        size_t hash;
        if (prev_hash[f.ori() + 1] == 0) {
            hash = f.hash();
            Ns += ns_in_fragment(f); // two times :(
        } else {
            char remove_char = f.raw_at(f.ori() == 1 ? -1 : anchor_size);
            char add_char = f.at(f.ori() == 1 ? -1 : 0);
            hash = reuse_hash(prev_hash[f.ori() + 1], anchor_size,
                              remove_char, add_char, f.ori() == 1);
            if (add_char == 'N') {
                Ns += 1;
            }
            if (remove_char == 'N') {
                Ns -= 1;
            }
        }
        prev_hash[f.ori() + 1] = hash;
        if (Ns == 0 && (add && filter.test_and_add(hash) ||
                    !add && filter.test(hash))) {
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

typedef std::map<std::string, Block*> StrToBlock;

static void find_blocks(SequencePtr s, size_t anchor_size, const Possible& p,
                        StrToBlock& str_to_block, int only_ori,
                        boost::mutex* mutex) {
    size_t prev_hash[3] = {0, 0, 0};
    Fragment f(s);
    s->make_first_fragment(f, anchor_size, only_ori);
    int Ns = 0; // number of 'N' inside the fragment (* 2 , ori)
    while (only_ori ? s->next_fragment_keeping_ori(f) : s->next_fragment(f)) {
        size_t hash;
        if (prev_hash[f.ori() + 1] == 0) {
            hash = f.hash();
            Ns += ns_in_fragment(f); // two times :(
        } else {
            char remove_char = f.raw_at(f.ori() == 1 ? -1 : anchor_size);
            char add_char = f.at(f.ori() == 1 ? -1 : 0);
            hash = reuse_hash(prev_hash[f.ori() + 1], anchor_size,
                              remove_char, add_char, f.ori() == 1);
            if (add_char == 'N') {
                Ns += 1;
            }
            if (remove_char == 'N') {
                Ns -= 1;
            }
        }
        prev_hash[f.ori() + 1] = hash;
        if (Ns == 0 && p.find(hash) != p.end()) {
            std::string key = f.str();
            Block* block;
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
                    block = str_to_block[key] = new Block();
                }
            }
            block->insert(new Fragment(f));
            if (mutex) {
                mutex->unlock();
            }
        }
    }
}

bool AnchorFinder::run_impl() const {
    boost::mutex* mutex = workers() == 1 ? 0 : new boost::mutex();
    size_t length_sum = 0;
    BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
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
        BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
            tasks.push_back(
                boost::bind(test_and_add, s, boost::ref(filter),
                            anchor_size_, boost::ref(possible_anchors),
                            add_ori_, only_ori_, mutex));
        }
        do_tasks(tasks_to_generator(tasks), workers());
    }
    StrToBlock str_to_block;
    Tasks tasks;
    BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
        tasks.push_back(
            boost::bind(find_blocks, s, anchor_size_,
                        boost::ref(possible_anchors),
                        boost::ref(str_to_block), only_ori_, mutex));
    }
    do_tasks(tasks_to_generator(tasks), workers());
    bool result = false;
    BOOST_FOREACH (const StrToBlock::value_type& key_and_block, str_to_block) {
        Block* block = key_and_block.second;
        if (block->size() >= 2) {
            block_set()->insert(block);
            result |= true;
        } else {
            delete block;
        }
    }
    delete mutex;
    return result;
}

const char* AnchorFinder::name_impl() const {
    return "Find anchors";
}

bool AnchorFinder::palindromes_elimination() const {
    return add_ori_ == -Sequence::FIRST_ORI;
}

void AnchorFinder::set_palindromes_elimination(bool eliminate) {
    add_ori_ = eliminate ? -Sequence::FIRST_ORI : Sequence::FIRST_ORI;
}

}

