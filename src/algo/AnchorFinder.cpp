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
#include <boost/scoped_ptr.hpp>
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

AnchorFinder::AnchorFinder() {
    add_opt("anchor-size", "anchor size", int(ANCHOR_SIZE));
    add_opt("no-palindromes", "eliminate palindromes", true);
    add_opt("only-ori",
            "consider only specified ori; 0 = consider both ori", 0);
    add_opt_rule("anchor-size > 0");
    add_opt_rule("only-ori >= -1");
    add_opt_rule("only-ori <= 1");
    declare_bs("target", "Blockset in which anchors are searched");
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
        if (Ns == 0 && ((add && filter.test_and_add(hash)) ||
                        (!add && filter.test(hash)))) {
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

void AnchorFinder::run_impl() const {
    int anchor_size = opt_value("anchor-size").as<int>();
    bool no_palindromes = opt_value("no-palindromes").as<bool>();
    int add_ori = no_palindromes ?
                  -Sequence::FIRST_ORI : Sequence::FIRST_ORI;
    int only_ori = opt_value("only-ori").as<int>();
    boost::mutex* mutex = workers() == 1 ? 0 : new boost::mutex();
    boost::scoped_ptr<boost::mutex> mutex_ptr(mutex);
    size_t length_sum = 0;
    BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
        length_sum += s->size();
    }
    if (std::log(length_sum) / std::log(4) > anchor_size) {
        length_sum = std::pow(double(4), double(anchor_size));
    }
    float error_prob = 1.0 / length_sum;
    std::set<size_t> possible_anchors;
    {
        BloomFilter filter(length_sum, error_prob);
        Tasks tasks;
        BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
            tasks.push_back(
                boost::bind(test_and_add, s, boost::ref(filter),
                            anchor_size, boost::ref(possible_anchors),
                            add_ori, only_ori, mutex));
        }
        do_tasks(tasks_to_generator(tasks), workers());
    }
    StrToBlock str_to_block;
    Tasks tasks;
    BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
        tasks.push_back(
            boost::bind(find_blocks, s, anchor_size,
                        boost::ref(possible_anchors),
                        boost::ref(str_to_block), only_ori, mutex));
    }
    do_tasks(tasks_to_generator(tasks), workers());
    BOOST_FOREACH (const StrToBlock::value_type& key_and_block, str_to_block) {
        Block* block = key_and_block.second;
        if (block->size() >= 2) {
            block_set()->insert(block);
        } else {
            delete block;
        }
    }
}

const char* AnchorFinder::name_impl() const {
    return "Find anchors";
}

}

