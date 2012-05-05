/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <map>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include "AnchorFinder.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BloomFilter.hpp"
#include "make_hash.hpp"

namespace bloomrepeats {

AnchorFinder::AnchorFinder():
    anchor_size_(ANCHOR_SIZE)
{ }

void AnchorFinder::add_sequnce(SequencePtr s) {
    seqs_.push_back(s);
}

typedef std::set<size_t> Possible;

const size_t HASH_MUL = 1484954565;

static void test_and_add(SequencePtr s, BloomFilter& filter,
                         size_t anchor_size, Possible& p) {
    for (size_t start = anchor_size;; start++) {
        size_t length = anchor_size;
        const char* data = s->get(start, length);
        if (length == anchor_size) {
            for (int ori = -1; ori <= 1; ori += 2) {
                if (filter.test_and_add(data, length, ori)) {
                    p.insert(make_hash(HASH_MUL, data, length, ori));
                }
            }
        } else {
            break;
        }
    }
}

typedef std::map<std::string, BlockPtr> StrToBlock;

static void find_blocks(SequencePtr s, size_t anchor_size,
                        const Possible& p, StrToBlock& str_to_block) {
    for (size_t start = anchor_size;; start++) {
        size_t length = anchor_size;
        const char* data = s->get(start, length);
        if (length == anchor_size) {
            for (int ori = -1; ori <= 1; ori += 2) {
                if (p.find(make_hash(HASH_MUL, data, length, ori)) != p.end()) {
                    FragmentPtr fragment = boost::make_shared<Fragment>(s,
                                           start, start + length - 1, ori);
                    std::string key = fragment->str();
                    BlockPtr block;
                    if (str_to_block.find(key) == str_to_block.end()) {
                        block = str_to_block[key] = boost::make_shared<Block>();
                    } else {
                        block = str_to_block[key];
                    }
                    block->insert(fragment);
                }
            }
        } else {
            break;
        }
    }
}

void AnchorFinder::run() {
    if (!anchor_handler_) {
        return;
    }
    size_t length_sum = 0;
    BOOST_FOREACH (SequencePtr s, seqs_) {
        length_sum += s->approximate_size();
    }
    length_sum *= 2; // ori = 1 | -1
    float error_prob = 1.0 / length_sum;
    std::set<size_t> possible_anchors;
    {
        BloomFilter filter(length_sum, error_prob);
        BOOST_FOREACH (SequencePtr s, seqs_) {
            test_and_add(s, filter, anchor_size_, possible_anchors);
        }
    }
    StrToBlock str_to_block;
    BOOST_FOREACH (SequencePtr s, seqs_) {
        find_blocks(s, anchor_size_, possible_anchors, str_to_block);
    }
    BOOST_FOREACH (const StrToBlock::value_type& key_and_block, str_to_block) {
        BlockPtr block = key_and_block.second;
        if (block->size() >= 2) {
            anchor_handler_(block);
        }
    }
}

}

