/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
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
#include "complement.hpp"
#include "simple_task.hpp"

namespace npge {

AnchorFinder::AnchorFinder() {
    add_gopt("anchor-size", "anchor size", "ANCHOR_SIZE");
    add_gopt("anchor-fp",
             "Probability of false positive in Bloom filter "
             "(first step of AnchorFinder)", "ANCHOR_FP");
    add_opt("no-palindromes", "eliminate palindromes", true);
    add_opt("only-ori",
            "consider only specified ori; "
            "0 = consider both ori", 0);
    add_opt_rule("anchor-size > 0");
    int max_anchor_size = sizeof(hash_t) * 8 / 2;
    add_opt_rule("anchor-size <= " + TO_S(max_anchor_size));
    add_opt_rule("only-ori >= -1");
    add_opt_rule("only-ori <= 1");
    declare_bs("target", "Blockset to search anchors in");
}

static size_t estimate_length(const BlockSet& bs) {
    typedef std::map<std::string, size_t> GenomeToLength;
    GenomeToLength gtl;
    size_t max_length = 0;
    BOOST_FOREACH (SequencePtr s, bs.seqs()) {
        size_t& genome_length = gtl[s->genome()];
        genome_length += s->size();
        max_length = std::max(max_length, genome_length);
    }
    if (gtl.size() == 1) {
        // consensuses
        // TODO why / 2
        return max_length / 2;
    } else {
        // TODO why * 1.5?
        return max_length / 2 * 3;
    }
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

typedef std::set<hash_t> Possible;

static void test_and_add(SequencePtr s, BloomFilter& filter,
                         size_t anchor_size,
                         Possible& p, int ori_to_add,
                         int only_ori,
                         boost::mutex* mutex) {
    bool prev[3] = {false, false, false};
    hash_t prev_hash[3] = {0, 0, 0};
    Fragment f(s);
    s->make_first_fragment(f, anchor_size, only_ori);
    int Ns = 0; // number of 'N' inside the fragment (* 2 , ori)
    while (only_ori ?
            s->next_fragment_keeping_ori(f) :
            s->next_fragment(f)) {
        bool add = only_ori || f.ori() == ori_to_add;
        hash_t hash;
        if (prev_hash[f.ori() + 1] == 0) {
            hash = f.hash();
            Ns += ns_in_fragment(f); // two times :(
        } else {
            char remove_char = f.raw_at((f.ori() == 1) ?
                                        -1 : anchor_size);
            char add_char = f.at(f.ori() == 1 ? -1 : 0);
            hash = reuse_hash(prev_hash[f.ori() + 1],
                              anchor_size,
                              remove_char, add_char,
                              f.ori() == 1);
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

typedef std::map<hash_t, Block*> HashToBlock;

static void find_blocks(SequencePtr s, size_t anchor_size,
                        const Possible& p,
                        HashToBlock& hash_to_block,
                        int only_ori,
                        boost::mutex* mutex) {
    hash_t prev_hash[3] = {0, 0, 0};
    Fragment f(s);
    s->make_first_fragment(f, anchor_size, only_ori);
    int Ns = 0; // number of 'N' inside the fragment (* 2 , ori)
    while (only_ori ?
            s->next_fragment_keeping_ori(f) :
            s->next_fragment(f)) {
        hash_t hash;
        if (prev_hash[f.ori() + 1] == 0) {
            hash = f.hash();
            Ns += ns_in_fragment(f); // two times :(
        } else {
            char remove_char = f.raw_at((f.ori() == 1) ?
                                        -1 : anchor_size);
            char add_char = f.at(f.ori() == 1 ? -1 : 0);
            hash = reuse_hash(prev_hash[f.ori() + 1],
                              anchor_size,
                              remove_char, add_char,
                              f.ori() == 1);
            if (add_char == 'N') {
                Ns += 1;
            }
            if (remove_char == 'N') {
                Ns -= 1;
            }
        }
        prev_hash[f.ori() + 1] = hash;
        if (Ns == 0 && p.find(hash) != p.end()) {
            hash_t key = hash;
            Block* block;
            if (mutex) {
                mutex->lock();
            }
            HashToBlock::iterator it = hash_to_block.find(key);
            if (it != hash_to_block.end()) {
                block = it->second;
            } else {
                hash_t yek = complement_hash(key, anchor_size);
                if (hash_to_block.find(yek) !=
                        hash_to_block.end() &&
                        !only_ori) {
                    if (mutex) {
                        mutex->unlock();
                    }
                    continue;
                } else {
                    block = hash_to_block[key] = new Block;
                }
            }
            block->insert(new Fragment(f));
            if (mutex) {
                mutex->unlock();
            }
        }
    }
}

static bool check_block(const Block* block,
                        int length, int size) {
    Fragments ff(block->begin(), block->end());
    for (int pos = 0; pos < length; pos++) {
        char c = ff[0]->raw_at(pos);
        for (int f_i = 1; f_i < size; f_i++) {
            if (ff[f_i]->raw_at(pos) != c) {
                return false;
            }
        }
    }
    return true;
}

void AnchorFinder::run_impl() const {
    int anchor_size = opt_value("anchor-size").as<int>();
    bool no_pal = opt_value("no-palindromes").as<bool>();
    int add_ori = no_pal ?
                  -Sequence::FIRST_ORI : Sequence::FIRST_ORI;
    int only_ori = opt_value("only-ori").as<int>();
    boost::mutex* mutex = (workers() == 1) ?
                          0 : new boost::mutex();
    boost::scoped_ptr<boost::mutex> mutex_ptr(mutex);
    size_t length_sum = estimate_length(*block_set());
    if (std::log(length_sum) / std::log(4) > anchor_size) {
        length_sum = std::pow(double(4), double(anchor_size));
    }
    Decimal ep_dec = opt_value("anchor-fp").as<Decimal>();
    double error_prob = ep_dec.to_d();
    Possible possible_anchors;
    {
        BloomFilter filter(length_sum, error_prob);
        Tasks tasks;
        BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
            tasks.push_back(
                boost::bind(test_and_add, s,
                            boost::ref(filter),
                            anchor_size,
                            boost::ref(possible_anchors),
                            add_ori, only_ori, mutex));
        }
        do_tasks(tasks_to_generator(tasks), workers());
    }
    HashToBlock hash_to_block;
    Tasks tasks;
    BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
        tasks.push_back(
            boost::bind(find_blocks, s, anchor_size,
                        boost::ref(possible_anchors),
                        boost::ref(hash_to_block), only_ori,
                        mutex));
    }
    do_tasks(tasks_to_generator(tasks), workers());
    typedef HashToBlock::value_type KeyAndBlock;
    BOOST_FOREACH (const KeyAndBlock& kab, hash_to_block) {
        Block* b = kab.second;
        int size = b->size();
        if (size >= 2 && check_block(b, anchor_size, size)) {
            block_set()->insert(b);
        } else {
            delete b;
        }
    }
}

const char* AnchorFinder::name_impl() const {
    return "Find anchors";
}

}

