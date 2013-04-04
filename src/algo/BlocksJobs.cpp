/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>

#include "BlocksJobs.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

typedef std::vector<Block*> BlocksVector;

void BlocksJobs::change_blocks(BlocksVector& blocks) const {
    change_blocks_impl(blocks);
}

bool BlocksJobs::apply_to_block(Block* block) const {
    return apply_to_block_impl(block);
}

typedef BlocksVector::const_iterator It;

static void process_blocks(It& it, const It& end, boost::mutex& mutex,
                           const BlocksJobs* jobs, bool& result) {
    bool changed = false;
    while (true) {
        Block* block;
        {
            boost::mutex::scoped_lock lock(mutex);
            if (it != end) {
                block = *it;
                ++it;
            } else {
                break;
            }
        }
        changed |= jobs->apply_to_block(block);
    }
    if (changed) {
        boost::mutex::scoped_lock lock(mutex);
        result = true;
    }
}

bool BlocksJobs::run_impl() const {
    bool result = false;
    BlocksVector bs(block_set()->begin(), block_set()->end());
    change_blocks(bs);
    It it = bs.begin();
    const It end = bs.end();
    boost::mutex mutex;
    boost::thread_group threads;
    for (int i = 1; i < workers(); i++) {
        using namespace boost;
        threads.create_thread(bind(process_blocks, ref(it), ref(end),
                                   ref(mutex), this, ref(result)));
    }
    process_blocks(it, end, mutex, this, result);
    threads.join_all();
    return result;
}

void BlocksJobs::change_blocks_impl(BlocksVector& blocks) const
{ }

}

