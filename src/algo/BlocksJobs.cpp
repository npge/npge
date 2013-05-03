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

BlocksJobs::BlocksJobs(const std::string& block_set_name):
    block_set_name_(block_set_name)
{ }

void BlocksJobs::change_blocks(BlocksVector& blocks) const {
    change_blocks_impl(blocks);
}

bool BlocksJobs::initialize_thread() const {
    return initialize_thread_impl();
}

bool BlocksJobs::apply_to_block(Block* block) const {
    return apply_to_block_impl(block);
}

bool BlocksJobs::finish_thread() const {
    return finish_thread_impl();
}

bool BlocksJobs::finish_work() const {
    return finish_work_impl();
}

typedef BlocksVector::const_iterator It;

static void process_blocks(It& it, const It& end, boost::mutex& mutex,
                           const BlocksJobs* jobs, bool& result) {
    bool changed = false;
    changed |= jobs->initialize_thread();
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
    changed |= jobs->finish_thread();
    if (changed) {
        boost::mutex::scoped_lock lock(mutex);
        result = true;
    }
}

bool BlocksJobs::run_impl() const {
    bool result = false;
    BlockSetPtr target = get_bs(block_set_name());
    BlocksVector bs(target->begin(), target->end());
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
    if (finish_work()) {
        result = true;
    }
    return result;
}

void BlocksJobs::change_blocks_impl(BlocksVector& blocks) const
{ }

bool BlocksJobs::initialize_thread_impl() const {
    return false;
}

bool BlocksJobs::finish_thread_impl() const {
    return false;
}

bool BlocksJobs::finish_work_impl() const {
    return false;
}

}

