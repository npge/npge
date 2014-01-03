/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/cast.hpp>

#include "BlocksJobs.hpp"
#include "BlockSet.hpp"
#include "thread_group.hpp"

namespace bloomrepeats {

typedef std::vector<Block*> BlocksVector;

class BlockGroup : public ThreadGroup {
public:
    bool changed_;

    BlockGroup(const BlocksJobs* jobs, const std::string& block_set_name):
        changed_(false), jobs_(jobs), bs_i_(0) {
        BlockSetPtr target = jobs->get_bs(block_set_name);
        BlocksVector _(target->begin(), target->end());
        bs_.swap(_);
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl();

    void perform_impl(int workers) {
        changed_ |= jobs_->change_blocks(bs_);
        ThreadGroup::perform_impl(workers);
        changed_ |= jobs_->finish_work();
    }

private:
    const BlocksJobs* jobs_;
    BlocksVector bs_;
    int bs_i_;
};

class BlockWorker : public ThreadWorker {
public:
    bool changed_;

    BlockWorker(const BlocksJobs* jobs, BlockGroup* group):
        ThreadWorker(group), jobs_(jobs), changed_(false)
    { }

    void work_impl() {
        changed_ |= jobs_->initialize_thread();
        ThreadWorker::work_impl();
        changed_ |= jobs_->finish_thread();
    }

    ~BlockWorker() {
        // this is called from main thread
        BlockGroup* tg;
        tg = boost::polymorphic_downcast<BlockGroup*>(thread_group());
        tg->changed_ |= changed_;
    }

private:
    const BlocksJobs* jobs_;
};

class BlockTask : public ThreadTask {
public:
    BlockTask(Block* block, const BlocksJobs* jobs, BlockWorker* worker):
        ThreadTask(worker), block_(block), jobs_(jobs)
    { }

    void run_impl() {
        bool changed = jobs_->apply_to_block(block_);
        BlockWorker* w = boost::polymorphic_downcast<BlockWorker*>(worker());
        w->changed_ |= changed;
    }

private:
    Block* block_;
    const BlocksJobs* jobs_;
};

ThreadTask* BlockGroup::create_task_impl(ThreadWorker* worker) {
    if (bs_i_ < bs_.size()) {
        Block* block = bs_[bs_i_];
        bs_i_ += 1;
        BlockWorker* w = boost::polymorphic_downcast<BlockWorker*>(worker);
        return new BlockTask(block, jobs_, w);
    } else {
        return 0;
    }
}

ThreadWorker* BlockGroup::create_worker_impl() {
    return new BlockWorker(jobs_, this);
}

BlocksJobs::BlocksJobs(const std::string& block_set_name):
    block_set_name_(block_set_name)
{ }

bool BlocksJobs::change_blocks(BlocksVector& blocks) const {
    change_blocks_impl(blocks);
}

bool BlocksJobs::initialize_thread() const {
    return initialize_thread_impl();
}

bool BlocksJobs::finish_thread() const {
    return finish_thread_impl();
}

bool BlocksJobs::finish_work() const {
    return finish_work_impl();
}

bool BlocksJobs::run_impl() const {
    BlockGroup block_group(this, block_set_name());
    block_group.perform(workers());
    return block_group.changed_;
}

bool BlocksJobs::change_blocks_impl(BlocksVector& blocks) const {
    return false;
}

bool BlocksJobs::initialize_thread_impl() const {
    return false;
}

bool BlocksJobs::apply_to_block_impl(Block* block) const {
    return false;
}

bool BlocksJobs::finish_thread_impl() const {
    return false;
}

bool BlocksJobs::finish_work_impl() const {
    return false;
}

}

