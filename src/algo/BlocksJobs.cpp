/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "BlocksJobs.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "thread_pool.hpp"

namespace npge {

typedef std::vector<Block*> BlocksVector;

ThreadData::ThreadData() {
}

ThreadData::~ThreadData() {
}

class BlockGroup : public ReusingThreadGroup {
public:
    BlockGroup(const BlocksJobs* jobs):
        jobs_(jobs), bs_i_(0) {
        std::string block_set_name = jobs->block_set_name();
        BlockSetPtr target = jobs->get_bs(block_set_name);
        BlocksVector _(target->begin(), target->end());
        bs_.swap(_);
        set_workers(jobs->workers());
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl();

    void perform_impl() {
        jobs_->change_blocks(bs_);
        jobs_->initialize_work();
        ReusingThreadGroup::perform_impl();
        jobs_->finish_work();
    }

private:
    const BlocksJobs* jobs_;
    BlocksVector bs_;
    int bs_i_;
};

class BlockWorker : public ThreadWorker {
public:
    ThreadData* data_;

    BlockWorker(const BlocksJobs* jobs, BlockGroup* group):
        ThreadWorker(group), jobs_(jobs) {
        data_ = jobs_->before_thread();
    }

    void work_impl() {
        jobs_->initialize_thread(data_);
        ThreadWorker::work_impl();
        jobs_->finish_thread(data_);
    }

    ~BlockWorker() {
        // this is called from main thread
        jobs_->after_thread(data_);
    }

private:
    const BlocksJobs* jobs_;
};

class BlockTask : public ThreadTask {
public:
    BlockTask(Block* block, const BlocksJobs* jobs, BlockWorker* worker):
        ThreadTask(worker), block_(block), jobs_(jobs) {
    }

    void run_impl() {
        BlockWorker* w = boost::polymorphic_downcast<BlockWorker*>(worker());
        jobs_->process_block(block_, w->data_);
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
    block_set_name_(block_set_name) {
}

struct BlockCompareName2 {
    bool operator()(const Block* b1, const Block* b2) const {
        typedef boost::tuple<int, int, const std::string&> Tie;
        return Tie(-b1->size(), -b1->alignment_length(), b1->name()) <
               Tie(-b2->size(), -b2->alignment_length(), b2->name());
    }
};

void BlocksJobs::sort_blocks(std::vector<Block*>& blocks) const {
    std::sort(blocks.begin(), blocks.end(), BlockCompareName2());
}

void BlocksJobs::change_blocks(BlocksVector& blocks) const {
    change_blocks_impl(blocks);
}

void BlocksJobs::initialize_work() const {
    initialize_work_impl();
}

ThreadData* BlocksJobs::before_thread() const {
    return before_thread_impl();
}

void BlocksJobs::initialize_thread(ThreadData* data) const {
    initialize_thread_impl(data);
}

void BlocksJobs::process_block(Block* block, ThreadData* data) const {
    check_interruption();
    process_block_impl(block, data);
}

void BlocksJobs::finish_thread(ThreadData* data) const {
    finish_thread_impl(data);
}

void BlocksJobs::after_thread(ThreadData* data) const {
    after_thread_impl(data);
    delete data;
}

void BlocksJobs::finish_work() const {
    finish_work_impl();
}

void BlocksJobs::run_impl() const {
    BlockGroup block_group(this);
    block_group.perform();
}

void BlocksJobs::change_blocks_impl(BlocksVector& blocks) const {
    sort_blocks(blocks);
}

void BlocksJobs::initialize_work_impl() const {
}

ThreadData* BlocksJobs::before_thread_impl() const {
    return 0;
}

void BlocksJobs::initialize_thread_impl(ThreadData* data) const {
}

void BlocksJobs::process_block_impl(Block* block, ThreadData* data) const {
}

void BlocksJobs::after_thread_impl(ThreadData* data) const {
}

void BlocksJobs::finish_thread_impl(ThreadData* data) const {
}

void BlocksJobs::finish_work_impl() const {
}

}

