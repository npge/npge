/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "BlocksJobs.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Meta.hpp"
#include "thread_pool.hpp"
#include "cast.hpp"

namespace npge {

typedef std::vector<Block*> BlocksVector;

WorkData::WorkData() {
}

WorkData::~WorkData() {
}

ThreadData::ThreadData():
    work_data_(0) {
}

ThreadData::~ThreadData() {
}

class BlockGroup : public ReusingThreadGroup {
public:
    BlockGroup(const BlocksJobs* jobs):
        jobs_(jobs), bs_i_(0), work_data_(0) {
        std::string block_set_name = jobs->block_set_name();
        BlockSetPtr target = jobs->get_bs(block_set_name);
        BlocksVector _(target->begin(), target->end());
        bs_.swap(_);
        set_workers(jobs->workers());
        const Meta* meta = jobs->meta();
        AnyAs big = meta->get_opt("BLOCKS_IN_GROUP", 1);
        blocks_in_group_ = big.as<int>();
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl();

    void perform_impl() {
        jobs_->change_blocks(bs_);
        jobs_->initialize_work();
        work_data_ = jobs_->before_work();
        ReusingThreadGroup::perform_impl();
        jobs_->finish_work();
        jobs_->after_work(work_data_);
        delete work_data_;
    }

private:
    const BlocksJobs* jobs_;
    WorkData* work_data_;
    BlocksVector bs_;
    int bs_i_;
    int blocks_in_group_;
};

class BlockWorker : public ThreadWorker {
public:
    ThreadData* data_;

    BlockWorker(const BlocksJobs* jobs, WorkData* work_data,
                BlockGroup* group):
        ThreadWorker(group),
        keeper_(jobs->meta()),
        jobs_(jobs),
        work_completed_(false) {
        data_ = jobs_->before_thread();
        if (data_) {
            data_->work_data_ = work_data;
        }
    }

    void work_impl() {
        jobs_->initialize_thread(data_);
        ThreadWorker::work_impl();
        jobs_->finish_thread(data_);
        work_completed_ = true;
    }

    ~BlockWorker() {
        // this is called from main thread
        if (work_completed_) {
            jobs_->after_thread(data_);
        }
    }

private:
    MetaThreadKeeper keeper_;
    const BlocksJobs* jobs_;
    bool work_completed_;
};

class BlockTask : public ThreadTask {
public:
    BlockTask(const BlocksJobs* jobs, BlockWorker* worker):
        ThreadTask(worker), jobs_(jobs) {
    }

    void run_impl() {
        BlockWorker* w = D_CAST<BlockWorker*>(worker());
        BOOST_FOREACH (Block* block, blocks_) {
            jobs_->process_block(block, w->data_);
        }
    }

    Blocks blocks_;
    const BlocksJobs* jobs_;
};

class OneBlockTask : public ThreadTask {
public:
    OneBlockTask(Block* block, const BlocksJobs* jobs,
                 BlockWorker* worker):
        ThreadTask(worker), block_(block), jobs_(jobs) {
    }

    void run_impl() {
        BlockWorker* w = D_CAST<BlockWorker*>(worker());
        jobs_->process_block(block_, w->data_);
    }

    Block* block_;
    const BlocksJobs* jobs_;
};

ThreadTask* BlockGroup::create_task_impl(ThreadWorker* worker) {
    if (bs_i_ < bs_.size()) {
        BlockWorker* w = D_CAST<BlockWorker*>(worker);
        if (workers() == 1) {
            // all blocks to one worker
            BlockTask* task = new BlockTask(jobs_, w);
            task->blocks_.swap(bs_);
            return task;
        }
        int tasks = bs_.size() - bs_i_;
        int taks_per_worker = tasks / workers();
        taks_per_worker = std::max(taks_per_worker, 1);
        int n = std::min(taks_per_worker, blocks_in_group_);
        if (n == 1) {
            Block* block = bs_[bs_i_];
            bs_i_ += 1;
            return new OneBlockTask(block, jobs_, w);
        } else {
            BlockTask* task = new BlockTask(jobs_, w);
            task->blocks_.reserve(n);
            for (int i = 0; i < n; i++) {
                Block* block = bs_[bs_i_];
                task->blocks_.push_back(block);
                bs_i_ += 1;
            }
            return task;
        }
    } else {
        return 0;
    }
}

ThreadWorker* BlockGroup::create_worker_impl() {
    return new BlockWorker(jobs_, work_data_, this);
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

WorkData* BlocksJobs::before_work() const {
    return before_work_impl();
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

void BlocksJobs::after_work(WorkData* work_data) const {
    after_work_impl(work_data);
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

WorkData* BlocksJobs::before_work_impl() const {
    return 0;
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

void BlocksJobs::after_work_impl(WorkData*) const {
}

}

