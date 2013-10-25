/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>

#include "thread_group.hpp"

namespace bloomrepeats {

// No threads beyond this file!

ThreadTask::ThreadTask(ThreadWorker* worker):
    worker_(worker), thread_group_(worker_->thread_group())
{ }

ThreadTask::~ThreadTask()
{ }

void ThreadTask::run() {
    run_impl();
}

ThreadWorker* ThreadTask::worker() const {
    return worker_;
}

ThreadGroup* ThreadTask::thread_group() const {
    return thread_group_;
}

ThreadWorker::ThreadWorker(ThreadGroup* thread_group):
    thread_group_(thread_group)
{ }

ThreadWorker::~ThreadWorker()
{ }

void ThreadWorker::work() {
    work_impl();
}

void ThreadWorker::run(ThreadTask* task) {
    run_impl(task);
}

ThreadGroup* ThreadWorker::thread_group() const {
    return thread_group_;
}

void ThreadWorker::work_impl() {
    while (true) {
        boost::scoped_ptr<ThreadTask> task(thread_group()->create_task(this));
        if (task) {
            run(task.get());
        } else {
            break;
        }
    }
}

void ThreadWorker::run_impl(ThreadTask* task) {
    task->run();
}

struct ThreadGroup::Impl {
    boost::mutex mutex_;
    int workers_;
};

ThreadGroup::ThreadGroup():
    impl_(new Impl)
{ }

ThreadGroup::~ThreadGroup() {
    delete impl_;
    impl_ = 0;
}

void ThreadGroup::perform(int workers) {
    perform_impl(workers);
}

ThreadTask* ThreadGroup::create_task(ThreadWorker* worker) {
    if (impl_->workers_ == 1) {
        return create_task_impl(worker);
    } else {
        boost::mutex::scoped_lock lock(impl_->mutex_);
        return create_task_impl(worker);
    }
}

ThreadWorker* ThreadGroup::create_worker() {
    return create_worker_impl();
}

void ThreadGroup::perform_impl(int workers) {
    impl_->workers_ = workers;
    boost::thread_group threads;
    typedef boost::shared_ptr<ThreadWorker> ThreadWorkerPtr;
    std::vector<ThreadWorkerPtr> workers_;
    for (int i = 1; i < workers; i++) {
        ThreadWorker* worker = create_worker();
        workers_.push_back(ThreadWorkerPtr(worker));
        threads.create_thread(boost::bind(&ThreadWorker::work, worker));
    }
    ThreadWorker* worker = create_worker();
    workers_.push_back(ThreadWorkerPtr(worker));
    worker->work();
    threads.join_all();
    // workers are deleted here
}

ThreadWorker* ThreadGroup::create_worker_impl() {
    return new ThreadWorker(this);
}

}

