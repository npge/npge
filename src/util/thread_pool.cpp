/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <set>
#include "boost-xtime.hpp"
#include <boost/asio/io_service.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

#include "thread_pool.hpp"

namespace npge {

// http://stackoverflow.com/a/19500405

typedef boost::mutex Mutex;
typedef boost::mutex::scoped_lock Lock;
typedef boost::condition_variable Condition;
typedef boost::asio::io_service IoService;
typedef boost::asio::io_service::work Work;
typedef std::set<ThreadWorker*> WorkersSet;

struct ThreadPoolImpl {
    IoService io_service_;
    boost::thread_group threads_;
    Work work_;

    WorkersSet finished_;
    Mutex finished_mutex_;
    Condition finished_condition_;

    ThreadPoolImpl():
        work_(io_service_) {
        int cores = boost::thread::hardware_concurrency();
        for (int i = 0; i < cores - 1; i++) {
            threads_.create_thread(boost::bind(&IoService::run,
                                               &io_service_));
        }
    }

    ~ThreadPoolImpl() {
        io_service_.stop();
        threads_.join_all();
    }
};

struct ThreadPool::Impl : public ThreadPoolImpl {
};

ThreadPool::ThreadPool():
    impl_(new Impl) {
}

ThreadPool::~ThreadPool() {
    delete impl_;
    impl_ = 0;
}

static void do_w(ThreadPoolImpl* impl, ThreadWorker* worker) {
    worker->perform();
    {
        Lock lock(impl->finished_mutex_);
        impl->finished_.insert(worker);
    }
    impl->finished_condition_.notify_all();
}

void ThreadPool::post(ThreadWorker* worker) {
    if (impl_->threads_.size() == 0) {
        // one core, no threads
        return;
    }
    impl_->io_service_.post(boost::bind(do_w, impl_, worker));
}

void ThreadPool::wait(ThreadWorker* worker) {
    if (impl_->threads_.size() == 0) {
        // one core, no threads
        return;
    }
    Lock lock(impl_->finished_mutex_);
    WorkersSet& finished = impl_->finished_;
    while (finished.find(worker) == finished.end()) {
        impl_->finished_condition_.wait(lock);
    }
    finished.erase(worker);
}

static ThreadPool* globalInstance_ = 0;

ThreadPool* ThreadPool::globalInstance() {
    // FIXME races and memory leak
    if (!globalInstance_) {
        globalInstance_ = new ThreadPool;
    }
    return globalInstance_;
}

struct GlobalThreadPoolDeleter {
    ~GlobalThreadPoolDeleter() {
        delete globalInstance_;
        globalInstance_ = 0;
    }
} gtpd;

ReusingThreadGroup::ReusingThreadGroup(ThreadPool* pool):
    pool_(pool) {
}

ReusingThreadGroup::~ReusingThreadGroup() {
}

ThreadPool* ReusingThreadGroup::pool() const {
    return pool_ ? : ThreadPool::globalInstance();
}

typedef boost::shared_ptr<ThreadWorker> ThreadWorkerPtr;
typedef std::vector<ThreadWorkerPtr> ThreadWorkers;

struct Waiter {
    const ThreadWorkers& ws_;
    ThreadPool* p_;

    Waiter(const ThreadWorkers& ws, ThreadPool* p):
        ws_(ws), p_(p) {
    }

    ~Waiter() {
        int workers = ws_.size();
        for (int i = 1; i < workers; i++) {
            p_->wait(ws_[i].get());
        }
    }
};

void ReusingThreadGroup::perform_impl() {
    ThreadWorkers ws;
    for (int i = 0; i < workers(); i++) {
        ws.push_back(ThreadWorkerPtr(create_worker()));
    }
    ThreadPool* p = pool();
    for (int i = 1; i < workers(); i++) {
        p->post(ws[i].get());
    }
    ThreadWorker* worker = ws[0].get();
    Waiter waiter(ws, p);
    worker->perform();
    // p->wait(worker) for each worker is called here
    // workers are deleted here
}

}

