/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_THREAD_GROUP_HPP_
#define NPGE_THREAD_GROUP_HPP_

#include <string>

namespace npge {

class ThreadTask;
class ThreadWorker;
class ThreadGroup;

/** ThreadTask to run */
class ThreadTask {
public:
    /** Constructor */
    ThreadTask(ThreadWorker* worker);

    /** Destructor */
    virtual ~ThreadTask();

    /** Perform the task */
    void run();

    /** Get ThreadWorker instance, created this ThreadTask */
    ThreadWorker* worker() const;

    /** Get ThreadGroup instance of this task */
    ThreadGroup* thread_group() const;

protected:
    /** Perform the task (implementation) */
    virtual void run_impl() = 0;

private:
    ThreadWorker* worker_;
    ThreadGroup* thread_group_;
};

/** Working thread */
class ThreadWorker {
public:
    /** Constructor */
    ThreadWorker(ThreadGroup* thread_group);

    /** Destructor.
    Calls check_worker.
    */
    virtual ~ThreadWorker();

    /** Perform.
    This methos is called from new thread.
    */
    void perform();

    /** Perform tasks */
    void work();

    /** Perform the task */
    void run(ThreadTask* task);

    /** Get ThreadGroup instance of this task */
    ThreadGroup* thread_group() const;

    /** Get error message */
    const std::string& error_message() const;

protected:
    /** Run work() under try-catch if workers() >= 2 */
    virtual void perform_impl();

    /** Perform tasks.
    Reimplement this to set thread's initializer and finalizer,
    which are called from new thread.
    */
    virtual void work_impl();

    /** Perform the task.
    Default implementation calls task().
    */
    virtual void run_impl(ThreadTask* task);

private:
    ThreadGroup* thread_group_;
    std::string error_message_;
};

/** Main class for running work */
class ThreadGroup {
public:
    /** Constructor */
    ThreadGroup();

    /** Destructor */
    virtual ~ThreadGroup();

    /** Perform tasks */
    void perform();

    /* With each call, return new task or empty function.
    Result=0 means "end" of task collection.
    Get mutes and call create_task_impl().
    Caller takes ownership.
    Calls check_worker.
    */
    ThreadTask* create_task(ThreadWorker* worker);

    /** Create new worker.
    Caller takes ownership.
    */
    ThreadWorker* create_worker();

    /** Check the worker for errors.
    This method is called from create_task() and ThreadWorker's
    destructor.
    */
    void check_worker(ThreadWorker* worker);

    /** Set number of workers */
    void set_workers(int workers);

    /** Return number of workers.
    Returns number of CPUs if workers == -1.
    */
    int workers() const;

protected:
    /* With each call, return new task or empty function.
    Result=0 means "end" of task collection.
    No simultaneous calls of one task generator are allowed.
    Call this under mutex.
    */
    virtual ThreadTask* create_task_impl(ThreadWorker* worker) = 0;

    /** Create new worker.
    This is called from just created thread.
    ThreadWorker's constructor and destructor are called from main thread.
    */
    virtual ThreadWorker* create_worker_impl();

    /** Check the worker for errors.
    If worker->error_message() is not empty, stores it
    so that create_task() returns 0 and perform_impl()
    throws Exception with error message afterwards.
    */
    virtual void check_worker_impl(ThreadWorker* worker);

    /** Perform tasks */
    virtual void perform_impl();

private:
    struct Impl;
    Impl* impl_;
};

}

#endif

