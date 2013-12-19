/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_THREAD_GROUP_HPP_
#define BR_THREAD_GROUP_HPP_

namespace bloomrepeats {

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

    /** Destructor */
    virtual ~ThreadWorker();

    /** Perform tasks */
    void work();

    /** Perform the task */
    void run(ThreadTask* task);

    /** Get ThreadGroup instance of this task */
    ThreadGroup* thread_group() const;

protected:
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
};

/** Main class for running work */
class ThreadGroup {
public:
    /** Constructor */
    ThreadGroup();

    /** Destructor */
    virtual ~ThreadGroup();

    /** Perform tasks */
    void perform(int workers);

    /* With each call, return new task or empty function.
    Result=0 means "end" of task collection.
    Get mutes and call create_task_impl().
    Caller takes ownership.
    */
    ThreadTask* create_task(ThreadWorker* worker);

    /** Create new worker.
    Caller takes ownership.
    */
    ThreadWorker* create_worker();

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

    /** Perform tasks */
    virtual void perform_impl(int workers);

private:
    struct Impl;
    Impl* impl_;
};

}

#endif

