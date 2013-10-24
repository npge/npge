/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_THREAD_GROUP_HPP_
#define BR_THREAD_GROUP_HPP_

namespace bloomrepeats {

class Task;
class Worker;
class ThreadGroup;

/** Task to run */
class Task {
public:
    /** Constructor */
    Task(Worker* worker);

    /** Destructor */
    virtual ~Task();

    /** Perform the task */
    void run();

    /** Get Worker instance, created this Task */
    Worker* worker() const;

    /** Get ThreadGroup instance of this task */
    ThreadGroup* thread_group() const;

protected:
    /** Perform the task (implementation) */
    virtual void run_impl() = 0;

private:
    Worker* worker_;
    ThreadGroup* thread_group_;
};

/** Working thread */
class Worker {
public:
    /** Constructor */
    Worker(ThreadGroup* thread_group);

    /** Destructor */
    virtual ~Worker();

    /** Perform tasks */
    void work();

    /** Perform the task */
    void run(Task* task);

    /** Get ThreadGroup instance of this task */
    ThreadGroup* thread_group() const;

protected:
    /** Perform tasks */
    virtual void work_impl();

    /** Perform the task.
    Default implementation calls task().
    */
    virtual void run_impl(Task* task);

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

    /** Perform one worker.
    This is called from each thread by perform_impl().
    */
    void perform_one();

    /* With each call, return new task or empty function.
    Result=0 means "end" of task collection.
    Get mutes and call create_task_impl().
    Caller takes ownership.
    */
    Task* create_task(Worker* worker);

    /** Create new worker.
    Caller takes ownership.
    */
    Worker* create_worker();

protected:
    /* With each call, return new task or empty function.
    Result=0 means "end" of task collection.
    No simultaneous calls of one task generator are allowed.
    Call this under mutex.
    */
    virtual Task* create_task_impl(Worker* worker) = 0;

    /** Create new worker.
    This is called from just created thread.
    Worker's constructor and destructor are thread's
    initializer and finalizer.
    */
    virtual Worker* create_worker_impl();

    /** Perform tasks */
    virtual void perform_impl(int workers);

    /** Perform one worker.
    Create worker, call its work() method, delete worker.
    */
    virtual void perform_one_impl();

private:
    struct Impl;
    Impl* impl_;
};

}

#endif

