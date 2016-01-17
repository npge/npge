/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_THREAD_POOL_HPP_
#define NPGE_THREAD_POOL_HPP_

#include "thread_group.hpp"

namespace npge {

/** Thread pool.
Number of threads is set to number of cores.
*/
class ThreadPool {
public:
    /** Constructor */
    ThreadPool();

    /** Destructor */
    ~ThreadPool();

    /** Run worker in thread */
    void post(ThreadWorker* worker);

    /** Block current thread until this worker has finished */
    void wait(ThreadWorker* worker);

    /** Return global thread pool */
    static ThreadPool* globalInstance();

private:
    class Impl;
    Impl* impl_;
};

/** Main class for running work (can reuse threads).
This subclass can use same thread multiple times
for new instances of ReusingThreadGroup.
*/
class ReusingThreadGroup : public ThreadGroup {
public:
    /** Constructor.
    If pool = 0, then the global instance is used.
    */
    ReusingThreadGroup(ThreadPool* pool = 0);

    /** Destructor */
    virtual ~ReusingThreadGroup();

    /** Return current thread pool */
    ThreadPool* pool() const;

protected:
    virtual void perform_impl();

private:
    ThreadPool* pool_;
};

}

#endif

