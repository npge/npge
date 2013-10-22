/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_THREAD_GROUP_HPP_
#define BR_THREAD_GROUP_HPP_

#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>

namespace bloomrepeats {

/** Task to run */
typedef boost::function<void()> Task;

/** Vector of tasks */
typedef std::vector<Task> Tasks;

/** Run tasks on thread group */
void do_tasks(Tasks& tasks, int workers);

}

#endif

