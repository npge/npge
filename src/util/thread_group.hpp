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

/** Task generator.
With each call, return new task or empty function.
Empty function means "end" of task collection.

No simultaneous calls of one task generator are allowed.
Call this under mutex.
*/
typedef boost::function<Task()> TaskGenerator;

/** Run tasks on thread group */
void do_tasks(TaskGenerator task_generator, int workers);

/** Vector of tasks */
typedef std::vector<Task> Tasks;

/** Create task generator operating on the tasks list */
TaskGenerator tasks_to_generator(Tasks& tasks);

}

#endif

