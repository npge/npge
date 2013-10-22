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

/** Run tasks on thread group.
\param task_generator Generator of tasks
\param workers Number of working thread, including main thread
\param thread_init Task which is run at thread start (if specified)
\param thread_finish Task which is run after the thread did the work
*/
void do_tasks(TaskGenerator task_generator, int workers,
              Task thread_init = Task(), Task thread_finish = Task());

/** Vector of tasks */
typedef std::vector<Task> Tasks;

/** Create task generator operating on the tasks list */
TaskGenerator tasks_to_generator(Tasks& tasks);

}

#endif

