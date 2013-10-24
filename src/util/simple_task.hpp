/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SIMPLE_TASK_HPP_
#define BR_SIMPLE_TASK_HPP_

#include <vector>
#include <boost/function.hpp>

namespace bloomrepeats {

typedef boost::function<void()> SimpleTask;

typedef boost::function<SimpleTask()> SimpleTaskGenerator;

/** Run tasks on thread group.
\param task_generator Generator of tasks
\param workers Number of working thread, including main thread
\param thread_init SimpleTask which is run at thread start (if specified)
\param thread_finish SimpleTask which is run after the thread did the work
*/
void do_tasks(SimpleTaskGenerator task_generator, int workers,
              SimpleTask thread_init = SimpleTask(),
              SimpleTask thread_finish = SimpleTask());

/** Vector of tasks */
typedef std::vector<SimpleTask> SimpleTasks;

/** Create task generator operating on the tasks list */
SimpleTaskGenerator tasks_to_generator(SimpleTasks& tasks);

}

#endif

