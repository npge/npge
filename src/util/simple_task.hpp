/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SIMPLE_TASK_HPP_
#define NPGE_SIMPLE_TASK_HPP_

#include <vector>
#include <boost/function.hpp>

namespace npge {

typedef boost::function<void()> Task;

typedef boost::function<Task()> TaskGenerator;

/** Run tasks on thread group.
\param task_generator Generator of tasks
\param workers Number of working thread, including main thread
\param thread_init Task which is run at thread start (if specified)
\param thread_finish Task which is run after the thread did the work
*/
void do_tasks(TaskGenerator task_generator, int workers,
              Task thread_init = Task(),
              Task thread_finish = Task());

/** Vector of tasks */
typedef std::vector<Task> Tasks;

/** Create task generator operating on the tasks list */
TaskGenerator tasks_to_generator(Tasks& tasks);

}

#endif

