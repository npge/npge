/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "simple_task.hpp"
#include "thread_pool.hpp"

namespace npge {

class Task_ : public ThreadTask {
public:
    Task_(Task task, ThreadWorker* worker):
        ThreadTask(worker), task_(task) {
    }

    void run_impl() {
        task_();
    }

private:
    Task task_;
};

class Worker_ : public ThreadWorker {
public:
    Worker_(Task thread_init, Task thread_finish,
            ReusingThreadGroup* thread_group):
        ThreadWorker(thread_group),
        thread_init_(thread_init),
        thread_finish_(thread_finish) {
    }

    void work_impl() {
        if (thread_init_) {
            thread_init_();
        }
        ThreadWorker::work_impl();
        if (thread_finish_) {
            thread_finish_();
        }
    }

private:
    Task thread_init_;
    Task thread_finish_;
};

class ThreadGroup_ : public ReusingThreadGroup {
public:
    ThreadGroup_(TaskGenerator task_generator,
                 Task thread_init, Task thread_finish):
        task_generator_(task_generator),
        thread_init_(thread_init),
        thread_finish_(thread_finish) {
    }

    ThreadTask* create_task_impl(ThreadWorker* worker) {
        Task simple_task = task_generator_();
        return simple_task ? new Task_(simple_task, worker) : 0;
    }

    ThreadWorker* create_worker_impl() {
        return new Worker_(thread_init_, thread_finish_, this);
    }

private:
    TaskGenerator task_generator_;
    Task thread_init_;
    Task thread_finish_;
};

void do_tasks(TaskGenerator task_generator, int workers,
              Task thread_init, Task thread_finish) {
    ThreadGroup_ thread_group(task_generator, thread_init, thread_finish);
    thread_group.set_workers(workers);
    thread_group.perform();
}

class VectorTaskGenerator {
public:
    VectorTaskGenerator(Tasks& tasks):
        tasks_(tasks) {
    }

    Task operator()() {
        Task task;
        if (!tasks_.empty()) {
            task = tasks_.back();
            tasks_.pop_back();
        }
        return task;
    }

private:
    Tasks& tasks_;
};

TaskGenerator tasks_to_generator(Tasks& tasks) {
    VectorTaskGenerator task_generator(tasks);
    return TaskGenerator(task_generator);
}

}

