/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "simple_task.hpp"
#include "thread_group.hpp"

namespace bloomrepeats {

class Task_ : public ThreadTask {
public:
    Task_(Task task, ThreadWorker* worker):
        ThreadTask(worker), task_(task)
    { }

    void run_impl() {
        task_();
    }

private:
    Task task_;
};

class Worker_ : public ThreadWorker {
public:
    Worker_(Task thread_init, Task thread_finish,
            ThreadGroup* thread_group):
        ThreadWorker(thread_group),
        thread_finish_(thread_finish) {
        if (thread_init) {
            thread_init();
        }
    }

    ~Worker_() {
        if (thread_finish_) {
            thread_finish_();
        }
    }

private:
    Task thread_finish_;
};

class ThreadGroup_ : public ThreadGroup {
public:
    ThreadGroup_(TaskGenerator task_generator,
                 Task thread_init, Task thread_finish):
        task_generator_(task_generator),
        thread_init_(thread_init),
        thread_finish_(thread_finish)
    { }

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
    thread_group.perform(workers);
}

class VectorTaskGenerator {
public:
    VectorTaskGenerator(Tasks& tasks):
        tasks_(tasks)
    { }

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

