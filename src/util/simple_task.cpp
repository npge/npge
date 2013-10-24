/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "simple_task.hpp"
#include "thread_group.hpp"

namespace bloomrepeats {

class SimpleTask_ : public Task {
public:
    SimpleTask_(SimpleTask task, Worker* worker):
        Task(worker), task_(task)
    { }

    void run_impl() {
        task_();
    }

private:
    SimpleTask task_;
};

class SimpleWorker_ : public Worker {
public:
    SimpleWorker_(SimpleTask thread_init, SimpleTask thread_finish,
                  ThreadGroup* thread_group):
        Worker(thread_group),
        thread_finish_(thread_finish) {
        if (thread_init) {
            thread_init();
        }
    }

    ~SimpleWorker_() {
        if (thread_finish_) {
            thread_finish_();
        }
    }

private:
    SimpleTask thread_finish_;
};

class SimpleThreadGroup_ : public ThreadGroup {
public:
    SimpleThreadGroup_(SimpleTaskGenerator task_generator,
                       SimpleTask thread_init, SimpleTask thread_finish):
        task_generator_(task_generator),
        thread_init_(thread_init),
        thread_finish_(thread_finish)
    { }

    Task* create_task_impl(Worker* worker) {
        SimpleTask simple_task = task_generator_();
        return simple_task ? new SimpleTask_(simple_task, worker) : 0;
    }

    Worker* create_worker_impl() {
        return new SimpleWorker_(thread_init_, thread_finish_, this);
    }

private:
    SimpleTaskGenerator task_generator_;
    SimpleTask thread_init_;
    SimpleTask thread_finish_;
};

void do_tasks(SimpleTaskGenerator task_generator, int workers,
              SimpleTask thread_init, SimpleTask thread_finish) {
    SimpleThreadGroup_ thread_group(task_generator, thread_init, thread_finish);
    thread_group.perform(workers);
}

class VectorTaskGenerator {
public:
    VectorTaskGenerator(SimpleTasks& tasks):
        tasks_(tasks)
    { }

    SimpleTask operator()() {
        SimpleTask task;
        if (!tasks_.empty()) {
            task = tasks_.back();
            tasks_.pop_back();
        }
        return task;
    }

private:
    SimpleTasks& tasks_;
};

SimpleTaskGenerator tasks_to_generator(SimpleTasks& tasks) {
    VectorTaskGenerator task_generator(tasks);
    return SimpleTaskGenerator(task_generator);
}

}

