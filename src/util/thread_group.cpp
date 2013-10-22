/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "thread_group.hpp"

namespace bloomrepeats {

static void process_some_seqs(Tasks& tasks, boost::mutex* mutex) {
    while (true) {
        Task task;
        if (mutex) {
            mutex->lock();
        }
        if (!tasks.empty()) {
            task = tasks.back();
            tasks.pop_back();
        }
        if (mutex) {
            mutex->unlock();
        }
        if (!task.empty()) {
            task();
        } else {
            break;
        }
    }
}

void do_tasks(Tasks& tasks, int workers) {
    boost::mutex mutex;
    boost::thread_group threads;
    for (int i = 1; i < workers; i++) {
        threads.create_thread(boost::bind(process_some_seqs, boost::ref(tasks),
                                          &mutex));
    }
    process_some_seqs(tasks, &mutex);
    threads.join_all();
}

}

