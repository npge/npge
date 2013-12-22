/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <climits>
#include <vector>
#include <string>
#include <boost/cast.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>

#include "block_hash.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "thread_group.hpp"

namespace bloomrepeats {

uint32_t block_hash(const Block* block) {
    std::vector<std::string> fragment_ids;
    BOOST_FOREACH (Fragment* f, *block) {
        fragment_ids.push_back(f->id());
    }
    std::sort(fragment_ids.begin(), fragment_ids.end());
    std::string joint = boost::algorithm::join(fragment_ids, " ");
    const int LOOP_SIZE = sizeof(uint32_t) * 2; // 2 = for * and for ^
    int new_size = ((joint.size() + LOOP_SIZE - 1) / LOOP_SIZE) * LOOP_SIZE;
    joint.resize(new_size, ' ');
    const uint32_t* value = reinterpret_cast<const uint32_t*>(joint.c_str());
    int loops = joint.size() / LOOP_SIZE;
    uint32_t a = 1;
    for (int i = 0; i < loops; i++) {
        a *= value[2 * i];
        a ^= value[2 * i + 1];
    }
    return a;
}

class HashTask;
class HashWorker;
class HashGroup;

class HashGroup : public ThreadGroup {
public:
    HashGroup(const BlockSet& block_set):
        it_(block_set.begin()), end_(block_set.end()), hash_(0)
    { }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl();

    BlockSet::const_iterator it_;
    BlockSet::const_iterator end_;

    uint32_t hash_;
};

class HashWorker : public ThreadWorker {
public:
    HashWorker(HashGroup* thread_group):
        ThreadWorker(thread_group), hash_(0)
    { }

    ~HashWorker() {
        HashGroup* g = boost::polymorphic_downcast<HashGroup*>(thread_group());
        g->hash_ ^= hash_;
    }

    uint32_t hash_;
};

class HashTask : public ThreadTask {
public:
    HashTask(const Block* block, HashWorker* worker):
        ThreadTask(worker), block_(block)
    { }

    void run_impl() {
        HashWorker* w = boost::polymorphic_downcast<HashWorker*>(worker());
        w->hash_ ^= block_hash(block_);
    }

private:
    const Block* block_;
};

ThreadTask* HashGroup::create_task_impl(ThreadWorker* worker) {
    while (it_ != end_ && (*it_)->size() <= 1) {
        it_++;
    }
    if (it_ == end_) {
        return 0;
    } else {
        HashWorker* w = boost::polymorphic_downcast<HashWorker*>(worker);
        HashTask* task = new HashTask(*it_, w);
        it_++;
        return task;
    }
}

ThreadWorker* HashGroup::create_worker_impl() {
    return new HashWorker(this);
}

uint32_t blockset_hash(const BlockSet& block_set, int workers) {
    HashGroup hash_group((block_set));
    hash_group.perform(workers);
    return hash_group.hash_;
}

}
