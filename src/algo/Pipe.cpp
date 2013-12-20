/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>

#include "Pipe.hpp"
#include "block_hash.hpp"

namespace bloomrepeats {

struct Pipe::Impl {
    std::vector<Processor*> processors_;
    int max_loops_;

    Impl():
        max_loops_(1)
    { }
};

Pipe::Pipe(BlockSetPtr other) {
    impl_ = new Impl;
    set_other(other);
}

Pipe::~Pipe() {
    delete impl_;
}

Pipe& Pipe::add(Processor* processor, const std::string& options) {
    processor->set_options(options, this);
    impl_->processors_.push_back(processor);
    processor->set_parent(this);
    return *this;
}

int Pipe::max_loops() const {
    return impl_->max_loops_;
}

void Pipe::set_max_loops(int max_loops) {
    impl_->max_loops_ = max_loops;
}

std::vector<Processor*> Pipe::processors() const {
    return impl_->processors_;
}

bool Pipe::run_impl() const {
    BOOST_FOREACH (Processor* processor, impl_->processors_) {
        processor->set_workers(workers());
    }
    uint32_t hash = blockset_hash(*block_set(), workers());
    for (int i = 0; i < max_loops() || max_loops() == -1; i++) {
        BOOST_FOREACH (Processor* processor, impl_->processors_) {
            processor->run();
        }
        uint32_t new_hash = blockset_hash(*block_set(), workers());
        if (new_hash == hash) {
            break;
        }
        hash = new_hash;
    }
    return true;
}

}

