/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <set>
#include <boost/foreach.hpp>

#include "Pipe.hpp"
#include "block_hash.hpp"

namespace npge {

struct Pipe::Impl {
    std::vector<Processor*> processors_;
    int max_iterations_;
    bool stopped_;

    Impl():
        max_iterations_(1),
        stopped_(false) {
    }
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

int Pipe::max_iterations() const {
    return impl_->max_iterations_;
}

void Pipe::set_max_iterations(int max_iterations) {
    impl_->max_iterations_ = max_iterations;
}

std::vector<Processor*> Pipe::processors() const {
    return impl_->processors_;
}

void Pipe::stop() {
    impl_->stopped_ = true;
}

void Pipe::run_impl() const {
    BOOST_FOREACH (Processor* processor, impl_->processors_) {
        processor->set_workers(workers());
    }
    std::set<hash_t> hashes;
    hashes.insert(blockset_hash(*block_set(), workers()));
    impl_->stopped_ = false;
    for (int i = 0; i < max_iterations() || max_iterations() == -1; i++) {
        BOOST_FOREACH (Processor* processor, impl_->processors_) {
            processor->run();
        }
        hash_t new_hash = blockset_hash(*block_set(), workers());
        if (hashes.find(new_hash) != hashes.end()) {
            break;
        }
        if (impl_->stopped_) {
            break;
        }
        hashes.insert(new_hash);
    }
}

}

