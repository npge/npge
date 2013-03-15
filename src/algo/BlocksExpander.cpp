/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <set>
#include <boost/foreach.hpp>

#include "BlocksExpander.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

BlocksExpander::BlocksExpander(int batch):
    ExpanderBase(batch)
{ }

void BlocksExpander::add_options_impl(po::options_description& desc) const {
    ExpanderBase::add_options_impl(desc);
}

void BlocksExpander::apply_options_impl(const po::variables_map& vm) {
    ExpanderBase::apply_options_impl(vm);
}

bool BlocksExpander::expand(Block* block) const {
    bool result = false;
    std::set<Block*> visited;
    std::vector<Fragment*> fragments(block->begin(), block->end());
    BOOST_FOREACH (Fragment* f, fragments) {
        for (int ori = 1; ori >= -1; ori -= 2) {
            Fragment* neighbor = f->neighbor(ori);
            if (neighbor) {
                FragmentDiff diff = neighbor->diff_to(*f);
                Block* b = neighbor->block();
                if (b && b != block && visited.find(b) == visited.end()) {
                    visited.insert(b);
                    BOOST_FOREACH (Fragment* fn, *b) {
                        Fragment candidate;
                        candidate.apply_coords(*fn);
                        candidate.patch(diff);
                        if (candidate.valid() &&
                                !block->common_positions(candidate) &&
                                aligned(*f, candidate)) {
                            Fragment* new_f = new Fragment();
                            new_f->apply_coords(candidate);
                            block->insert(new_f);
                            new_f->find_place(fn);
                            result = true;
                        }
                    }
                }
            }
        }
    }
    return result;
}

bool BlocksExpander::run_impl() const {
    bool result = false;
    BOOST_FOREACH (Block* block, *block_set()) {
        BOOST_ASSERT(block);
        result |= expand(block);
    }
    return result;
}

}

