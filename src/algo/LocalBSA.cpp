/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>

#include "LocalBSA.hpp"
#include "FindBSA.hpp"
#include "FragmentCollection.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "throw_assert.hpp"
#include "bsa_algo.hpp"
#include "tree.hpp"
#include "block_hash.hpp"
#include "global.hpp"

namespace npge {

LocalBSA::LocalBSA() {
    declare_bs("target", "Target blockset");
    declare_bs("other", "Blockset with global blocks");
}

void LocalBSA::run_impl() const {
    BlockSet& bs = *block_set();
    VectorFc fc;
    fc.add_bs(bs);
    fc.prepare();
    int genomes = genomes_number(*block_set());
    BOOST_FOREACH (Block* global_block, *other()) {
        Fragments ff;
        BOOST_FOREACH (Fragment* f, *global_block) {
            fc.find_overlap_fragments(ff, f);
        }
        BSA rows;
        BOOST_FOREACH (Fragment* f, ff) {
            rows[f->seq()].fragments.push_back(f);
        }
        bsa_sort(rows);
        boost::scoped_ptr<TreeNode> tree((bsa_make_tree(rows)));
        BSA& aln = block_set()->bsa(global_block->name());
        bsa_make_aln_by_tree(aln, rows, tree.get(), genomes);
        bsa_orient(aln);
    }
}

const char* LocalBSA::name_impl() const {
    return "Build blockset alignments for global blocks";
}

}

