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

static void rows2Bsa(BSA& aln, const BSA& rows, int genomes) {
    boost::scoped_ptr<TreeNode> tree(bsa_make_tree(rows));
    bool try_inverse = false;
    bsa_make_aln_by_tree(aln, rows, tree.get(),
                         genomes, try_inverse);
}

struct FragmentCompareD {
    bool operator()(const Fragment* f1,
                    const Fragment* f2) const {
        return *f1 < *f2;
    }
};

static void makeBsa(
    Block* master_block, BSA& aln, const VectorFc& fc,
    int genomes) {
    BSA rows;
    BOOST_FOREACH (Fragment* f, *master_block) {
        Sequence* seq = f->seq();
        Fragments& ff = rows[seq].fragments;
        fc.find_overlap_fragments(ff, f);
        std::sort(ff.begin(), ff.end(),
                  FragmentCompareD());
        if (f->ori() == -1) {
            std::reverse(ff.begin(), ff.end());
        }
        rows[seq].ori = f->ori();
    }
    rows2Bsa(aln, rows, genomes);
}

void LocalBSA::run_impl() const {
    BlockSet& bs = *block_set();
    VectorFc fc;
    fc.add_bs(bs);
    fc.prepare();
    int genomes = genomes_number(*block_set());
    BOOST_FOREACH (Block* master_block, *other()) {
        ASSERT_FALSE(has_repeats(master_block));
        BSA& aln = block_set()->bsa(master_block->name());
        makeBsa(master_block, aln, fc, genomes);
    }
}

const char* LocalBSA::name_impl() const {
    return "Build blockset alignments for global blocks";
}

}

