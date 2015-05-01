/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
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

typedef std::map<Sequence*, int> Indexes;

static bool isEnded(const BSA& rows, const Indexes& indexes) {
    bool all_ended = true;
    bool any_ended = false;
    BOOST_FOREACH (const BSA::value_type& seq_row, rows) {
        Sequence* seq = seq_row.first;
        const BSRow& row = seq_row.second;
        int index = indexes.find(seq)->second;
        bool ended = index >= row.fragments.size();
        if (ended) {
            any_ended = true;
        } else {
            all_ended = false;
        }
    }
    ASSERT_EQ(all_ended, any_ended);
    return all_ended;
}

static bool isStem(Fragment* f, int genomes) {
    ASSERT_TRUE(f);
    Block* b = f->block();
    ASSERT_TRUE(b);
    return is_exact_stem(b, genomes);
}

static void moveStem(BSA& aln, const BSA& rows,
                     Indexes& indexes, int genomes) {
    BOOST_FOREACH (const BSA::value_type& seq_row, rows) {
        Sequence* seq = seq_row.first;
        const BSRow& row = seq_row.second;
        int& index = indexes[seq];
        ASSERT_LT(index, row.fragments.size());
        Fragment* f = row.fragments[index];
        ASSERT_TRUE(isStem(f, genomes));
        aln[seq].fragments.push_back(f);
        index += 1;
    }
}

static void moveNonStem(BSA& aln, const BSA& rows,
                        Indexes& indexes, int genomes) {
    BSA rows2;
    BOOST_FOREACH (const BSA::value_type& seq_row, rows) {
        Sequence* seq = seq_row.first;
        const BSRow& row = seq_row.second;
        int& index = indexes[seq];
        BSRow& row2 = rows2[seq];
        row2.ori = 1;
        while (!isStem(row.fragments[index], genomes)) {
            row2.fragments.push_back(row.fragments[index]);
            index += 1;
        }
    }
    BSA aln2;
    rows2Bsa(aln2, rows2, genomes);
    BOOST_FOREACH (BSA::value_type& seq_row, aln) {
        Sequence* seq = seq_row.first;
        BSRow& row = seq_row.second;
        const BSRow& row2 = aln2[seq];
        ASSERT_EQ(row2.ori, 1);
        BOOST_FOREACH (Fragment* f, row2.fragments) {
            row.fragments.push_back(f);
        }
    }
}

static void bsaForGlobal(
    BSA& aln, const BSA& rows, int genomes) {
    Indexes indexes;
    BOOST_FOREACH (const BSA::value_type& seq_row, rows) {
        Sequence* seq = seq_row.first;
        const BSRow& row = seq_row.second;
        ASSERT_GTE(row.fragments.size(), 1);
        aln[seq].ori = row.ori;
        indexes[seq] = 0;
    }
    moveStem(aln, rows, indexes, genomes);
    while (!isEnded(rows, indexes)) {
        // find fragments of non-stem blocks
        moveNonStem(aln, rows, indexes, genomes);
        // add fragments of stem-blocks
        moveStem(aln, rows, indexes, genomes);
    }
}

static void makeBsa(
    Block* master_block, BSA& aln, const VectorFc& fc,
    int genomes, char type) {
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
    if (type == 'g') {
        bsaForGlobal(aln, rows, genomes);
    } else {
        rows2Bsa(aln, rows, genomes);
    }
    // check BSA is valid
    BOOST_FOREACH (const BSA::value_type& seq_row, aln) {
        Sequence* seq = seq_row.first;
        const BSRow& row = seq_row.second;
        int ori = row.ori;
        const Fragments& ff = row.fragments;
        Fragment* prev = 0;
        BOOST_FOREACH (Fragment* f, ff) {
            if (prev && f) {
                ASSERT_MSG(fc.neighbor(prev, ori) == f,
                        (prev->id() + " " + f->id()).c_str());
            }
            if (f) {
                prev = f;
            }
        }
    }
}

void LocalBSA::run_impl() const {
    BlockSet& bs = *block_set();
    VectorFc fc;
    fc.add_bs(bs);
    fc.prepare();
    int genomes = genomes_number(*block_set());
    BOOST_FOREACH (Block* master_block, *other()) {
        ASSERT_GTE(master_block->name().size(), 1);
        char type = master_block->name()[0];
        ASSERT_MSG(type == 'g' || type == 'i',
                   "Master block must be either global "
                   "or intermediate");
        ASSERT_FALSE(has_repeats(master_block));
        BSA& aln = block_set()->bsa(master_block->name());
        makeBsa(master_block, aln, fc, genomes, type);
    }
}

const char* LocalBSA::name_impl() const {
    return "Build blockset alignments for global blocks";
}

}

