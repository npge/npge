/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <algorithm>
#include <ostream>
#include <boost/foreach.hpp>

#include "block_set_alignment.hpp"
#include "GeneralAligner.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

BSRow::BSRow():
    ori(1)
{ }

int bsa_length(const BSA& bsa) {
    if (bsa.empty()) {
        return 0;
    } else {
        int length = bsa.begin()->second.fragments.size();
        BOOST_FOREACH (const BSA::value_type& seq_and_row, bsa) {
            const BSRow& row = seq_and_row.second;
            BOOST_ASSERT(row.fragments.size() == length);
        }
        return length;
    }
}

static struct FragmentCompareA {
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        return *f1 < *f2;
    }
} fragment_compare_a;

void bsa_make_rows(BSA& rows, const BlockSet& bs,
                   const std::string& chr) {
    if (rows.empty()) {
        BOOST_FOREACH (const SequencePtr& seq, bs.seqs()) {
            if (chr.empty() || seq->chromosome() == chr) {
                rows[seq.get()] = BSRow();
            }
        }
    }
    BOOST_FOREACH (const Block* block, bs) {
        BOOST_FOREACH (Fragment* fragment, *block) {
            Sequence* seq = fragment->seq();
            BOOST_ASSERT(seq);
            BSA::iterator it = rows.find(seq);
            if (it != rows.end()) {
                BSRow& row = it->second;
                row.fragments.push_back(fragment);
            }
        }
    }
    BOOST_FOREACH (BSA::value_type& seq_and_row, rows) {
        BSRow& row = seq_and_row.second;
        Fragments& fragments = row.fragments;
        std::sort(fragments.begin(), fragments.end(),
                  fragment_compare_a);
    }
}

void bsa_inverse(BSA& aln) {
    BOOST_FOREACH (BSA::value_type& seq_and_row, aln) {
        BSRow& row = seq_and_row.second;
        row.ori *= (-1);
        Fragments& fragments = row.fragments;
        std::reverse(fragments.begin(), fragments.end());
    }
}

struct BSContents {
    const BSA* first_;
    const BSA* second_;

    BSContents():
        first_(0), second_(0)
    { }

    BSContents(const BSA& first, const BSA& second):
        first_(&first), second_(&second)
    { }

    int first_size() const {
        return bsa_length(*first_);
    }

    int second_size() const {
        return bsa_length(*second_);
    }

    int substitution(int row, int col) const {
        typedef std::pair<Block*, int> BlockOri;
        typedef std::set<BlockOri> BlockOriSet;
        BlockOriSet bos;
        BOOST_FOREACH (const BSA::value_type& seq_and_row, *first_) {
            const BSRow& bs_row = seq_and_row.second;
            Fragment* fragment = bs_row.fragments[row];
            if (fragment) {
                Block* block = fragment->block();
                int ori = bs_row.ori * fragment->ori();
                bos.insert(BlockOri(block, ori));
            }
        }
        BOOST_FOREACH (const BSA::value_type& seq_and_row, *second_) {
            const BSRow& bs_row = seq_and_row.second;
            Fragment* fragment = bs_row.fragments[col];
            if (fragment) {
                Block* block = fragment->block();
                int ori = bs_row.ori * fragment->ori();
                if (bos.find(BlockOri(block, ori)) != bos.end()) {
                    return 0;
                }
            }
        }
        return 1;
    }
};

typedef GeneralAligner<BSContents>::PairAlignment PairAlignment;

void bsa_align(BSA& both, int& score,
               const BSA& first, const BSA& second) {
    GeneralAligner<BSContents> ga;
    ga.set_max_errors(-1); // unlimited errors
    int first_size = bsa_length(first);
    int second_size = bsa_length(second);
    ga.set_gap_range(std::max(first_size, second_size));
    // ^^ FIXME GeneralAligner full matrix
    ga.set_contents(BSContents(first, second));
    int first_last, second_last;
    ga.align(first_last, second_last);
    first_last = first_size - 1;
    second_last = second_size - 1;
    score = ga.at(first_last, second_last);
    PairAlignment alignment;
    ga.export_alignment(first_last, second_last, alignment);
    typedef std::pair<int, int> Match;
    both.clear();
    std::vector<const BSA*> bsas;
    bsas.push_back(&first);
    bsas.push_back(&second);
    BOOST_FOREACH (const BSA* orig_aln, bsas) {
        bool is_first = (orig_aln == &first);
        BOOST_FOREACH (const BSA::value_type& seq_and_row, *orig_aln) {
            Sequence* seq = seq_and_row.first;
            const BSRow& orig_row = seq_and_row.second;
            const Fragments& orig_f = orig_row.fragments;
            BSRow& new_row = both[seq];
            new_row.ori = orig_row.ori;
            Fragments& new_f = new_row.fragments;
            BOOST_FOREACH (const Match& m, alignment) {
                int in_orig = is_first ? m.first : m.second;
                new_f.push_back((in_orig != -1) ? orig_f[in_orig] : 0);
            }
        }
    }
}

void bsa_make_aln(BSA& aln, const BSA& rows) {
    aln.clear();
    if (rows.empty()) {
        return;
    }
    BSA::const_iterator first_it = rows.begin();
    Sequence* first_seq = first_it->first;
    aln[first_seq] = first_it->second;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, rows) {
        Sequence* seq = seq_and_row.first;
        const BSRow& row = seq_and_row.second;
        if (seq != first_seq) {
            BSA both_direct, both_inverse;
            int score_direct, score_inverse;
            {
                BSA second;
                second[seq] = row;
                bsa_align(both_direct, score_direct, aln, second);
            }
            {
                BSA second;
                second[seq] = row;
                bsa_inverse(second);
                bsa_align(both_inverse, score_inverse, aln, second);
            }
            bool use_direct = (score_direct < score_inverse);
            BSA& both = use_direct ? both_direct : both_inverse;
            aln.swap(both);
        }
    }
}

void bsa_print(std::ostream& out, const BSA& aln, bool blocks) {
    BOOST_FOREACH (const BSA::value_type& seq_and_row, aln) {
        Sequence* seq = seq_and_row.first;
        const BSRow& row = seq_and_row.second;
        out << ((row.ori == 1) ? '+' : '-');
        out << seq->genome();
        BOOST_FOREACH (Fragment* fragment, row.fragments) {
            out << '\t';
            if (fragment) {
                if (blocks) {
                    out << fragment->block()->name();
                } else {
                    out << fragment->id();
                }
            } else {
                out << '-';
            }
        }
        out << "\n";
    }
}

}

