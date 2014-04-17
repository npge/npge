/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <set>
#include <algorithm>
#include <ostream>
#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "bsa_algo.hpp"
#include "GeneralAligner.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "tree.hpp"
#include "block_hash.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace bloomrepeats {

typedef std::pair<Block*, int> BlockOri;
typedef std::set<BlockOri> BlockOriSet;

static struct FragmentCompareA {
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        return *f1 < *f2;
    }
} fragment_compare_a;

void bsa_make_rows(BSA& rows, const BlockSet& bs) {
    if (rows.empty()) {
        BOOST_FOREACH (const SequencePtr& seq, bs.seqs()) {
            rows[seq.get()] = BSRow();
        }
    }
    BOOST_FOREACH (const Block* block, bs) {
        BOOST_FOREACH (Fragment* fragment, *block) {
            Sequence* seq = fragment->seq();
            ASSERT_TRUE(seq);
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
    int genomes_;

    BSContents():
        first_(0), second_(0), genomes_(0)
    { }

    BSContents(const BSA& first, const BSA& second, int genomes):
        first_(&first), second_(&second), genomes_(genomes)
    { }

    int first_size() const {
        return bsa_length(*first_);
    }

    int second_size() const {
        return bsa_length(*second_);
    }

    bool is_stem(const Block* block) const {
        return is_exact_stem(block, genomes_);
    }

    int substitution(int row, int col) const {
        BlockOriSet bos;
        BOOST_FOREACH (const BSA::value_type& seq_and_row, *first_) {
            const BSRow& bs_row = seq_and_row.second;
            ASSERT_LT(row, bs_row.fragments.size());
            Fragment* fragment = bs_row.fragments[row];
            if (fragment) {
                Block* block = fragment->block();
                int ori = bs_row.ori * fragment->ori();
                bos.insert(BlockOri(block, ori));
            }
        }
        BOOST_FOREACH (const BSA::value_type& seq_and_row, *second_) {
            const BSRow& bs_row = seq_and_row.second;
            ASSERT_LT(col, bs_row.fragments.size());
            Fragment* fragment = bs_row.fragments[col];
            if (fragment) {
                Block* block = fragment->block();
                int ori = bs_row.ori * fragment->ori();
                if (bos.find(BlockOri(block, ori)) != bos.end()) {
                    int score = 2 * log(block->alignment_length());
                    if (is_stem(block)) {
                        score *= 2;
                    }
                    return -1 - score;
                }
            }
        }
        return 1;
    }
};

void bsa_align(BSA& both, int& score,
               const BSA& first, const BSA& second, int genomes) {
    int gap_penalty = 5;
    BSContents bsc((first), second, genomes);
    typedef ContentsProxy<BSContents> BSProxy;
    BSProxy proxy((bsc));
    PairAlignment alignment;
    bool c = bsa_is_circular(first) && bsa_is_circular(second);
    score = find_aln(alignment, proxy, gap_penalty, c);
    typedef std::pair<int, int> Match;
    both.clear();
    std::vector<const BSA*> bsas;
    bsas.push_back(&first);
    bsas.push_back(&second);
    BOOST_FOREACH (const BSA* bsa, bsas) {
        bool is_first = (bsa == &first);
        BOOST_FOREACH (const BSA::value_type& seq_and_row, *bsa) {
            Sequence* seq = seq_and_row.first;
            const BSRow& orig_row = seq_and_row.second;
            const Fragments& orig_f = orig_row.fragments;
            BSRow& new_row = both[seq];
            new_row.ori = orig_row.ori;
            Fragments& new_f = new_row.fragments;
            BOOST_FOREACH (const Match& m, alignment) {
                int in_orig = is_first ? m.first : m.second;
                if (in_orig == -1) {
                    new_f.push_back(0);
                } else {
                    ASSERT_LT(in_orig, orig_f.size());
                    new_f.push_back(orig_f[in_orig]);
                }
            }
        }
    }
}

void bsa_make_aln(BSA& aln, const BSAs& parts, int genomes) {
    aln.clear();
    if (parts.empty()) {
        return;
    }
    aln = parts[0];
    for (int i = 1; i < parts.size(); i++) {
        BSA both_direct, both_inverse;
        int score_direct, score_inverse;
        {
            const BSA& second = parts[i];
            bsa_align(both_direct, score_direct, aln,
                      second, genomes);
        }
        {
            BSA second = parts[i];
            bsa_inverse(second);
            bsa_align(both_inverse, score_inverse, aln,
                      second, genomes);
        }
        bool use_direct = (score_direct < score_inverse);
        BSA& both = use_direct ? both_direct : both_inverse;
        aln.swap(both);
        bsa_move_fragments(aln);
        bsa_remove_pure_gaps(aln);
    }
}

void bsa_make_aln(BSA& aln, const BSA& rows, int genomes) {
    BSAs parts;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, rows) {
        Sequence* seq = seq_and_row.first;
        const BSRow& row = seq_and_row.second;
        parts.push_back(BSA());
        parts.back()[seq] = row;
    }
    bsa_make_aln(aln, parts, genomes);
}

class SequenceLeaf : public LeafNode {
public:
    SequenceLeaf(Sequence* seq, const BSRow* bsrow):
        seq_(seq), bsrow_(bsrow)
    { }

    Sequence* seq() const {
        return seq_;
    }

    const BSRow* bsrow() const {
        return bsrow_;
    }

protected:
    TreeNode* clone_impl() const {
        return new SequenceLeaf(seq_, bsrow_);
    }

    double distance_to_impl(const LeafNode* leaf) const {
        const SequenceLeaf* seq_leaf;
        seq_leaf = dynamic_cast<const SequenceLeaf*>(leaf);
        if (!seq_leaf) {
            throw Exception("Bad leaf type");
        }
        BlockOriSet bos;
        BOOST_FOREACH (Fragment* f, bsrow()->fragments) {
            bos.insert(BlockOri(f->block(), f->ori() * bsrow()->ori));
        }
        int in_both = 0;
        BOOST_FOREACH (Fragment* f, seq_leaf->bsrow()->fragments) {
            BlockOri bo(f->block(), f->ori() * bsrow()->ori);
            if (bos.find(bo) != bos.end()) {
                in_both += 1;
            }
        }
        int this_size = bsrow()->fragments.size();
        int other_size = seq_leaf->bsrow()->fragments.size();
        int total = this_size + other_size - in_both + 1;
        // +1 not to divide by 0
        return double(in_both) / double(total);
    }

    std::string name_impl() const {
        return seq_->name();
    }

private:
    Sequence* seq_;
    const BSRow* bsrow_;
};

static void bsa_make_aln_by_tree(BSA& aln, const TreeNode* tree,
                                 int genomes) {
    const SequenceLeaf* seq_leaf;
    seq_leaf = dynamic_cast<const SequenceLeaf*>(tree);
    if (seq_leaf) {
        aln[seq_leaf->seq()] = *seq_leaf->bsrow();
    } else {
        BSAs parts;
        BOOST_FOREACH (TreeNode* child, tree->children()) {
            parts.push_back(BSA());
            bsa_make_aln_by_tree(parts.back(), child, genomes);
        }
        bsa_make_aln(aln, parts, genomes);
    }
}

void bsa_make_aln_by_tree(BSA& aln, const BSA& rows,
                          const TreeNode* tree0, int genomes) {
    boost::scoped_ptr<TreeNode> tree((bsa_convert_tree(rows, tree0)));
    bsa_make_aln_by_tree(aln, tree.get(), genomes);
}

void bsa_remove_pure_gaps(BSA& aln) {
    std::vector<Fragments*> fragments_v;
    BOOST_FOREACH (BSA::value_type& seq_and_row, aln) {
        BSRow& row = seq_and_row.second;
        fragments_v.push_back(&row.fragments);
    }
    int length = bsa_length(aln);
    int dst = 0;
    for (int src = 0; src < length; src++) {
        bool pure_gap = true;
        BOOST_FOREACH (Fragments* fragments, fragments_v) {
            if ((*fragments)[src]) {
                pure_gap = false;
                break;
            }
        }
        if (!pure_gap) {
            if (src != dst) {
                BOOST_FOREACH (Fragments* fragments, fragments_v) {
                    (*fragments)[dst] = (*fragments)[src];
                }
            }
            dst += 1;
        }
    }
    BOOST_FOREACH (Fragments* fragments, fragments_v) {
        fragments->resize(dst);
    }
}

typedef std::vector<BSRow*> BSRows;

static double count_block_ori(const BSRows& bsrows, int col,
                              Block* block, int ori) {
    double result = 0.0;
    BOOST_FOREACH (const BSRow* bsrow, bsrows) {
        Fragment* f = bsrow->fragments[col];
        if (f && f->block() == block) {
            if (f->ori() * bsrow->ori == ori) {
                result += 1.0;
            } else {
                result += 0.5;
            }
        }
    }
    return result;
}

static bool move_f(const BSRows& bsrows, BSRow& bsrow, int col) {
    Fragment* f = bsrow.fragments[col];
    if (!f) {
        return false;
    }
    Block* block = f->block();
    int ori = f->ori() * bsrow.ori;
    double score = count_block_ori(bsrows, col, block, ori) - 1;
    // -1 because f itself was counted
    ASSERT_GTE(score, 0);
    int best_score = score;
    int best_col = col;
    for (int i = col - 1; i >= 0; i--) {
        if (bsrow.fragments[i]) {
            break;
        }
        double ascore = count_block_ori(bsrows, i, block, ori);
        if (ascore > best_score) {
            best_col = i;
            best_score = ascore;
        }
    }
    for (int i = col + 1; i < bsrow.fragments.size(); i++) {
        if (bsrow.fragments[i]) {
            break;
        }
        double ascore = count_block_ori(bsrows, i, block, ori);
        if (ascore > best_score) {
            best_col = i;
            best_score = ascore;
        }
    }
    if (best_col != col) {
        bsrow.fragments[best_col] = f;
        bsrow.fragments[col] = 0;
        return true;
    }
    return false;
}

void bsa_move_fragments(BSA& aln) {
    BSRows bsrows;
    BOOST_FOREACH (BSA::value_type& seq_and_row, aln) {
        BSRow& row = seq_and_row.second;
        bsrows.push_back(&row);
    }
    int length = bsa_length(aln);
    bool goon = true;
    while (goon) {
        goon = false;
        BOOST_FOREACH (BSRow* bsrow, bsrows) {
            for (int col = 0; col < length; col++) {
                goon |= move_f(bsrows, *bsrow, col);
            }
        }
    }
}

static void append_col(BSRows& new_bsrows, const BSRows& bsrows,
                       int col, int size) {
    for (int i = 0; i < size; i++) {
        BSRow* bsrow = bsrows[i];
        Fragment* fragment = bsrow->fragments[col];
        BSRow* new_bsrow = new_bsrows[i];
        new_bsrow->fragments.push_back(fragment);
    }
}

void bsa_unwind(BSA& aln) {
    BSA new_aln;
    BSRows bsrows, new_bsrows;
    BOOST_FOREACH (BSA::value_type& seq_and_row, aln) {
        Sequence* seq = seq_and_row.first;
        BSRow& row = seq_and_row.second;
        bsrows.push_back(&row);
        BSRow& new_row = new_aln[seq];
        new_row.ori = row.ori;
        new_bsrows.push_back(&new_row);
    }
    int length = bsa_length(aln);
    int size = bsrows.size();
    ASSERT_EQ(size, bsrows.size());
    ASSERT_EQ(size, new_bsrows.size());
    ASSERT_EQ(size, aln.size());
    ASSERT_EQ(size, new_aln.size());
    for (int col = 0; col < length; col++) {
        BlockOriSet bos;
        bool gap = false;
        BOOST_FOREACH (BSRow* bsrow, bsrows) {
            Fragment* fragment = bsrow->fragments[col];
            if (fragment) {
                Block* block = fragment->block();
                BOOST_ASSERT(block);
                int ori = fragment->ori() * bsrow->ori;
                bos.insert(BlockOri(block, ori));
            } else {
                gap = true;
            }
        }
        if (!gap || bos.size() <= 1) {
            // no changes
            append_col(new_bsrows, bsrows, col, size);
        } else {
            // split
            BOOST_FOREACH (const BlockOri& bo, bos) {
                for (int i = 0; i < size; i++) {
                    BSRow* bsrow = bsrows[i];
                    Fragment* fragment = bsrow->fragments[col];
                    if (fragment) {
                        Block* block = fragment->block();
                        BOOST_ASSERT(block);
                        int ori = fragment->ori() * bsrow->ori;
                        if (BlockOri(block, ori) != bo) {
                            fragment = 0;
                        }
                    }
                    BSRow* new_bsrow = new_bsrows[i];
                    new_bsrow->fragments.push_back(fragment);
                }
            }
        }
    }
    aln.swap(new_aln);
}

static void apply_shift(BSA& bsa, int shift) {
    int L = bsa_length(bsa);
    ASSERT_GTE(shift, 0);
    ASSERT_LT(shift, L);
    BOOST_FOREACH (BSA::value_type& seq_and_row, bsa) {
        BSRow& old_bsrow = seq_and_row.second;
        Fragments& old_ff = old_bsrow.fragments;
        Fragments new_ff;
        for (int i = shift; i < L; i++) {
            new_ff.push_back(old_ff[i]);
        }
        for (int i = 0; i < shift; i++) {
            new_ff.push_back(old_ff[i]);
        }
        old_ff.swap(new_ff);
    }
}

static int column_size(int col, const BSRows& bsrows) {
    int occupied = 0;
    BOOST_FOREACH (const BSRow* bsrow, bsrows) {
        Fragment* fragment = bsrow->fragments[col];
        if (fragment) {
            occupied += 1;
        }
    }
    return occupied;
}

typedef std::set<Sequence*> SequenceSet;

static void add_seqs_to_set(SequenceSet& seq_set, int col,
                            const BSRows& bsrows) {
    BOOST_FOREACH (const BSRow* bsrow, bsrows) {
        Fragment* fragment = bsrow->fragments[col];
        if (fragment) {
            Sequence* seq = fragment->seq();
            BOOST_ASSERT(seq);
            seq_set.insert(seq);
        }
    }
}

static bool has_seq(int col, const SequenceSet& seq_set,
                    const BSRows& bsrows) {
    BOOST_FOREACH (const BSRow* bsrow, bsrows) {
        Fragment* fragment = bsrow->fragments[col];
        if (fragment) {
            Sequence* seq = fragment->seq();
            BOOST_ASSERT(seq);
            if (seq_set.find(seq) != seq_set.end()) {
                return true;
            }
        }
    }
    return false;
}

void bsa_move_columns(BSA& aln) {
    BSA new_aln;
    BSRows bsrows, new_bsrows;
    BOOST_FOREACH (BSA::value_type& seq_and_row, aln) {
        Sequence* seq = seq_and_row.first;
        BSRow& row = seq_and_row.second;
        bsrows.push_back(&row);
        BSRow& new_row = new_aln[seq];
        new_row.ori = row.ori;
        new_bsrows.push_back(&new_row);
    }
    bool circular = bsa_is_circular(aln);
    int length = bsa_length(aln);
    int size = bsrows.size();
    ASSERT_EQ(size, bsrows.size());
    ASSERT_EQ(size, new_bsrows.size());
    ASSERT_EQ(size, aln.size());
    ASSERT_EQ(size, new_aln.size());
    for (int col = 0; col < length; col++) {
        if (column_size(col, bsrows) == size) {
            apply_shift(aln, col);
            break;
        }
    }
    typedef std::set<int> IntSet;
    // sorted
    IntSet columns;
    for (int col = 0; col < length; col++) {
        columns.insert(col);
    }
    while (!columns.empty()) {
        SequenceSet occupied;
        int best_col = -1;
        int best_score = -1;
        BOOST_FOREACH (int c, columns) {
            if (!has_seq(c, occupied, bsrows)) {
                int score = column_size(c, bsrows);
                if (score > best_score) {
                    best_col = c;
                    best_score = score;
                }
            }
            add_seqs_to_set(occupied, c, bsrows);
        }
        ASSERT_NE(best_col, -1);
        ASSERT_NE(best_score, -1);
        append_col(new_bsrows, bsrows, best_col, size);
        columns.erase(best_col);
    }
    aln.swap(new_aln);
}

TreeNode* bsa_make_tree(const BSA& rows) {
    TreeNode* tree = new TreeNode;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, rows) {
        Sequence* seq = seq_and_row.first;
        const BSRow& row = seq_and_row.second;
        tree->add_child(new SequenceLeaf(seq, &row));
    }
    tree->upgma();
    return tree;
}

typedef std::map<std::string, Sequence*> Genome2Seq;

static TreeNode* bsa_convert_tree(const BSA& rows,
                                  const Genome2Seq& g2s,
                                  const TreeNode* tree) {
    const LeafNode* leaf = dynamic_cast<const LeafNode*>(tree);
    if (leaf) {
        Genome2Seq::const_iterator it_s = g2s.find(leaf->name());
        BOOST_ASSERT(it_s != g2s.end());
        Sequence* seq = it_s->second;
        ASSERT_TRUE(seq);
        BSA::const_iterator it_b = rows.find(seq);
        BOOST_ASSERT(it_b != rows.end());
        const BSRow& row = it_b->second;
        return new SequenceLeaf(seq, &row);
    } else {
        TreeNode* new_node = new TreeNode;
        BOOST_FOREACH (TreeNode* child, tree->children()) {
            new_node->add_child(bsa_convert_tree(rows, g2s, child));
        }
        return new_node;
    }
}

TreeNode* bsa_convert_tree(const BSA& rows, const TreeNode* tree) {
    Genome2Seq g2s;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, rows) {
        Sequence* seq = seq_and_row.first;
        if (g2s[seq->genome()]) {
            throw Exception("Two sequences from same genome: " +
                            seq->genome());
        }
        g2s[seq->genome()] = seq;
        if (g2s[seq->name()]) {
            throw Exception("Two sequences from same name: " +
                            seq->name());
        }
        g2s[seq->name()] = seq;
    }
    Leafs leafs;
    tree->all_leafs(leafs);
    if (g2s.size() != leafs.size() * 2) {
        throw Exception("Size of tree leafs != number of rows");
    }
    BOOST_FOREACH (LeafNode* leaf, leafs) {
        if (!g2s[leaf->name()]) {
            throw Exception("Leaf name " + leaf->name() +
                            " not found in rows");
        }
    }
    return bsa_convert_tree(rows, g2s, tree);
}

typedef std::vector<int> Starts;

static int sum_of_starts(const Starts& starts, int shift, int L) {
    int s = 0;
    BOOST_FOREACH (int start, starts) {
        start -= shift;
        if (start < 0) {
            start += L;
        } else if (start >= L) {
            start -= L;
        }
        s += start;
    }
    return s;
}

static void find_best_shift(BSA& bsa) {
    // build list of starts
    Starts starts;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, bsa) {
        const BSRow& bsrow = seq_and_row.second;
        for (int i = 0; i < bsrow.fragments.size(); i++) {
            Fragment* f = bsrow.fragments[i];
            if (f) {
                if (!f->prev() && !f->next()) {
                    // fragment is not connected
                    return;
                }
                int end_ori = (!f->prev()) ? -1 :
                              (!f->next()) ? 1 : 0;
                end_ori *= bsrow.ori;
                if (end_ori == -1) {
                    starts.push_back(i);
                }
            }
        }
    }
    // find best shift
    int L = bsa_length(bsa);
    int best_shift = 0;
    int best_sum = sum_of_starts(starts, 0, L);
    BOOST_FOREACH (int shift, starts) {
        int this_sum = sum_of_starts(starts, shift, L);
        if (this_sum < best_sum) {
            best_shift = shift;
            best_sum = this_sum;
        }
    }
    // apply best shift
    apply_shift(bsa, best_shift);
}

void bsa_orient(BSA& bsa) {
    int direct = 0, inverse = 0;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, bsa) {
        const BSRow& bsrow = seq_and_row.second;
        if (bsrow.ori == 1) {
            direct += 1;
        } else {
            inverse += 1;
        }
    }
    if (inverse > direct) {
        bsa_inverse(bsa);
    }
    if (bsa_is_circular(bsa)) {
        find_best_shift(bsa);
    }
}

void bsa_filter_exact_stem(BSA& bsa, int genomes) {
    BOOST_FOREACH (BSA::value_type& seq_and_row, bsa) {
        BSRow& row = seq_and_row.second;
        BOOST_FOREACH (Fragment*& fragment, row.fragments) {
            if (fragment) {
                Block* block = fragment->block();
                BOOST_ASSERT(block);
                if (!is_exact_stem(block, genomes)) {
                    fragment = 0;
                }
            }
        }
    }
}

void bsa_filter_long(BSA& bsa, int min_length) {
    BOOST_FOREACH (BSA::value_type& seq_and_row, bsa) {
        BSRow& row = seq_and_row.second;
        BOOST_FOREACH (Fragment*& fragment, row.fragments) {
            if (fragment) {
                Block* block = fragment->block();
                BOOST_ASSERT(block);
                if (block->alignment_length() < min_length) {
                    fragment = 0;
                }
            }
        }
    }
}

void bsa_print(std::ostream& out, const BSA& aln,
               const std::string& name,
               bool blocks, bool orientation) {
    typedef std::map<std::string, Sequence*> Name2Seq;
    Name2Seq n2s;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, aln) {
        Sequence* seq = seq_and_row.first;
        n2s[seq->name()] = seq;
    }
    // map is sorted
    BOOST_FOREACH (const Name2Seq::value_type& name_and_seq, n2s) {
        Sequence* seq = name_and_seq.second;
        out << name << '\t';
        const BSRow& row = aln.find(seq)->second;
        out << ((row.ori == 1) ? '+' : '-');
        out << seq->name();
        BOOST_FOREACH (Fragment* fragment, row.fragments) {
            out << '\t';
            if (fragment) {
                if (blocks) {
                    out << fragment->block()->name();
                } else {
                    out << fragment->id();
                }
                if (orientation) {
                    int ori = row.ori * fragment->ori();
                    out << ' ' << ((ori == 1) ? '>' : '<');
                }
            } else {
                out << '-';
            }
        }
        out << "\n";
    }
    out.flush();
}

void bsa_print_conservative(std::ostream& out, const BSA& aln,
                            const std::string& name) {
    if (aln.empty()) {
        return;
    }
    out << '#' << name;
    out << '\t' << "conservative";
    int length = bsa_length(aln);
    std::vector<BlockOri> conservative;
    const BSRow& first_row = aln.begin()->second;
    for (int col = 0; col < length; col++) {
        Fragment* fragment = first_row.fragments[col];
        Block* block = 0;
        int ori = 0;
        if (fragment) {
            block = fragment->block();
            BOOST_ASSERT(block);
            ori = fragment->ori() * first_row.ori;
        }
        conservative.push_back(BlockOri(block, ori));
    }
    BOOST_FOREACH (const BSA::value_type& seq_and_row, aln) {
        const BSRow& bsrow = seq_and_row.second;
        for (int col = 0; col < length; col++) {
            Fragment* fragment = bsrow.fragments[col];
            Block* block = 0;
            int ori = 0;
            if (fragment) {
                block = fragment->block();
                BOOST_ASSERT(block);
                ori = fragment->ori() * bsrow.ori;
            }
            if (conservative[col] != BlockOri(block, ori)) {
                conservative[col] = BlockOri();
            }
        }
    }
    for (int col = 0; col < length; col++) {
        out << '\t';
        Block* block = conservative[col].first;
        int ori = conservative[col].second;
        if (block) {
            out << block->name();
            out << ' ' << ((ori == 1) ? '>' : '<');
        } else {
            out << '-';
        }
    }
    out << "\n";
    out.flush();
}

const int BASE_INDEX = 2;

static bool match_parts(int shift, const Fragments& ff_orig,
                        const std::vector<std::string>& parts) {
    ASSERT_LT(shift, ff_orig.size());
    int orig_index = shift;
    for (int i = BASE_INDEX; i < parts.size(); i++) {
        const std::string& part = parts[i];
        if (part != "-") {
            if (orig_index >= ff_orig.size()) {
                orig_index -= ff_orig.size();
            }
            ASSERT_LT(orig_index, ff_orig.size());
            Fragment* f = ff_orig[orig_index];
            if (f->id() != part && f->block()->name() != part) {
                return false;
            }
            orig_index += 1;
        }
    }
    if (orig_index % ff_orig.size() != shift) {
        // not all fragments occured
        return false;
    }
    return true;
}

static void read_parts(int shift, const Fragments& ff_orig,
                       Fragments& ff_new,
                       const std::vector<std::string>& parts) {
    ASSERT_LT(shift, ff_orig.size());
    int orig_index = shift;
    for (int i = BASE_INDEX; i < parts.size(); i++) {
        const std::string& part = parts[i];
        if (part == "-") {
            ff_new.push_back(0);
        } else {
            if (orig_index >= ff_orig.size()) {
                orig_index -= ff_orig.size();
            }
            ASSERT_LT(orig_index, ff_orig.size());
            Fragment* f = ff_orig[orig_index];
            BOOST_ASSERT(f->id() == part ||
                         f->block()->name() == part);
            ff_new.push_back(f);
            orig_index += 1;
        }
    }
    ASSERT_EQ(orig_index % ff_orig.size(), shift);
    ASSERT_EQ(parts.size(), ff_new.size() + BASE_INDEX);
}

static void remove_orientation(Strings& parts) {
    BOOST_FOREACH (std::string& part, parts) {
        if (part.length() >= 3) {
            char last = part[part.length() - 1];
            char prelast = part[part.length() - 2];
            if (prelast == ' ' && (last == '>' || last == '<')) {
                part.resize(part.length() - 2);
            }
        }
    }
}

void bsa_input(BlockSet& bs, std::istream& in) {
    std::map<std::string, Sequence*> name2seq;
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        name2seq[seq->name()] = seq.get();
    }
    BSA rows;
    bsa_make_rows(rows, bs);
    for (std::string line; std::getline(in, line);) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        using namespace boost::algorithm;
        Strings parts;
        split(parts, line, is_any_of("\t"));
        remove_orientation(parts);
        if (parts.size() < 3) {
            continue;
        }
        const std::string& name = parts[0];
        BSA& bsa = bs.bsa(name);
        const std::string& ori_seq = parts[1];
        ASSERT_GTE(ori_seq.size(), 2);
        BOOST_ASSERT(ori_seq[0] == '+' || ori_seq[0] == '-');
        int ori = (ori_seq[0] == '+') ? (1) : (-1);
        std::string seq = ori_seq.substr(1);
        Sequence* s = name2seq[seq];
        ASSERT_TRUE(s);
        BOOST_ASSERT(rows.find(s) != rows.end());
        BSRow& bsrow_orig = rows[s];
        BSRow& bsrow_new = bsa[s];
        bsrow_new.ori = ori;
        Fragments& ff_orig = bsrow_orig.fragments;
        Fragments& ff_new = bsrow_new.fragments;
        ff_new.clear();
        if (ori == -1) {
            std::reverse(ff_orig.begin(), ff_orig.end());
        }
        if (s->circular()) {
            bool ok = false;
            for (int shift = 0; shift < ff_orig.size(); shift++) {
                if (match_parts(shift, ff_orig, parts)) {
                    read_parts(shift, ff_orig, ff_new, parts);
                    ok = true;
                    break;
                }
            }
            BOOST_ASSERT_MSG(ok, "bad match block set alignment");
        } else {
            BOOST_ASSERT(match_parts(0, ff_orig, parts));
            read_parts(0, ff_orig, ff_new, parts);
        }
    }
}

}

