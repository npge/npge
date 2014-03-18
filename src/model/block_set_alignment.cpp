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

#include "block_set_alignment.hpp"
#include "GeneralAligner.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "tree.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

typedef std::pair<Block*, int> BlockOri;
typedef std::set<BlockOri> BlockOriSet;

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

void bsa_make_rows(BSA& rows, const BlockSet& bs) {
    if (rows.empty()) {
        BOOST_FOREACH (const SequencePtr& seq, bs.seqs()) {
            rows[seq.get()] = BSRow();
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
                    return -1 - log(block->alignment_length());
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
    ga.set_local(true);
    int first_size = bsa_length(first);
    int second_size = bsa_length(second);
    ga.set_gap_penalty(5);
    ga.set_gap_range(std::max(first_size, second_size));
    // ^^ FIXME GeneralAligner full matrix
    ga.set_contents(BSContents(first, second));
    int first_last, second_last;
    ga.align(first_last, second_last);
    score = ga.opt_score();
    PairAlignment alignment;
    first_last = first_size - 1;
    second_last = second_size - 1;
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

void bsa_make_aln(BSA& aln, const BSAs& parts) {
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
            bsa_align(both_direct, score_direct, aln, second);
        }
        {
            BSA second = parts[i];
            bsa_inverse(second);
            bsa_align(both_inverse, score_inverse, aln, second);
        }
        bool use_direct = (score_direct < score_inverse);
        BSA& both = use_direct ? both_direct : both_inverse;
        aln.swap(both);
    }
}

void bsa_make_aln(BSA& aln, const BSA& rows) {
    BSAs parts;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, rows) {
        Sequence* seq = seq_and_row.first;
        const BSRow& row = seq_and_row.second;
        parts.push_back(BSA());
        parts.back()[seq] = row;
    }
    bsa_make_aln(aln, parts);
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

static void bsa_make_aln_by_tree(BSA& aln, const TreeNode* tree) {
    const SequenceLeaf* seq_leaf;
    seq_leaf = dynamic_cast<const SequenceLeaf*>(tree);
    if (seq_leaf) {
        aln[seq_leaf->seq()] = *seq_leaf->bsrow();
    } else {
        BSAs parts;
        BOOST_FOREACH (TreeNode* child, tree->children()) {
            parts.push_back(BSA());
            bsa_make_aln_by_tree(parts.back(), child);
        }
        bsa_make_aln(aln, parts);
    }
}

void bsa_make_aln_by_tree(BSA& aln, const BSA& rows,
                          const TreeNode* tree0) {
    boost::scoped_ptr<TreeNode> tree((bsa_convert_tree(rows, tree0)));
    bsa_make_aln_by_tree(aln, tree.get());
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

static int count_block_ori(const BSRows& bsrows, int col,
                           Block* block, int ori) {
    int result = 0;
    BOOST_FOREACH (const BSRow* bsrow, bsrows) {
        Fragment* f = bsrow->fragments[col];
        if (f && f->block() == block && f->ori() * bsrow->ori == ori) {
            result += 1;
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
    int score = count_block_ori(bsrows, col, block, ori) - 1;
    // -1 because f itself was counted
    BOOST_ASSERT(score >= 0);
    int best_score = score;
    int best_col = col;
    for (int i = col - 1; i >= 0; i--) {
        if (bsrow.fragments[i]) {
            break;
        }
        int ascore = count_block_ori(bsrows, i, block, ori);
        if (ascore > best_score) {
            best_col = i;
            best_score = ascore;
        }
    }
    for (int i = col + 1; i < bsrow.fragments.size(); i++) {
        if (bsrow.fragments[i]) {
            break;
        }
        int ascore = count_block_ori(bsrows, i, block, ori);
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

TreeNode* bsa_make_tree(const BSA& rows) {
    TreeNode* tree = new TreeNode;
    BOOST_FOREACH (const BSA::value_type& seq_and_row, rows) {
        Sequence* seq = seq_and_row.first;
        const BSRow& row = seq_and_row.second;
        tree->add_child(new SequenceLeaf(seq, &row));
    }
    tree->neighbor_joining();
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
        BOOST_ASSERT(seq);
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

void bsa_print(std::ostream& out, const BSA& aln,
               const std::string& name,
               bool blocks) {
    BOOST_FOREACH (const BSA::value_type& seq_and_row, aln) {
        out << name << '\t';
        Sequence* seq = seq_and_row.first;
        const BSRow& row = seq_and_row.second;
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
            } else {
                out << '-';
            }
        }
        out << "\n";
    }
    out.flush();
}

void bsa_input(BlockSet& bs, std::istream& in) {
    BSA rows;
    bsa_make_rows(rows, bs);
    for (std::string line; std::getline(in, line);) {
        using namespace boost::algorithm;
        std::vector<std::string> parts;
        split(parts, line, is_any_of("\t"));
        if (parts.size() < 3) {
            continue;
        }
        const std::string& name = parts[0];
        BSA& bsa = bs.bsa(name);
        const std::string& ori_seq = parts[1];
        BOOST_ASSERT(ori_seq.size() >= 2);
        BOOST_ASSERT(ori_seq[0] == '+' || ori_seq[0] == '-');
        int ori = (ori_seq[0] == '+') ? (1) : (-1);
        std::string seq = ori_seq.substr(1);
        SequencePtr s = bs.seq_from_name(seq);
        BOOST_ASSERT(s);
        BOOST_ASSERT(rows.find(s.get()) != rows.end());
        BSRow& bsrow_orig = rows[s.get()];
        BSRow& bsrow_new = bsa[s.get()];
        bsrow_new.ori = ori;
        Fragments& ff_orig = bsrow_orig.fragments;
        Fragments& ff_new = bsrow_new.fragments;
        ff_new.clear();
        if (ori == -1) {
            std::reverse(ff_orig.begin(), ff_orig.end());
        }
        int orig_index = 0;
        for (int i = 2; i < parts.size(); i++) {
            const std::string& part = parts[i];
            if (part == "-") {
                ff_new.push_back(0);
            } else {
                BOOST_ASSERT(orig_index < ff_orig.size());
                Fragment* f = ff_orig[orig_index];
                BOOST_ASSERT(f->id() == part ||
                             f->block()->name() == part);
                ff_new.push_back(f);
                orig_index += 1;
            }
        }
        BOOST_ASSERT(orig_index == ff_orig.size());
        BOOST_ASSERT(parts.size() == ff_new.size() + 2);
    }
}

}

