/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <memory>
#include <algorithm>
#include <boost/cast.hpp>
#include <boost/foreach.hpp>

#include "MultipleAligner.hpp"
#include "ExpanderBase.hpp"
#include "PairAligner.hpp"
#include "tree.hpp"
#include "throw_assert.hpp"

namespace npge {

class SequenceNode : public TreeNode {
public:
    std::string consensus_;
    PairAlignment aln_;
    int level_;
    int index_in_seqs_;
};

static SequenceNode* align_two(PairAligner* pa,
                               SequenceNode* n1, SequenceNode* n2) {
    std::auto_ptr<SequenceNode> result((new SequenceNode));
    pa->set_first(n1->consensus_.c_str(), n1->consensus_.size());
    pa->set_second(n2->consensus_.c_str(), n2->consensus_.size());
    int first_last, second_last;
    std::string first_aln, second_aln;
    pa->align(first_last, second_last,
              &first_aln, &second_aln,
              &(result->aln_), '-');
    ASSERT_EQ(first_last, n1->consensus_.size() - 1);
    ASSERT_EQ(second_last, n2->consensus_.size() - 1);
    int length = result->aln_.size();
    ASSERT_GTE(length, n1->consensus_.size());
    ASSERT_GTE(length, n2->consensus_.size());
    ASSERT_EQ(first_aln.length(), length);
    ASSERT_EQ(second_aln.length(), length);
    if (n2->level_ > n1->level_) {
        first_aln.swap(second_aln);
    }
    for (int i = 0; i < length; i++) {
        if (first_aln[i] == '-') {
            first_aln[i] = second_aln[i];
            ASSERT_NE(first_aln[i], '-');
        }
    }
    result->consensus_.swap(first_aln);
    result->level_ = std::max(n1->level_, n2->level_) + 1;
    result->add_child(n1);
    result->add_child(n2);
    return result.release();
}

static TreeNode* build_trivial_tree(const Strings& seqs) {
    std::auto_ptr<TreeNode> tree((new TreeNode));
    for (int i = 0; i < seqs.size(); i++) {
        SequenceNode* leaf = new SequenceNode;
        leaf->consensus_ = seqs[i];
        leaf->level_ = 0;
        leaf->index_in_seqs_ = i;
        tree->add_child(leaf);
    }
    return tree.release();
}

static void build_tree(TreeNode* tree, PairAligner* pa) {
    typedef std::multimap<int, SequenceNode*> Level2Node;
    Level2Node l2n;
    BOOST_FOREACH (TreeNode* node, tree->children()) {
        SequenceNode* leaf;
        leaf = boost::polymorphic_downcast<SequenceNode*>(node);
        l2n.insert(std::make_pair(0, leaf));
    }
    while (l2n.size() >= 2) {
        // take two nodes with lowest level
        Level2Node::iterator begin = l2n.begin();
        SequenceNode* n1 = begin->second;
        l2n.erase(begin);
        tree->detach_child(n1);
        //
        begin = l2n.begin();
        SequenceNode* n2 = begin->second;
        l2n.erase(begin);
        tree->detach_child(n2);
        //
        SequenceNode* new_node = align_two(pa, n1, n2);
        l2n.insert(std::make_pair(new_node->level_, new_node));
        tree->add_child(new_node);
    }
}

static void change_seq(std::string& seq, const SequenceNode* sn) {
    TreeNode* parent = sn->parent();
    if (!parent) {
        return;
    }
    SequenceNode* sn_parent = dynamic_cast<SequenceNode*>(parent);
    if (!sn_parent) {
        return;
    }
    ASSERT_EQ(sn_parent->children().size(), 2);
    bool first = (sn_parent->children().front() == sn);
    const PairAlignment& aln = sn_parent->aln_;
    ASSERT_GTE(aln.size(), seq.size());
    // copy parts of sequence to higher positions according to aln
    int seq_pos = seq.size() - 1;
    seq.resize(aln.size());
    for (int aln_pos = aln.size() - 1; aln_pos >= 0; aln_pos--) {
        if (aln_pos == seq_pos) {
            // optimization
            seq_pos = -1;
            break;
        }
        char c = '-';
        const AlignmentPair& pair = aln[aln_pos];
        int pos = first ? pair.first : pair.second;
        if (pos != -1) {
            ASSERT_LTE(pos, seq_pos);
            ASSERT_GTE(pos, 0);
            c = seq[seq_pos];
            seq_pos -= 1;
        }
        seq[aln_pos] = c;
    }
    ASSERT_EQ(seq_pos, -1);
    ASSERT_EQ(seq.size(), aln.size());
    change_seq(seq, sn_parent);
}

static void change_seqs(Strings& seqs, const TreeNode* tree) {
    ASSERT_EQ(tree->children().size(), 1);
    TreeNode* top = tree->children().front();
    SequenceNode* top_sn = dynamic_cast<SequenceNode*>(top);
    ASSERT_TRUE(top_sn);
    int length = top_sn->consensus_.length();
    Nodes leafs;
    tree->all_end_nodes(leafs);
    BOOST_FOREACH (TreeNode* leaf, leafs) {
        SequenceNode* sn;
        sn = boost::polymorphic_downcast<SequenceNode*>(leaf);
        ASSERT_LT(sn->index_in_seqs_, seqs.size());
        std::string& seq = seqs[sn->index_in_seqs_];
        seq.reserve(length);
        change_seq(seq, sn);
    }
}

void MultipleAligner::multiple_aligner(Strings& seqs,
                                       PairAligner* pa) {
    pa->set_max_errors(-1);
    pa->set_no_tail(false);
    std::auto_ptr<TreeNode> tree((build_trivial_tree(seqs)));
    build_tree(tree.get(), pa);
    change_seqs(seqs, tree.get());
}

MultipleAligner::MultipleAligner() {
    add_expander_options(this);
}

void MultipleAligner::align_seqs_impl(Strings& seqs) const {
    PairAligner pa;
    apply_pair_aligner_options(&pa, this);
    multiple_aligner(seqs, &pa);
}

const char* MultipleAligner::name_impl() const {
    return "Internal aligner";
}

}

