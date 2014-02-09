/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <vector>
#include <algorithm>
#include <boost/cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>

#include "ConsensusTree.hpp"
#include "BlocksJobs.hpp"
#include "PrintTree.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "block_stat.hpp"
#include "Exception.hpp"
#include "tree.hpp"

namespace bloomrepeats {

class BranchData : public ThreadData {
public:
    BranchTable table;
};

class BranchGenerator : public BlocksJobs {
public:
    mutable BranchTable table;

    BranchGenerator() {
        print_tree_  = new PrintTree;
        print_tree_->set_parent(this);
        add_opt("min-noident",
                "Minimal number of nonidentical positions in block", 5);
        add_opt("log",
                "block weight is log(number of nonidentical positions) "
                "(otherwse linear)", true);
    }

    bool initialize_work() const {
        table.clear();
    }

    ThreadData* before_thread_impl() const {
        return new BranchData;
    }

    struct GenomeNameCompare {
        bool operator()(const LeafNode* l1, const LeafNode* l2) const {
            const FragmentLeaf* fl1, *fl2;
            fl1 = boost::polymorphic_downcast<const FragmentLeaf*>(l1);
            fl2 = boost::polymorphic_downcast<const FragmentLeaf*>(l2);
            const Fragment* f1 = fl1->fragment();
            const Fragment* f2 = fl2->fragment();
            return f1->seq()->genome() < f2->seq()->genome();
        }
    };

    bool process_block_impl(Block* block, ThreadData* data) const {
        BranchData* d = boost::polymorphic_downcast<BranchData*>(data);
        AlignmentStat stat;
        make_stat(stat, block);
        double block_weight = stat.noident_nogap() + stat.noident_gap();
        if (block_weight < opt_value("min-noident").as<int>()) {
            return false;
        }
        if (opt_value("log").as<bool>()) {
            double block_weight = log(block_weight);
        }
        boost::scoped_ptr<TreeNode> tree(print_tree_->make_tree(block));
        Leafs leafs;
        tree->all_leafs(leafs);
        std::sort(leafs.begin(), leafs.end(), GenomeNameCompare());
        tree->branch_table(d->table, leafs, block_weight);
    }

    bool after_thread_impl(ThreadData* data) const {
        BranchData* d = boost::polymorphic_downcast<BranchData*>(data);
        BOOST_FOREACH (const BranchTable::value_type& branch_length, d->table) {
            table[branch_length.first] += branch_length.second;
        }
    }

private:
    PrintTree* print_tree_;
};

ConsensusTree::ConsensusTree() {
    branch_generator_ = new BranchGenerator;
    branch_generator_->set_parent(this);
}

bool ConsensusTree::run_impl() const {
    branch_generator_->apply(block_set());
    typedef std::pair<double, std::string> Weight_Branch;
    typedef std::vector<Weight_Branch> BranchVector;
    BranchVector branch_vector;
    BOOST_FOREACH (const BranchTable::value_type& branch_length,
            branch_generator_->table) {
        branch_vector.push_back(Weight_Branch(branch_length.second,
                    branch_length.first));
    }
    std::sort(branch_vector.rbegin(), branch_vector.rend()); // reverse
    BranchVector compatible_branches;
    BOOST_FOREACH (const Weight_Branch& branch, branch_vector) {
        bool compatible = true;
        BOOST_FOREACH (const Weight_Branch& branch0, compatible_branches) {
            if (!TreeNode::branches_compatible(branch.second, branch0.second)) {
                compatible = false;
                break;
            }
        }
        if (compatible) {
            compatible_branches.push_back(branch);
            std::cerr << branch.second << "\n";
            // TODO
        }
    }
}

}

