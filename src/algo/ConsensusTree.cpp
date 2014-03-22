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
#include <boost/algorithm/string/join.hpp>

#include "ConsensusTree.hpp"
#include "Union.hpp"
#include "Stem.hpp"
#include "BlocksJobs.hpp"
#include "PrintTree.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "block_stat.hpp"
#include "Exception.hpp"
#include "tree.hpp"
#include "boost-assert.hpp"

namespace bloomrepeats {

typedef std::map<std::string, double> LeafLength;
typedef std::vector<Block*> Blocks;
typedef std::map<std::string, Blocks> BranchBlocks;

class BranchData : public ThreadData {
public:
    BranchTable table;
    BranchBlocks branch_blocks;
    LeafLength leaf_length;
};

class BranchGenerator : public BlocksJobs {
public:
    mutable BranchTable table;
    mutable BranchBlocks branch_blocks;
    mutable LeafLength leaf_length;

    BranchGenerator() {
        print_tree_  = new PrintTree;
        print_tree_->set_parent(this);
        add_opt("min-noident",
                "Minimal number of nonidentical positions in block", 5);
        add_opt("log",
                "block weight is log(number of nonidentical positions) "
                "(otherwse linear)", true);
    }

    bool initialize_work_impl() const {
        table.clear();
        branch_blocks.clear();
        leaf_length.clear();
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

    static void add_table(BranchTable& dst, const BranchTable& src) {
        BOOST_FOREACH (const BranchTable::value_type& branch_length,
                      src) {
            dst[branch_length.first] += branch_length.second;
        }
    }

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
        BranchTable t;
        tree->branch_table(t, leafs, block_weight);
        add_table(d->table, t);
        BOOST_FOREACH (const BranchTable::value_type& branch_length,
                      t) {
            d->branch_blocks[branch_length.first].push_back(block);
        }
        BOOST_FOREACH (LeafNode* leaf, leafs) {
            FragmentLeaf* fl;
            fl = boost::polymorphic_downcast<FragmentLeaf*>(leaf);
            std::string genome = fl->fragment()->seq()->genome();
            d->leaf_length[genome] += leaf->length() * block_weight;
        }
    }

    bool after_thread_impl(ThreadData* data) const {
        BranchData* d = boost::polymorphic_downcast<BranchData*>(data);
        add_table(table, d->table);
        BOOST_FOREACH (const LeafLength::value_type& ll,
                      d->leaf_length) {
            leaf_length[ll.first] += ll.second;
        }
        BOOST_FOREACH (const BranchBlocks::value_type& bb,
                      d->branch_blocks) {
            const std::string& branch_str = bb.first;
            const Blocks& blocks = bb.second;
            Blocks& dst_blocks = branch_blocks[branch_str];
            BOOST_FOREACH (Block* block, blocks) {
                dst_blocks.push_back(block);
            }
        }
    }

private:
    PrintTree* print_tree_;
};

ConsensusTree::ConsensusTree():
    file_writer_(this, "out-consensus-tree", "Output file with statistics") {
    branch_generator_ = new BranchGenerator;
    branch_generator_->set_parent(this);
}

static std::vector<std::string> genomes_list(BlockSetPtr bs) {
    if (bs->empty()) {
        return std::vector<std::string>();
    }
    std::set<std::string> genomes;
    BOOST_FOREACH (Fragment* f, *bs->front()) {
        genomes.insert(f->seq()->genome());
    }
    std::vector<std::string> genomes_v(genomes.begin(),
                                       genomes.end());
    std::sort(genomes_v.begin(), genomes_v.end()); // useless
    return genomes_v;
}

class GenomeLeaf : public LeafNode {
public:
    GenomeLeaf(const std::string& g):
        genome_(g)
    { }

    double distance_to_impl(const LeafNode* leaf) const {
        // useless
        return -1.0;
    }

    std::string name_impl() const {
        return genome_;
    }

    TreeNode* clone_impl() const {
        return new GenomeLeaf(genome_);
    }

private:
    std::string genome_;
};

typedef std::map<std::string, GenomeLeaf*> Genome2Leaf;

static int branch_size(const std::string& branch) {
    int c0 = 0, c1 = 0;
    BOOST_FOREACH (char c, branch) {
        if (c == '0') {
            c0 += 1;
        } else {
            c1 += 1;
        }
    }
    return std::min(c0, c1);
}

typedef std::pair<double, std::string> Weight_Branch;

struct BranchCompare {
    bool operator()(const Weight_Branch& b1,
                    const Weight_Branch& b2) {
        return branch_size(b1.second) < branch_size(b2.second);
    }
};

static TreeNode* ancestor(TreeNode* node, TreeNode* tree) {
    while (node->parent() != tree) {
        BOOST_ASSERT(node->parent());
        node = node->parent();
    }
    return node;
}

bool ConsensusTree::run_impl() const {
    // convert to pipe
    std::ostream& out = file_writer_.output();
    Union copy(block_set());
    copy.run();
    copy.block_set()->add_sequences(block_set()->seqs());
    Stem stem;
    stem.set_opt_value("exact", true);
    stem.apply(copy.block_set());
    std::vector<std::string> genomes_v = genomes_list(copy.block_set());
    branch_generator_->apply(copy.block_set());
    typedef std::vector<Weight_Branch> BranchVector;
    BranchVector branch_vector;
    BOOST_FOREACH (const BranchTable::value_type& branch_length,
                  branch_generator_->table) {
        branch_vector.push_back(Weight_Branch(branch_length.second,
                                              branch_length.first));
    }
    std::sort(branch_vector.rbegin(), branch_vector.rend()); // reverse
    TreeNode cons_tree;
    Leafs cons_leafs;
    LeafLength& leaf_length = branch_generator_->leaf_length;
    Genome2Leaf g2f;
    BOOST_FOREACH (std::string genome, genomes_v) {
        GenomeLeaf* leaf = new GenomeLeaf(genome);
        leaf->set_length(leaf_length[genome]);
        cons_tree.add_child(leaf);
        g2f[genome] = leaf;
        cons_leafs.push_back(leaf);
    }
    BranchVector compatible_branches;
    BOOST_FOREACH (const Weight_Branch& branch, branch_vector) {
        bool compatible = true;
        BOOST_FOREACH (const Weight_Branch& branch0, compatible_branches) {
            if (!TreeNode::branches_compatible(branch.second, branch0.second)) {
                compatible = false;
                break;
            }
        }
        Blocks& blocks = branch_generator_->branch_blocks[branch.second];
        std::vector<std::string> block_names;
        BOOST_FOREACH (Block* block, blocks) {
            block_names.push_back(block->name());
        }
        using namespace boost::algorithm;
        std::string blocks_str = join(block_names, ",");
        if (compatible) {
            compatible_branches.push_back(branch);
            out
                    << TreeNode::branch_as_sets(cons_leafs, branch.second)
                    << " weight=" << branch.first << "\n";
        } else {
            out << "Incompatible branch: "
                << TreeNode::branch_as_sets(cons_leafs, branch.second)
                << " weight=" << branch.first << "\n";
        }
        out << "blocks: " << blocks_str << "\n";
    }
    std::sort(compatible_branches.begin(), compatible_branches.end(),
              BranchCompare());
    BOOST_FOREACH (const Weight_Branch& branch, compatible_branches) {
        std::set<TreeNode*> nodes0, nodes1;
        double length = branch.first;
        const std::string& branch_str = branch.second;
        BOOST_ASSERT(branch_str.size() == genomes_v.size());
        for (int i = 0; i < branch_str.size(); i++) {
            char c = branch_str[i];
            std::set<TreeNode*>& nodes = (c == '0') ? nodes0 : nodes1;
            nodes.insert(ancestor(g2f[genomes_v[i]], &cons_tree));
        }
        bool n0 = (nodes0.size() < nodes1.size());
        std::set<TreeNode*>& nodes = n0 ? nodes0 : nodes1;
        TreeNode* branch_node = new TreeNode;
        cons_tree.add_child(branch_node);
        branch_node->set_length(length);
        BOOST_FOREACH (TreeNode* node, nodes) {
            branch_node->add_child(node);
        }
    }
    out << cons_tree.newick() << "\n";
}

}

