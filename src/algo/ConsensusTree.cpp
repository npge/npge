/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <vector>
#include <algorithm>
#include <boost/bind.hpp>
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
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

typedef std::map<std::string, double> LeafLength;
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
                "Minimum number of nonidentical positions in block", 5);
        add_opt("log",
                "block weight is log(number of nonidentical positions) "
                "(otherwse linear)", true);
    }

    void initialize_work_impl() const {
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

    void process_block_impl(Block* block, ThreadData* data) const {
        BranchData* d = boost::polymorphic_downcast<BranchData*>(data);
        AlignmentStat stat;
        make_stat(stat, block);
        double block_weight = stat.noident_nogap() + stat.noident_gap();
        if (block_weight < opt_value("min-noident").as<int>()) {
            return;
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

    void after_thread_impl(ThreadData* data) const {
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

static bool check_bootstrap_values(std::string& message,
                                   Processor* p) {
    std::string v = p->opt_value("bootstrap-values").as<std::string>();
    if (v != "blocks" && v != "length" &&
            v != "diagnostic-positions") {
        message = "bad bootstrap-values";
        return false;
    }
    return true;
}

static bool check_bootstrap_print(std::string& message,
                                  Processor* p) {
    std::string v = p->opt_value("bootstrap-print").as<std::string>();
    if (v != "no" && v != "in-braces" && v != "before-length") {
        message = "bad bootstrap-print";
        return false;
    }
    return true;
}

ConsensusTree::ConsensusTree():
    branch_(this, "out-branch",
            "Output file with branches"),
    tre_(this, "out-consensus-tree",
         "Output file with consensus tree") {
    branch_generator_ = new BranchGenerator;
    branch_generator_->set_parent(this);
    add_opt("bootstrap-values", "What to use as bootstrap values "
            "('blocks', 'length', 'diagnostic-positions')",
            std::string("diagnostic-positions"));
    add_opt("bootstrap-diagnostic-stem", "If only stem blocks are used "
            "to calculate number of diagnostic positions", true);
    add_opt("bootstrap-diagnostic-min-block", "Minimum block used "
            "to calculate number of diagnostic positions", 4);
    add_opt_check(boost::bind(check_bootstrap_values, _1, this));
    add_opt("bootstrap-percent", "If bootstrap is percentage", false);
    add_opt("bootstrap-print", "How to print bootstrap values "
            "('no', 'in-braces', 'before-length')",
            std::string("in-braces"));
    add_opt_check(boost::bind(check_bootstrap_print, _1, this));
    add_opt("tree-pseudo-leafs", "Convert each leaf to short branch "
            "(workaround to make viewer programs show bootstrap "
            "values of leafs)", false);
    add_opt("print-branches", "If branches are logged", true);
    declare_bs("target", "Target blockset");
}

static Strings genomes_list(BlockSetPtr bs) {
    if (bs->empty()) {
        return Strings();
    }
    std::set<std::string> genomes;
    BOOST_FOREACH (Fragment* f, *bs->front()) {
        genomes.insert(f->seq()->genome());
    }
    Strings genomes_v(genomes.begin(), genomes.end());
    std::sort(genomes_v.begin(), genomes_v.end()); // useless
    return genomes_v;
}

class GenomeLeaf : public LeafNode {
public:
    GenomeLeaf(const std::string& g):
        genome_(g) {
    }

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
        ASSERT_TRUE(node->parent());
        node = node->parent();
    }
    return node;
}

static TreeNode::ShowBootstrap parse_bp(const std::string& bp) {
    if (bp == "no") {
        return TreeNode::NO_BOOTSTRAP;
    } else if (bp == "in-braces") {
        return TreeNode::BOOTSTRAP_IN_BRACES;
    } else if (bp == "before-length") {
        return TreeNode::BOOTSTRAP_BEFORE_LENGTH;
    } else {
        throw Exception("Bad bootstrap-print");
    }
}

enum BootstrapValues {
    BLOCKS,
    LENGTH,
    DIAGNOSTIC_POSITIONS
};

static BootstrapValues parse_bv(const std::string& bv) {
    if (bv == "blocks") {
        return BLOCKS;
    } else if (bv == "length") {
        return LENGTH;
    } else if (bv == "diagnostic-positions") {
        return DIAGNOSTIC_POSITIONS;
    } else {
        throw Exception("Bad bootstrap-values");
    }
}

class BootstrapDiagnosticPositions : public BlocksJobs {
public:
    typedef std::set<std::string> StringSet;
    typedef std::pair<TreeNode*, StringSet> Clade;
    typedef std::vector<Clade> Clades;
    typedef std::map<TreeNode*, double> Bootstrap;

    void set_tree(TreeNode* cons_tree) {
        cons_tree_ = cons_tree;
    }

    void set_min_block_size(int min_block_size) {
        min_block_size_ = min_block_size;
    }

private:
    mutable Nodes nodes_;
    mutable Clades clades_;
    mutable Bootstrap bootstrap_;

    TreeNode* cons_tree_;
    int min_block_size_;

    struct BPSData : public ThreadData {
        Bootstrap bootstrap_;
    };

protected:
    void initialize_work_impl() const {
        cons_tree_->all_descendants(nodes_);
        BOOST_FOREACH (TreeNode* node, nodes_) {
            Leafs leafs;
            node->all_leafs_and_this(leafs);
            StringSet genomes;
            BOOST_FOREACH (LeafNode* leaf, leafs) {
                genomes.insert(leaf->name());
            }
            Clade clade((node), genomes);
            clades_.push_back(clade);
        }
    }

    ThreadData* before_thread_impl() const {
        return new BPSData;
    }

    void process_block_impl(Block* block, ThreadData* data) const {
        if (block->size() < min_block_size_) {
            return;
        }
        int length = block->alignment_length();
        BPSData* d = boost::polymorphic_downcast<BPSData*>(data);
        Bootstrap& b = d->bootstrap_;
        BOOST_FOREACH (const Clade& clade, clades_) {
            TreeNode* node = clade.first;
            const StringSet& genomes = clade.second;
            Fragments ff, other;
            BOOST_FOREACH (Fragment* f, *block) {
                ASSERT_TRUE(f->seq());
                std::string genome = f->seq()->genome();
                if (genomes.find(genome) != genomes.end()) {
                    ff.push_back(f);
                } else {
                    other.push_back(f);
                }
            }
            if (!ff.empty() && !other.empty()) {
                int diagnostic_cols = 0;
                for (int col = 0; col < length; col++) {
                    if (is_diagnostic(col, ff, other)) {
                        diagnostic_cols += 1;
                    }
                }
                b[node] += diagnostic_cols;
            }
        }
    }

    void after_thread_impl(ThreadData* data) const {
        BPSData* d = boost::polymorphic_downcast<BPSData*>(data);
        Bootstrap& b = d->bootstrap_;
        BOOST_FOREACH (TreeNode* node, nodes_) {
            bootstrap_[node] += b[node];
        }
    }

    void finish_work_impl() const {
        BOOST_FOREACH (TreeNode* node, nodes_) {
            node->set_bootstrap(bootstrap_[node]);
        }
    }
};

void ConsensusTree::run_impl() const {
    // TODO convert to pipe
    std::string bv = opt_value("bootstrap-values").as<std::string>();
    bool print_branches = opt_value("print-branches").as<bool>();
    BootstrapValues bsv = parse_bv(bv);
    std::ostream& branch_out = branch_.output();
    std::ostream& out = tre_.output();
    Union copy(block_set());
    copy.run();
    copy.block_set()->add_sequences(block_set()->seqs());
    Stem stem;
    stem.set_opt_value("exact", true);
    stem.apply(copy.block_set());
    Strings genomes_v = genomes_list(copy.block_set());
    branch_generator_->apply(copy.block_set());
    typedef std::vector<Weight_Branch> BranchVector;
    BranchVector branch_vector;
    BOOST_FOREACH (const BranchTable::value_type& branch_length,
                  branch_generator_->table) {
        branch_vector.push_back(Weight_Branch(branch_length.second,
                                              branch_length.first));
    }
    std::sort(branch_vector.rbegin(), branch_vector.rend()); // reverse
    boost::shared_ptr<TreeNode> cons_tree((new TreeNode));
    Leafs cons_leafs;
    LeafLength& leaf_length = branch_generator_->leaf_length;
    BranchBlocks& branch_blocks = branch_generator_->branch_blocks;
    Genome2Leaf g2f;
    BOOST_FOREACH (std::string genome, genomes_v) {
        GenomeLeaf* leaf = new GenomeLeaf(genome);
        leaf->set_length(leaf_length[genome]);
        cons_tree->add_child(leaf);
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
        Blocks& blocks = branch_blocks[branch.second];
        if (compatible) {
            compatible_branches.push_back(branch);
            if (print_branches) {
                branch_out
                        << TreeNode::branch_as_sets(cons_leafs, branch.second)
                        << " weight=" << branch.first << "\n";
            }
        } else if (print_branches) {
            branch_out
                    << "Incompatible branch: "
                    << TreeNode::branch_as_sets(cons_leafs, branch.second)
                    << " weight=" << branch.first << "\n";
        }
        if (print_branches) {
            Strings block_names;
            BOOST_FOREACH (Block* block, blocks) {
                block_names.push_back(block->name());
            }
            using namespace boost::algorithm;
            std::string blocks_str = join(block_names, ",");
            branch_out
                    << "blocks (" << blocks.size()
                    << "): " << blocks_str << "\n";
        }
    }
    std::sort(compatible_branches.begin(), compatible_branches.end(),
              BranchCompare());
    BOOST_FOREACH (const Weight_Branch& branch, compatible_branches) {
        std::set<TreeNode*> nodes0, nodes1;
        double length = branch.first;
        const std::string& branch_str = branch.second;
        ASSERT_EQ(branch_str.size(), genomes_v.size());
        for (int i = 0; i < branch_str.size(); i++) {
            char c = branch_str[i];
            std::set<TreeNode*>& nodes = (c == '0') ? nodes0 : nodes1;
            nodes.insert(ancestor(g2f[genomes_v[i]], cons_tree.get()));
        }
        bool n0 = (nodes0.size() < nodes1.size());
        std::set<TreeNode*>& nodes = n0 ? nodes0 : nodes1;
        TreeNode* branch_node = new TreeNode;
        cons_tree->add_child(branch_node);
        branch_node->set_length(length);
        if (bsv == BLOCKS) {
            const Blocks& blocks = branch_blocks[branch.second];
            branch_node->set_bootstrap(blocks.size());
        } else if (bsv == LENGTH) {
            const Blocks& blocks = branch_blocks[branch.second];
            int sum = 0;
            BOOST_FOREACH (Block* block, blocks) {
                sum += block->alignment_length();
            }
            branch_node->set_bootstrap(sum);
        }
        BOOST_FOREACH (TreeNode* node, nodes) {
            branch_node->add_child(node);
        }
    }
    if (bsv == DIAGNOSTIC_POSITIONS) {
        bool stem_only = opt_value("bootstrap-diagnostic-stem").as<bool>();
        int min_block = opt_value("bootstrap-diagnostic-min-block").as<int>();
        BlockSetPtr diagnostic_bs = stem_only ? copy.block_set() : block_set();
        BootstrapDiagnosticPositions bdp;
        bdp.set_block_set(diagnostic_bs);
        bdp.set_tree(cons_tree.get());
        bdp.set_min_block_size(min_block);
        bdp.set_workers(workers());
        bdp.run();
    }
    bool bootstrap_percent = opt_value("bootstrap-percent").as<bool>();
    if (bootstrap_percent) {
        double sum = 0;
        BOOST_FOREACH (TreeNode* node, cons_tree->children()) {
            sum += node->bootstrap();
        }
        if (sum > 0.0001) {
            double factor = 100.0 / sum;
            Nodes all_nodes;
            cons_tree->all_nodes(all_nodes);
            BOOST_FOREACH (TreeNode* node, all_nodes) {
                node->set_bootstrap(node->bootstrap() * factor);
            }
        }
    }
    bool pseudo_leafs = opt_value("tree-pseudo-leafs").as<bool>();
    if (pseudo_leafs) {
        TreeNode* clone = cons_tree->clone_with_pseudo_leafs();
        cons_tree.reset(clone);
    }
    bool lengthes = true;
    std::string bp = opt_value("bootstrap-print").as<std::string>();
    TreeNode::ShowBootstrap sbs = parse_bp(bp);
    out << cons_tree->newick(lengthes, sbs) << "\n";
}

const char* ConsensusTree::name_impl() const {
    return "Print consensus tree";
}

}

