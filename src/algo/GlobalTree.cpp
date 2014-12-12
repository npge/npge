/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <boost/bind.hpp>

#include "GlobalTree.hpp"
#include "ConsensusTree.hpp"
#include "FragmentDistance.hpp"
#include "RemoveNonStem.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "block_stat.hpp"
#include "tree.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"
#include "global.hpp"

namespace npge {

GlobalTree::GlobalTree():
    tre_(this, "out-global-tree",
         "Output file with global tree") {
    add_opt("bootstrap-print", "How to print bootstrap values "
            "('no', 'in-braces', 'before-length')",
            std::string("before-length"));
    add_opt_check(boost::bind(check_bootstrap_print, _1, this));
    add_opt("tree-pseudo-leafs", "Convert each leaf to short branch "
            "(workaround to make viewer programs show bootstrap "
            "values of leafs)", false);
    declare_bs("target", "Target blockset");
}

typedef std::map<std::string, double> Genome2Double;
typedef std::map<std::string, Genome2Double> Dist;

class GenomeLeaf : public LeafNode {
public:
    GenomeLeaf(const std::string& g, const Dist* dist):
        genome_(g), dist_(dist) {
    }

    double distance_to_impl(const LeafNode* l0) const {
        const GenomeLeaf* leaf = D_CAST<const GenomeLeaf*>(l0);
        Dist& dist = *(const_cast<Dist*>(dist_));
        return dist[genome_][leaf->genome_];
    }

    std::string name_impl() const {
        return genome_;
    }

    TreeNode* clone_impl() const {
        return new GenomeLeaf(genome_, dist_);
    }

private:
    std::string genome_;
    const Dist* dist_;
};

void add_dist(Dist& dist, FragmentDistance* d, Block* block) {
    Fragments ff(block->begin(), block->end());
    for (int i = 0; i < ff.size(); i++) {
        Fragment* f1 = ff[i];
        std::string genome1 = f1->seq()->genome();
        for (int j = 0; j < i; j++) {
            Fragment* f2 = ff[j];
            std::string genome2 = f2->seq()->genome();
            FragmentDistance::Distance di;
            di = d->fragment_distance(f1, f2);
            int mutations = di.penalty;
            dist[genome1][genome2] += mutations;
            dist[genome2][genome1] += mutations;
        }
    }
}

void GlobalTree::run_impl() const {
    std::ostream& out = tre_.output();
    BlockSetPtr copy = block_set()->clone();
    RemoveNonStem stem;
    stem.set_opt_value("exact", true);
    stem.apply(copy);
    Dist dist;
    FragmentDistance d;
    BOOST_FOREACH (Block* block, *copy) {
        add_dist(dist, &d, block);
    }
    boost::shared_ptr<TreeNode> tree(new TreeNode);
    Strings genomes_v = genomes_list(copy);
    BOOST_FOREACH (std::string genome, genomes_v) {
        GenomeLeaf* leaf = new GenomeLeaf(genome, &dist);
        tree->add_child(leaf);
    }
    tree->neighbor_joining();
    //
    add_diagnostic(tree.get(), copy,
                   genomes_v.size(), workers());
    //
    if (opt_value("tree-pseudo-leafs").as<bool>()) {
        TreeNode* clone = tree->clone_with_pseudo_leafs();
        tree.reset(clone);
    }
    bool lengthes = true;
    std::string bp;
    bp = opt_value("bootstrap-print").as<std::string>();
    TreeNode::ShowBootstrap sbs = parse_bp(bp);
    out << tree->newick(lengthes, sbs) << "\n";
}

const char* GlobalTree::name_impl() const {
    return "Print global tree of genomes";
}

}

