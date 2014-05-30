/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>

#include "FindBSA.hpp"
#include "bsa_algo.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "block_hash.hpp"
#include "tree.hpp"
#include "Connector.hpp"
#include "global.hpp"

namespace npge {

typedef Strings SeqGroups;

FindBSA::FindBSA() {
    add_opt("bsa-name", "Name of new block set alignment",
            std::string(""));
    SeqGroups seq_groups;
    add_opt("bsa-seqs", "List of sequences used for alignment, "
            "sequence groups can be selected by genome name or "
            "chromosome name, 'all' means all sequences of block set",
            seq_groups);
    declare_bs("target", "Target blockset");
}

void FindBSA::run_impl() const {
    std::string name = opt_value("bsa-name").as<std::string>();
    SeqGroups seq_groups = opt_value("bsa-seqs").as<SeqGroups>();
    BSA rows;
    BlockSet& bs = *block_set();
    std::map<std::string, SequencePtr> name2seq;
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        name2seq[seq->name()] = seq;
    }
    BOOST_FOREACH (std::string seq_group, seq_groups) {
        SequencePtr seq = name2seq[seq_group];
        if (seq) {
            rows[seq.get()] = BSRow();
        } else {
            BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
                if (seq->genome() == seq_group ||
                        seq->chromosome() == seq_group ||
                        seq_group == "all") {
                    rows[seq.get()] = BSRow();
                }
            }
        }
    }
    bsa_make_rows(rows, *block_set());
    boost::scoped_ptr<TreeNode> tree((bsa_make_tree(rows)));
    BSA& aln = block_set()->bsa(name);
    int genomes = genomes_number(*block_set());
    bsa_make_aln_by_tree(aln, rows, tree.get(), genomes);
    Connector c;
    c.apply(block_set());
    bsa_orient(aln);
}

const char* FindBSA::name_impl() const {
    return "Build block set alignment";
}

}

