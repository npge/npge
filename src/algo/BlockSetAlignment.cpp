/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>

#include "BlockSetAlignment.hpp"
#include "block_set_alignment.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "tree.hpp"

namespace bloomrepeats {

typedef std::vector<std::string> SeqGroups;

BlockSetAlignment::BlockSetAlignment() {
    add_opt("bsa-name", "Name of new block set alignment",
            std::string(""));
    SeqGroups seq_groups;
    add_opt("bsa-seqs", "List of sequences used for alignment, "
            "sequence groups can be selected by genome name or "
            "chromosome name, 'all' means all sequences of block set",
            seq_groups);
}

bool BlockSetAlignment::run_impl() const {
    std::string name = opt_value("bsa-name").as<std::string>();
    SeqGroups seq_groups = opt_value("bsa-seqs").as<SeqGroups>();
    BSA rows;
    BlockSet& bs = *block_set();
    BOOST_FOREACH (std::string seq_group, seq_groups) {
        SequencePtr seq = bs.seq_from_name(seq_group);
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
    bsa_make_aln_by_tree(aln, rows, tree.get());
    return false;
}

const char* BlockSetAlignment::name_impl() const {
    return "Build block set alignment";
}

}

