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
#include "tree.hpp"

namespace bloomrepeats {

BlockSetAlignment::BlockSetAlignment() {
    add_opt("bs-aln-chr", "Chromosome used for block set alignment",
            std::string("chr1"));
}

bool BlockSetAlignment::run_impl() const {
    std::string chr = opt_value("bs-aln-chr").as<std::string>();
    BSA rows;
    bsa_make_rows(rows, *block_set(), chr);
    boost::scoped_ptr<TreeNode> tree((bsa_make_tree(rows)));
    BSA& aln = block_set()->bsa(chr);
    bsa_make_aln_by_tree(aln, rows, tree.get());
    return false;
}

const char* BlockSetAlignment::name_impl() const {
    return "Build block set alignment";
}

}

