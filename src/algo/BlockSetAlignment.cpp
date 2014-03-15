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
#include "Fragment.hpp"
#include "Block.hpp"
#include "Sequence.hpp"
#include "tree.hpp"

namespace bloomrepeats {

BlockSetAlignment::BlockSetAlignment():
    file_writer_(this, "out-bs-aln",
                 "Output file with block set alignment") {
    add_opt("bs-aln-chr", "Chromosome used for block set alignment",
            std::string("chr1"));
    add_opt("bs-aln-blocks", "Print block names in alignment "
            "(else fragments)", true);
}

bool BlockSetAlignment::run_impl() const {
    std::string chr = opt_value("bs-aln-chr").as<std::string>();
    bool blocks = opt_value("bs-aln-blocks").as<bool>();
    std::ostream& out = file_writer_.output();
    BSA rows, aln;
    bsa_make_rows(rows, *block_set(), chr);
    boost::scoped_ptr<TreeNode> tree((bsa_make_tree(rows)));
    bsa_make_aln_by_tree(aln, rows, tree.get());
    bsa_print(out, aln, blocks);
    return false;
}

}

