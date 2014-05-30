/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>

#include "AddBlocks.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "SeqStorage.hpp"
#include "RowStorage.hpp"
#include "read_block_set.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"
#include "to_s.hpp"
#include "global.hpp"

namespace npge {

AddBlocks::AddBlocks():
    file_reader_(this, "in-blocks", "input fasta file(s) with blocks") {
    add_seq_storage_options(this);
    add_row_storage_options(this);
    declare_bs("target", "Default blockset where blocks are added");
}

void AddBlocks::run_impl() const {
    Strings block_sets;
    get_block_sets(block_sets);
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
        BlockSetFastaReader reader(*block_set(), input_file,
                                   row_type(this),
                                   seq_type(this));
        BOOST_FOREACH (const std::string& bs_name, block_sets) {
            reader.set_block_set(bs_name, get_bs(bs_name).get());
        }
        reader.read_all_sequences();
    }
    BOOST_FOREACH (const std::string& bs_name, block_sets) {
        BOOST_FOREACH (const Block* block, *get_bs(bs_name)) {
            test_block(block);
        }
    }
}

const char* AddBlocks::name_impl() const {
    return "Input block set";
}

}

