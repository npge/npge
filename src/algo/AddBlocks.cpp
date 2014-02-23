/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
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
#include "throw_assert.hpp"

namespace bloomrepeats {

AddBlocks::AddBlocks():
    file_reader_(this, "in-blocks", "input fasta file(s) with blocks") {
    add_seq_storage_options(this);
    add_row_storage_options(this);
}

bool AddBlocks::run_impl() const {
    int size_before = block_set()->size();
    std::vector<std::string> block_sets;
    get_block_sets(block_sets);
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
        BlockSetFastaReader reader(*block_set(), input_file,
                                   import_alignment(this),
                                   row_type(this),
                                   seq_type(this));
        BOOST_FOREACH (const std::string& bs_name, block_sets) {
            reader.set_block_set(bs_name, get_bs(bs_name).get());
        }
        reader.read_all_sequences();
    }
    if (import_alignment(this)) {
        BOOST_FOREACH (const std::string& bs_name, block_sets) {
            BOOST_FOREACH (Block* block, *get_bs(bs_name)) {
                BOOST_FOREACH (Fragment* f, *block) {
                    BOOST_ASSERT(f->row());
                    BOOST_ASSERT(f->length() <= f->row()->length());
                }
            }
        }
    }
    return block_set()->size() > size_before;
}

const char* AddBlocks::name_impl() const {
    return "Input block set";
}

}

