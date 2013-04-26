/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>

#include "AddBlocks.hpp"
#include "BlockSet.hpp"
#include "read_block_set.hpp"

namespace bloomrepeats {

AddBlocks::AddBlocks(bool keep_alignment, RowType row_type,
                     SequenceType seq_type):
    row_storage_(keep_alignment, row_type),
    seq_storage_(seq_type)
{ }

void AddBlocks::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("in-blocks", po::value<Files>()->multitoken()->required(),
     "input fasta file(s) with blocks")
   ;
    row_storage_.add_options_impl(desc);
    seq_storage_.add_options_impl(desc);
}

void AddBlocks::apply_options_impl(const po::variables_map& vm) {
    if (vm.count(prefixed("in-blocks"))) {
        set_input_files(vm[prefixed("in-blocks")].as<Files>());
    }
    row_storage_.apply_options_impl(vm);
    seq_storage_.apply_options_impl(vm);
}

bool AddBlocks::run_impl() const {
    int size_before = block_set()->size();
    BOOST_FOREACH (std::istream& input_file, *this) {
        BlockSetFastaReader reader(*block_set(), input_file,
                                   row_storage_.keep_alignment(),
                                   row_storage_.row_type(),
                                   seq_storage_.seq_type());
        std::vector<std::string> block_sets;
        get_block_sets(block_sets);
        BOOST_FOREACH (const std::string& bs_name, block_sets) {
            reader.set_block_set(bs_name, get_bs(bs_name).get());
        }
        reader.read_all_sequences();
    }
    return block_set()->size() > size_before;
}

const char* AddBlocks::name_impl() const {
    return "Input block set";
}

}

