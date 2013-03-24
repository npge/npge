/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>
#include <boost/foreach.hpp>

#include "AddBlocks.hpp"
#include "Connector.hpp"
#include "Rest.hpp"
#include "BlockSet.hpp"
#include "read_block_set.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

AddBlocks::AddBlocks(bool keep_alignment, RowType row_type,
                     SequenceType seq_type):
    RowStorage(keep_alignment, row_type),
    SeqStorage(seq_type)
{ }

void AddBlocks::add_options_impl(po::options_description& desc) const {
    bloomrepeats::add_unique_options(desc)
    ("in-blocks", po::value<Files>()->required(),
     "input fasta file(s) with blocks")
   ;
    RowStorage::add_options_impl(desc);
    SeqStorage::add_options_impl(desc);
}

void AddBlocks::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("in-blocks")) {
        set_input_files(vm["in-blocks"].as<Files>());
    }
    RowStorage::apply_options_impl(vm);
    SeqStorage::apply_options_impl(vm);
    if (vm.count("in-seqs") && !vm["in-seqs"].as<Files>().empty()) {
        set_seq_type(NO_SEQUENCE);
    }
}

bool AddBlocks::run_impl() const {
    int size_before = block_set()->size();
    BOOST_FOREACH (std::string file_name, input_files()) {
        std::ifstream input_file(file_name.c_str());
        BlockSetFastaReader reader(*block_set(), input_file,
                                   keep_alignment(), row_type(), seq_type());
        reader.read_all_sequences();
    }
#ifndef NDEBUG
    if (seq_type() != NO_SEQUENCE) {
        Connector c;
        c.apply(block_set());
        Rest r(block_set());
        BlockSetPtr rest = new_bs();
        r.apply(rest);
        BOOST_ASSERT_MSG(rest->empty(), "Sequences were not covered entirely "
                "by fragments, please pass --in-seqs");
    }
#endif
    return block_set()->size() > size_before;
}

const char* AddBlocks::name_impl() const {
    return "Input block set";
}

}

