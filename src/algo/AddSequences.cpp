/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>
#include <boost/foreach.hpp>

#include "AddSequences.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

AddSequences::AddSequences(const std::string& storage):
    SeqStorage(storage)
{ }

void AddSequences::add_options_impl(po::options_description& desc) const {
    SeqStorage::add_options_impl(desc);
    add_unique_options(desc)
    ("in-seqs,i", po::value<Files>()->required(),
     "input fasta file(s)")
   ;
}

void AddSequences::apply_options_impl(const po::variables_map& vm) {
    SeqStorage::apply_options_impl(vm);
    set_input_files(vm["in-seqs"].as<Files>());
}

static void read_all_seqs(const AddSequences* self, std::istream& input,
                          std::vector<SequencePtr>& seqs) {
    while (true) {
        SequencePtr seq = self->create_sequence();
        seq->read_from_file(input);
        if (seq->size() > 0) {
            seqs.push_back(seq);
        } else {
            break;
        }
    }
}

bool AddSequences::run_impl() const {
    std::vector<SequencePtr> seqs;
    BOOST_FOREACH (std::string file_name, input_files()) {
        std::ifstream input_file(file_name.c_str());
        read_all_seqs(this, input_file, seqs);
        // TODO memorize name of input file for each sequence
    }
    block_set()->add_sequences(seqs);
    return !seqs.empty();
}

const char* AddSequences::name_impl() const {
    return "Input sequences";
}

}

