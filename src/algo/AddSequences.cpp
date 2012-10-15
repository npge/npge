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
    storage_(storage)
{ }

void AddSequences::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("input-file,i", po::value<Files>()->required(),
     "input fasta file(s)")
    ("seq-storage", po::value<std::string>()->default_value(storage()),
     "way of storing sequence in memory ('asis' or 'compact')");
   ;
}

void AddSequences::apply_options_impl(const po::variables_map& vm) {
    set_files(vm["input-file"].as<Files>());
    std::string storage = vm["seq-storage"].as<std::string>();
    if (storage != "asis" && storage != "compact") {
        throw Exception("'storage' must be 'asis' or 'compact'");
    }
    set_storage(storage);
}

enum StorageMode {
    IN_MEMORY,
    COMPACT
};

static void read_all_seqs(std::istream& input, std::vector<SequencePtr>& seqs,
                          StorageMode mode) {
    while (true) {
        SequencePtr seq;
        if (mode == COMPACT) {
            seq = SequencePtr(new CompactSequence(input));
        } else {
            // use this method by default
            seq = SequencePtr(new InMemorySequence(input));
        }
        if (seq->size() > 0) {
            seqs.push_back(seq);
        } else {
            break;
        }
    }
}

bool AddSequences::run_impl() const {
    std::vector<SequencePtr> seqs;
    StorageMode mode = (storage() == "asis") ? IN_MEMORY : COMPACT;
    BOOST_FOREACH (std::string file_name, files()) {
        std::ifstream input_file(file_name.c_str());
        read_all_seqs(input_file, seqs, mode);
        // TODO memorize name of input file for each sequence
    }
    block_set()->add_sequences(seqs);
    return !seqs.empty();
}

}

