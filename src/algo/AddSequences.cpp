/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "AddSequences.hpp"
#include "SeqStorage.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

AddSequences::AddSequences():
        file_reader_(this, "in-seqs", "input fasta file(s)") {
    add_seq_storage_options(this);
}

static void read_all_seqs(const AddSequences* self, std::istream& input,
                          std::vector<SequencePtr>& seqs) {
    while (true) {
        SequencePtr seq = create_sequence(self);
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
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
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

