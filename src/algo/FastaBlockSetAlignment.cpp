/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "FastaBlockSetAlignment.hpp"
#include "block_set_alignment.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "write_fasta.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

FastaBlockSetAlignment::FastaBlockSetAlignment():
    file_writer_(this, "bsa-fasta",
                 "Output fasta file with block set alignment") {
    add_opt("bsa-name", "Name of block set alignment.",
            std::string(""), true);
}

struct SeqCmp {
    bool operator()(const Sequence* s1, const Sequence* s2) const {
        return s1->name() < s2->name();
    }
};

void FastaBlockSetAlignment::run_impl() const {
    std::ostream& out = file_writer_.output();
    std::string bsa_name = opt_value("bsa-name").as<std::string>();
    if (!block_set()->has_bsa(bsa_name)) {
        throw Exception("No such bsa: " + bsa_name);
    }
    const BSA& bsa = block_set()->bsa(bsa_name);
    std::map<Sequence*, std::string> seq2text;
    int length = bsa_length(bsa);
    for (int col = 0; col < length; col++) {
        int max_l = 0;
        BOOST_FOREACH (const BSA::value_type& seq_and_bsrow, bsa) {
            const BSRow& bsrow = seq_and_bsrow.second;
            BOOST_ASSERT(bsrow.fragments.size() > col);
            Fragment* f = bsrow.fragments[col];
            int l = f ? f->alignment_length() : 0;
            max_l = std::max(max_l, l);
        }
        BOOST_FOREACH (const BSA::value_type& seq_and_bsrow, bsa) {
            Sequence* seq = seq_and_bsrow.first;
            const BSRow& bsrow = seq_and_bsrow.second;
            Fragment* f = bsrow.fragments[col];
            std::string text;
            if (f) {
                text = f->str();
            }
            text.resize(max_l, '-');
            seq2text[seq] += text;
        }
    }
    std::vector<Sequence*> seqs;
    BOOST_FOREACH (const BSA::value_type& seq_and_bsrow, bsa) {
        Sequence* seq = seq_and_bsrow.first;
        seqs.push_back(seq);
    }
    std::sort(seqs.begin(), seqs.end(), SeqCmp());
    BOOST_FOREACH (Sequence* seq, seqs) {
        write_fasta(out, seq->name(), seq->description(),
                    seq2text[seq], 60);
    }
}

const char* FastaBlockSetAlignment::name_impl() const {
    return "Print block set alignment as fasta";
}

}

