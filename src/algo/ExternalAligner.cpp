/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "ExternalAligner.hpp"
#include "FastaReader.hpp"
#include "write_fasta.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "to_s.hpp"

namespace npge {

ExternalAligner::ExternalAligner() {
    add_gopt("aligner-cmd",
             "Template of command for external aligner",
             "EXTERNAL_ALIGNER_CMD");
}

void ExternalAligner::align_seqs_impl(Strings& seqs) const {
    std::string input = tmp_file();
    ASSERT_FALSE(input.empty());
    std::string output = tmp_file();
    ASSERT_FALSE(output.empty());
    {
        boost::shared_ptr<std::ostream> file = name_to_ostream(input);
        std::ostream& out = *file;
        for (int i = 0; i < seqs.size(); i++) {
            write_fasta(out, TO_S(i), "", seqs[i], 60);
        }
    }
    align_file(input, output);
    Strings rows;
    read_alignment(rows, output);
    ASSERT_EQ(rows.size(), seqs.size());
    seqs.swap(rows);
    remove_file(input);
    remove_file(output);
}

void ExternalAligner::align_file(const std::string& input,
                                 const std::string& output) const {
    TimeIncrementer ti(this);
    std::string cmd = opt_value("aligner-cmd").as<std::string>();
    std::string cmd_string = str(boost::format(cmd) % input % output);
    int r = system(cmd_string.c_str());
    if (r) {
        throw Exception("external aligner failed with code " +
                        TO_S(r));
    }
}

class AlignmentReader : public FastaReader {
public:
    AlignmentReader(Strings& rows,
                    std::istream& input):
        rows_(rows),
        FastaReader(input) {
    }

    void new_sequence(const std::string& name,
                      const std::string& description) {
        rows_.push_back("");
    }

    void grow_sequence(const std::string& data) {
        ASSERT_FALSE(rows_.empty());
        rows_.back() += data;
    }

    Strings& rows_;
};

void ExternalAligner::read_alignment(Strings& rows,
                                     const std::string& file) const {
    TimeIncrementer ti(this);
    boost::shared_ptr<std::istream> aligned = name_to_istream(file);
    AlignmentReader reader(rows, *aligned);
    reader.read_all_sequences();
}

const char* ExternalAligner::name_impl() const {
    return "External aligner";
}

MafftAligner::MafftAligner() {
    set_opt_value("aligner-cmd", std::string("$MAFFT_CMD"));
}

const char* MafftAligner::name_impl() const {
    return "Mafft aligner";
}

}

