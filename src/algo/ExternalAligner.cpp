/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>

#include "ExternalAligner.hpp"
#include "FastaReader.hpp"
#include "write_fasta.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "cast.hpp"
#include "opts_lib.hpp"

namespace npge {

ExternalAligner::ExternalAligner() {
    add_opt("aligner-cmd",
            "Template of command for external aligner",
            std::string(), true);
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
    if (!go("NPGE_DEBUG").as<bool>()) {
        remove_file(input);
        remove_file(output);
    }
}

void ExternalAligner::align_file(const std::string& input,
                                 const std::string& output) const {
    TimeIncrementer ti(this);
    std::string input_esc = escape_backslash(input);
    std::string output_esc = escape_backslash(output);
    std::string cmd = opt_value("aligner-cmd").as<std::string>();
    std::string cmd_string = str(boost::format(cmd) %
                                 input_esc % output_esc);
    int r = system(cmd_string.c_str());
    if (r) {
        throw Exception("external aligner failed with code " +
                        TO_S(r) + ". Command: " + cmd_string);
    }
}

class AlignmentReader : public FastaReader {
public:
    AlignmentReader(Strings& rows,
                    std::istream& input):
        FastaReader(input),
        rows_(rows), i_(0) {
    }

    void new_sequence(const std::string& name,
                      const std::string& description) {
        i_ = L_CAST<int>(name);
        if (i_ >= rows_.size()) {
            rows_.resize(i_ + 1);
        }
        ASSERT_LT(i_, rows_.size());
    }

    void grow_sequence(const std::string& data) {
        ASSERT_FALSE(rows_.empty());
        ASSERT_LT(i_, rows_.size());
        rows_[i_] += data;
    }

    Strings& rows_;
    int i_;
};

void ExternalAligner::read_alignment(Strings& rows,
                                     const std::string& file) const {
    TimeIncrementer ti(this);
    boost::shared_ptr<std::istream> aligned = name_to_istream(file);
    AlignmentReader reader(rows, *aligned);
    reader.read_all_sequences();
}

std::string ExternalAligner::aligner_type() const {
    return "external";
}

const char* ExternalAligner::name_impl() const {
    return "External aligner";
}

static std::string aligner_cmd(const Processor* p,
                               const std::string& name) {
    return make_external_cmd(p->meta(), name);
}

MafftAligner::MafftAligner() {
    set_opt_getter("aligner-cmd",
                   boost::bind(aligner_cmd, this, "MAFFT"));
    set_opt_prefix("mafft-");
}

std::string MafftAligner::aligner_type() const {
    return "mafft";
}

const char* MafftAligner::name_impl() const {
    return "Mafft aligner";
}

MuscleAligner::MuscleAligner() {
    set_opt_getter("aligner-cmd",
                   boost::bind(aligner_cmd, this, "MUSCLE"));
    set_opt_prefix("muscle-");
}

std::string MuscleAligner::aligner_type() const {
    return "muscle";
}

const char* MuscleAligner::name_impl() const {
    return "Muscle aligner";
}

}

