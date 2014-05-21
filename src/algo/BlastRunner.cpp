/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "BlastRunner.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "to_s.hpp"

namespace bloomrepeats {

BlastRunner::BlastRunner():
    file_reader_(this, "in-consensus", "Input files with consensuses"),
    file_writer_(this, "out-hits",
                 "Output file with blast hits", true) {
    add_gopt("blast-plus", "Use blast+ (otherwise blast)",
             "BLAST_PLUS");
    add_gopt("evalue", "Max acceptable e-value of hit",
             "BLAST_EVALUE");
    add_gopt("skip-low-complexity-regions",
             "Tell blast not to search in "
             "low complexity regions", "BLAST_DUST");
}

struct BlastDeleter {
    std::string bank;

    ~BlastDeleter() {
        remove_file(bank);
        remove_file(bank + ".nhr");
        remove_file(bank + ".nin");
        remove_file(bank + ".nsq");
    }
};

const char* FORMATDB = "formatdb -l /dev/null -p F -i {in} -n {bank}";
const char* FORMATDB_PLUS = "makeblastdb -dbtype nucl "
                            "-out {bank} -in {in} -logfile /dev/null";

const char* BLASTN = "blastall -p blastn -m 8 -d {bank} "
                     "-i {in} -e {evalue} -a {workers} -F {F} > {out}";
const char* BLASTN_PLUS = "blastn -task blastn -outfmt 6 -db {bank} "
                          "-query {in} -evalue {evalue} "
                          "-num_threads {workers} -dust {F} > {out}";

void BlastRunner::run_impl() const {
    std::string output_file = file_writer_.output_file();
    ASSERT_MSG(!output_file.empty(),
               "BlastRunner, empty output_file");
    Strings inputs = file_reader_.input_files();
    std::string input = boost::algorithm::join(inputs, " ");
    std::string bank = tmp_file();
    BlastDeleter bd;
    bd.bank = bank;
    bool blast_plus = opt_value("blast-plus").as<bool>();
    std::string cmd1 = blast_plus ? FORMATDB_PLUS : FORMATDB;
    using namespace boost::algorithm;
    replace_first(cmd1, "{in}", input);
    replace_first(cmd1, "{bank}", bank);
    int r = system(cmd1.c_str());
    if (r) {
        throw Exception("formatdb failed with code " + TO_S(r));
    }
    bool slcr = opt_value("skip-low-complexity-regions").as<bool>();
    std::string F = slcr ? "T" : "F";
    if (blast_plus) {
        F = slcr ? "yes" : "no";
    }
    double evalue = opt_value("evalue").as<double>();
    std::string cmd2 = blast_plus ? BLASTN_PLUS : BLASTN;
    replace_first(cmd2, "{bank}", bank);
    replace_first(cmd2, "{in}", input);
    replace_first(cmd2, "{evalue}", TO_S(evalue));
    replace_first(cmd2, "{workers}", TO_S(workers()));
    replace_first(cmd2, "{F}", F);
    replace_first(cmd2, "{out}", output_file);
    r = system(cmd2.c_str());
    if (r) {
        throw Exception("blastall failed with code " + TO_S(r));
    }
}

const char* BlastRunner::name_impl() const {
    return "Blast runner";
}

}

