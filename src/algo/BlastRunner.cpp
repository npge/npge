/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
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

namespace npge {

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
        if (!bank.empty()) {
            remove_file(bank);
            remove_file(bank + ".nhr");
            remove_file(bank + ".nin");
            remove_file(bank + ".nsq");
        }
    }
};

const char* FORMATDB = "formatdb -l {nul} -p F -i {in} "
                       "-n {bank}";
const char* FORMATDB_PLUS = "makeblastdb -dbtype nucl "
                            "-out {bank} -in {in} "
                            "-logfile {nul}";

const char* BLASTN = "blastall -p blastn -m 8 -d {bank} "
                     "-i {in} -e {evalue} -a {workers} -F {F} > {out}";
const char* BLASTN_PLUS = "blastn -task blastn -outfmt 6 -db {bank} "
                          "-query {in} -evalue {evalue} "
                          "-num_threads {workers} -dust {F} > {out}";

static std::string name_in_cmd(const std::string& cmd) {
    int space_in_cmd = cmd.find(' ');
    if (space_in_cmd == std::string::npos) {
        return cmd;
    } else {
        return cmd.substr(0, space_in_cmd);
    }
}

void BlastRunner::run_impl() const {
    std::string output_file = file_writer_.output_file();
    ASSERT_MSG(!output_file.empty(),
               "BlastRunner, empty output_file");
    Strings inputs = file_reader_.input_files();
    std::string input = boost::algorithm::join(inputs, " ");
    std::string bank = tmp_file();
    BlastDeleter bd;
    if (!go("NPGE_DEBUG").as<bool>()) {
        bd.bank = bank;
    }
    bool blast_plus = opt_value("blast-plus").as<bool>();
    std::string cmd1 = blast_plus ? FORMATDB_PLUS : FORMATDB;
    using namespace boost::algorithm;
    replace_first(cmd1, "{in}", escape_backslash(input));
    replace_first(cmd1, "{bank}", escape_backslash(bank));
    replace_first(cmd1, "{nul}", go("DEV_NULL").to_s());
    int r = system(cmd1.c_str());
    if (r) {
        std::string c = name_in_cmd(cmd1);
        throw Exception(c + " failed with code " + TO_S(r) +
                        ". Command: " + cmd1);
    }
    bool slcr = opt_value("skip-low-complexity-regions").as<bool>();
    std::string F = slcr ? "T" : "F";
    if (blast_plus) {
        F = slcr ? "yes" : "no";
    }
    double evalue = opt_value("evalue").as<double>();
    std::string cmd2 = blast_plus ? BLASTN_PLUS : BLASTN;
    replace_first(cmd2, "{bank}", escape_backslash(bank));
    replace_first(cmd2, "{in}", escape_backslash(input));
    replace_first(cmd2, "{evalue}", TO_S(evalue));
    replace_first(cmd2, "{workers}", TO_S(workers()));
    replace_first(cmd2, "{F}", F);
    replace_first(cmd2, "{out}", escape_backslash(output_file));
    r = system(cmd2.c_str());
    if (r) {
        std::string c = name_in_cmd(cmd2);
        throw Exception(c + " failed with code " + TO_S(r) +
                        ". Command: " + cmd1);
    }
}

const char* BlastRunner::name_impl() const {
    return "Blast runner";
}

}

