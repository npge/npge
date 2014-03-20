/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <boost/algorithm/string/join.hpp>

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
    add_opt("evalue", "Max acceptable e-value of hit", 0.001);
    add_opt("skip-low-complexity-regions",
            "Tell blast not to search in low complexity regions",
            false);
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

bool BlastRunner::run_impl() const {
    std::string output_file = file_writer_.output_file();
    BOOST_ASSERT_MSG(!output_file.empty(), "BlastRunner, empty output_file");
    std::string input = boost::algorithm::join(file_reader_.input_files(), " ");
    std::string bank = tmp_file();
    BlastDeleter bd;
    bd.bank = bank;
    bool slcr = opt_value("skip-low-complexity-regions").as<bool>();
    std::string F = slcr ? " -F T " : " -F F ";
    int r = system(("formatdb -l /dev/null -p F -i " + input +
                    " -n " + bank).c_str());
    if (r) {
        throw Exception("formatdb failed with code " + TO_S(r));
    }
    double evalue = opt_value("evalue").as<double>();
    r = system(("blastall -p blastn -m 8 -d " + bank + " -i " + input +
                " -e " + boost::lexical_cast<std::string>(evalue) +
                " -a " + boost::lexical_cast<std::string>(workers()) +
                // TODO measure ^^
                F +
                " > " + output_file).c_str());
    if (r) {
        throw Exception("blastall failed with code " + TO_S(r));
    }
    return true;
}

const char* BlastRunner::name_impl() const {
    return "Blast runner";
}

}

