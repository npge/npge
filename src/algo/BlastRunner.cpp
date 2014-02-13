/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <boost/algorithm/string/join.hpp>

#include "BlastRunner.hpp"
#include "temp_file.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

BlastRunner::BlastRunner():
    file_reader_(this, "in-consensus", "Input files with consensuses"),
    file_writer_(this, "out-hits", "Output file with blast hits", true),
    evalue_(0.001), skip_low_complexity_regions_(false)
    // FIXME add options
{ }

bool BlastRunner::run_impl() const {
    std::string output_file = file_writer_.output_file();
    BOOST_ASSERT_MSG(!output_file.empty(), "BlastRunner, empty output_file");
    std::string input = boost::algorithm::join(file_reader_.input_files(), " ");
    std::string bank = temp_file();
    std::string F = skip_low_complexity_regions() ? " -F T " : " -F F ";
    system(("formatdb -l /dev/null -p F -i " + input + " -n " + bank).c_str());
    system(("blastall -p blastn -m 8 -d " + bank + " -i " + input +
            " -e " + boost::lexical_cast<std::string>(evalue()) +
            " -a " + boost::lexical_cast<std::string>(workers()) + // TODO measure
            F +
            " > " + output_file).c_str());
    remove_file(bank);
    remove_file(bank + ".nhr");
    remove_file(bank + ".nin");
    remove_file(bank + ".nsq");
    return true;
}

const char* BlastRunner::name_impl() const {
    return "Blast runner";
}

}

