/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <cstdio>
#include <boost/algorithm/string/join.hpp>

#include "BlastRunner.hpp"
#include "temp_file.hpp"

namespace bloomrepeats {

BlastRunner::BlastRunner():
    evalue_(0.001)
{ }

void BlastRunner::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("in-consensus", po::value<Files>()->required(),
     "Input files with consensuses")
    ("out-hits", po::value<std::string>()->required(),
     "Output file with blast hits")
   ;
}

void BlastRunner::apply_options_impl(const po::variables_map& vm) {
    set_input_files(vm["in-consensus"].as<Files>());
    set_output_file(vm["out-hits"].as<std::string>());
}

bool BlastRunner::run_impl() const {
    std::string input = boost::algorithm::join(input_files(), " ");
    std::string bank = temp_file();
    system(("formatdb -p F -i " + input + " -n " + bank).c_str());
    system(("blastall -p blastn -m 8 -d " + bank + " -i " + input +
            " -e " + boost::lexical_cast<std::string>(evalue()) +
            " > " + output_file()).c_str());
    remove(bank.c_str());
    return true;
}

const char* BlastRunner::name_impl() const {
    return "Blast runner";
}

}

