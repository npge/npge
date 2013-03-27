/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>

#include "process.hpp"
#include "meta_pipe.hpp"
#include "Meta.hpp"
#include "read_file.hpp"

using namespace bloomrepeats;

int main(int argc, char** argv) {
    std::string script = read_stdin();
    Meta meta;
    ProcessorPtr p = parse_script(script, &meta);
    return process(argc, argv, p, p->name());
}

