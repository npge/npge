/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <string>
#include <streambuf>

#include "process.hpp"
#include "meta_pipe.hpp"
#include "Meta.hpp"

using namespace bloomrepeats;

int main(int argc, char** argv) {
    std::string script((std::istreambuf_iterator<char>(std::cin)),
                        std::istreambuf_iterator<char>());
    Meta meta;
    ProcessorPtr p = parse_script(script, &meta);
    return process(argc, argv, p, p->name());
}

