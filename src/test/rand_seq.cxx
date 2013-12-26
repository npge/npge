/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <boost/lexical_cast.hpp>

int main(int argc, char** argv) {
    int length = 1e6;
    if (argc >= 2) {
        length = boost::lexical_cast<int>(argv[1]);
    }
    std::string name = "test";
    if (argc >= 3) {
        name = argv[2];
    }
    const int LINE = 60;
    std::srand(std::time(NULL));
    std::cout << ">" << name << std::endl;
    for (int line = 0; line < length / LINE + 1; line++) {
        for (int i = 0; i < LINE && line * LINE + i < length; i++) {
            std::cout << ("ATGC"[std::rand() % 4]);
        }
        std::cout << std::endl;
    }
}

