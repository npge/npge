/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <ctime>
#include <iostream>

int main() {
    std::srand(std::time(NULL));
    std::cout << ">test" << std::endl;
    for (int line = 0; line < 100000; line++) {
        for (int i = 0; i < 60; i++) {
            std::cout << ("atgc"[std::rand() % 4]);
        }
        std::cout << std::endl;
    }
}

