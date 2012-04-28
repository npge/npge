/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <cstdlib>
#include <ctime>

#include "BloomFilter.hpp"

namespace br {

typedef std::size_t size_t;

BloomFilter::BloomFilter(size_t bits, size_t hashes):
    bits_(bits), hash_mul_(hashes) {
    std::srand(std::time(NULL));
    for (size_t i = 0; i < hashes; i++) {
        hash_mul_[i] = std::rand();
    }
}

size_t BloomFilter::bits() const {
    return bits_.size();
}

size_t BloomFilter::hashes() const {
    return hash_mul_.size();
}

void BloomFilter::add(const char* start, size_t length) {
    for (size_t hash = 0; hash < hashes(); hash++) {
        bits_[make_index(hash, start, length)] = true;
    }
}

void BloomFilter::add(const std::string& member) {
    add(member.c_str(), member.length());
}

bool BloomFilter::test(const char* start, size_t length) const {
    for (size_t hash = 0; hash < hashes(); hash++) {
        if (!bits_[make_index(hash, start, length)]) {
            return false;
        }
    }
    return true;
}

bool BloomFilter::test(const std::string& member) const {
    return test(member.c_str(), member.length());
}

static size_t char_to_size(char c) {
    if (c == 'a') {
        return 0;
    } else if (c == 't') {
        return 1;
    } else if (c == 'g') {
        return 2;
    } else { // if (c == 'c') {
        return 3;
    }
}

size_t BloomFilter::make_index(size_t hash, const char* start,
                               size_t length) const {
    size_t hash_mul = hash_mul_[hash];
    size_t result = 1;
    const char* end = start + length;
    for (const char* i = start; i < end; i++) {
        size_t value = char_to_size(*i);
        result ^= value;
        result *= hash_mul;
    }
    return result % bits();
}

}

