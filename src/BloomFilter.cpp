/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <boost/math/constants/constants.hpp>

#include "BloomFilter.hpp"

namespace bloomrepeats {

const float ln_two = boost::math::constants::ln_two<float>();

BloomFilter::BloomFilter()
{ }

BloomFilter::BloomFilter(size_t members, float error_prob) {
    set_members(members, error_prob);
    set_optimal_hashes(members);
}

void BloomFilter::set_members(size_t members, float error_prob) {
    set_bits(optimal_bits(members, error_prob));
}

size_t BloomFilter::bits() const {
    return bits_.size();
}

void BloomFilter::set_bits(size_t bits) {
    bits_.resize(0);
    bits_.resize(bits);
}

void BloomFilter::set_optimal_hashes(size_t members) {
    set_hashes(optimal_hashes(members, bits()));
}

size_t BloomFilter::hashes() const {
    return hash_mul_.size();
}

void BloomFilter::set_hashes(size_t hashes) {
    hash_mul_.resize(0);
    hash_mul_.resize(hashes);
    std::srand(std::time(NULL));
    for (size_t i = 0; i < hashes; i++) {
        hash_mul_[i] = std::rand();
    }
}

bool BloomFilter::test_and_add(const char* start, size_t length) {
    bool result = true;
    for (size_t hash = 0; hash < hashes(); hash++) {
        size_t index = make_index(hash, start, length);
        if (!bits_[index]) {
            result = false;
        }
        bits_[index] = true;
    }
    return result;
}

bool BloomFilter::test_and_add(const std::string& member) {
    return test_and_add(member.c_str(), member.length());
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

size_t BloomFilter::optimal_bits(size_t members, float error_prob) {
    return round(members * (-std::log(error_prob) / (ln_two * ln_two)));
}

size_t BloomFilter::optimal_hashes(size_t members, size_t bits) {
    return round(ln_two * bits / members);
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
        result *= hash_mul;
        result ^= value;
    }
    return result % bits();
}

}

