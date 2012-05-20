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
#include "Fragment.hpp"
#include "make_hash.hpp"

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

bool BloomFilter::test_and_add(const char* start, size_t length, int ori) {
    bool result = true;
    for (size_t hash = 0; hash < hashes(); hash++) {
        size_t index = make_index(hash, start, length, ori);
        if (!bits_[index]) {
            result = false;
        }
        bits_[index] = true;
    }
    return result;
}

bool BloomFilter::test_and_add(const std::string& member, int ori) {
    const char* start = ori == 1 ? member.c_str() :
                        member.c_str() + member.length() - 1;
    return test_and_add(start, member.length(), ori);
}

bool BloomFilter::test_and_add(const Fragment& member) {
    return test_and_add(member.begin(), member.length(), member.ori());
}

void BloomFilter::add(const char* start, size_t length, int ori) {
    for (size_t hash = 0; hash < hashes(); hash++) {
        bits_[make_index(hash, start, length, ori)] = true;
    }
}

void BloomFilter::add(const std::string& member, int ori) {
    const char* start = ori == 1 ? member.c_str() :
                        member.c_str() + member.length() - 1;
    add(start, member.length(), ori);
}

void BloomFilter::add(const Fragment& member) {
    add(member.begin(), member.length(), member.ori());
}

bool BloomFilter::test(const char* start, size_t length, int ori) const {
    for (size_t hash = 0; hash < hashes(); hash++) {
        if (!bits_[make_index(hash, start, length, ori)]) {
            return false;
        }
    }
    return true;
}

bool BloomFilter::test(const std::string& member, int ori) const {
    const char* start = ori == 1 ? member.c_str() :
                        member.c_str() + member.length() - 1;
    return test(start, member.length(), ori);
}

bool BloomFilter::test(const Fragment& member) const {
    return test(member.begin(), member.length(), member.ori());
}

size_t BloomFilter::optimal_bits(size_t members, float error_prob) {
    int result = members * (-std::log(error_prob) / (ln_two * ln_two)) + 0.5;
    if (result % 2 == 0) {
        result += 1;
    }
    return result;
}

size_t BloomFilter::optimal_hashes(size_t members, size_t bits) {
    return round(ln_two * bits / members);
}

size_t BloomFilter::make_index(size_t hash, const char* start,
                               size_t length, int ori) const {
    return make_hash(hash_mul_[hash], start, length, ori) % bits();
}

}

