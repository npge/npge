/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
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

namespace npge {

const double ln_two = boost::math::constants::ln_two<double>();

BloomFilter::BloomFilter() {
}

BloomFilter::BloomFilter(size_t members, double error_prob) {
    set_members(members, error_prob);
    set_optimal_hashes(members);
}

void BloomFilter::set_members(size_t members, double error_prob) {
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
    return hash_parameter_.size();
}

void BloomFilter::set_hashes(size_t hashes) {
    hash_parameter_.resize(0);
    hash_parameter_.resize(hashes);
    std::srand(std::time(NULL));
    for (size_t i = 0; i < hashes; i++) {
        hash_parameter_[i] = std::rand();
    }
}

bool BloomFilter::test_and_add(hash_t hash) {
    bool result = true;
    for (size_t i = 0; i < hashes(); i++) {
        size_t index = make_index(i, hash);
        if (!bits_[index]) {
            result = false;
        }
        bits_[index] = true;
    }
    return result;
}

bool BloomFilter::test_and_add(const char* start, size_t length, int ori) {
    return test_and_add(make_hash(start, length, ori));
}

bool BloomFilter::test_and_add(const std::string& member, int ori) {
    const char* start = ori == 1 ? member.c_str() :
                        member.c_str() + member.length() - 1;
    return test_and_add(start, member.length(), ori);
}

bool BloomFilter::test_and_add(const Fragment& member) {
    return test_and_add(member.str());
}

void BloomFilter::add(hash_t hash) {
    for (size_t i = 0; i < hashes(); i++) {
        bits_[make_index(i, hash)] = true;
    }
}

void BloomFilter::add(const char* start, size_t length, int ori) {
    add(make_hash(start, length, ori));
}

void BloomFilter::add(const std::string& member, int ori) {
    const char* start = ori == 1 ? member.c_str() :
                        member.c_str() + member.length() - 1;
    add(start, member.length(), ori);
}

void BloomFilter::add(const Fragment& member) {
    add(member.str());
}

bool BloomFilter::test(hash_t hash) const {
    for (size_t i = 0; i < hashes(); i++) {
        if (!bits_[make_index(i, hash)]) {
            return false;
        }
    }
    return true;
}

bool BloomFilter::test(const char* start, size_t length, int ori) const {
    return test(make_hash(start, length, ori));
}

bool BloomFilter::test(const std::string& member, int ori) const {
    const char* start = ori == 1 ? member.c_str() :
                        member.c_str() + member.length() - 1;
    return test(start, member.length(), ori);
}

bool BloomFilter::test(const Fragment& member) const {
    return test(member.str());
}

size_t BloomFilter::true_bits() const {
    size_t result = 0;
    for (size_t i = 0; i < bits_.size(); i++) {
        if (bits_[i]) {
            result += 1;
        }
    }
    return result;
}

size_t BloomFilter::optimal_bits(size_t members, double error_prob) {
    int result = members * (-std::log(error_prob) / (ln_two * ln_two)) + 0.5;
    if (result % 2 == 0) {
        result += 1;
    }
    if (result < 1) {
        result = 1;
    }
    return result;
}

size_t BloomFilter::optimal_hashes(size_t members, size_t bits) {
    int result = round(ln_two * bits / members);
    if (result < 1) {
        result = 1;
    }
    return result;
}

size_t BloomFilter::make_index(size_t hash_index,
                               hash_t hash) const {
    hash_t xored = (hash ^ hash_parameter_[hash_index]);
    return xored % hash_t(bits());
}

}

