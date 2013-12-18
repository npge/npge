/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <climits>
#include <vector>
#include <string>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>

#include "block_hash.hpp"
#include "Fragment.hpp"
#include "Block.hpp"

namespace bloomrepeats {

uint32_t block_hash(const Block* block) {
    std::vector<std::string> fragment_ids;
    BOOST_FOREACH (Fragment* f, *block) {
        fragment_ids.push_back(f->id());
    }
    std::sort(fragment_ids.begin(), fragment_ids.end());
    std::string joint = boost::algorithm::join(fragment_ids, " ");
    const int LOOP_SIZE = sizeof(uint32_t) * 2; // 2 = for * and for ^
    int new_size = ((joint.size() + LOOP_SIZE - 1) / LOOP_SIZE) * LOOP_SIZE;
    joint.resize(new_size, ' ');
    const uint32_t* value = reinterpret_cast<const uint32_t*>(joint.c_str());
    int loops = joint.size() / LOOP_SIZE;
    uint32_t a = 1;
    for (int i = 0; i < loops; i++) {
        a *= value[2 * i];
        a ^= value[2 * i + 1];
    }
    return a;
}

}

