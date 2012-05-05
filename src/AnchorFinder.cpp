/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "AnchorFinder.hpp"
#include "Sequence.hpp"
#include "BloomFilter.hpp"

namespace bloomrepeats {

AnchorFinder::AnchorFinder():
    anchor_size_(ANCHOR_SIZE)
{ }

void AnchorFinder::add_sequnce(SequencePtr sequence) {
    seqs_.push_back(sequence);
}

void AnchorFinder::run() {
    if (!anchor_handler_) {
        return;
    }
    size_t length_sum = 0;
    BOOST_FOREACH (SequencePtr sequence, seqs_) {
        length_sum += sequence->approximate_size();
    }
    length_sum *= 2; // ori = 1 | -1
    float error_prob = 1.0 / length_sum;
    BloomFilter filter(length_sum, error_prob);
    BOOST_FOREACH (SequencePtr sequence, seqs_) {
        test_and_add(sequence, filter);
    }
}

void AnchorFinder::test_and_add(SequencePtr sequence, BloomFilter& filter) {
    for (size_t start = anchor_size_;; start++) {
        size_t length = anchor_size_;
        const char* data = sequence->get(start, length);
        if (length == anchor_size_) {
            if (filter.test_and_add(data, length, 1)) {
                anchor_handler_(sequence, start, length);
            }
            if (filter.test_and_add(data, length, -1)) {
                anchor_handler_(sequence, start, length);
            }
        } else {
            break;
        }
    }
}

}

