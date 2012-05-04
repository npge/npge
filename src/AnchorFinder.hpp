/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ANCHOR_HINDER_HPP_
#define BR_ANCHOR_HINDER_HPP_

#include <vector>
#include <boost/function.hpp>

#include "Sequence.hpp"
#include "global.hpp"

namespace bloomrepeats {

/** Finder of possible short anchors.

*/
class AnchorFinder {
public:
    /** Function, called with an anchor.
    Parameters: sequence, start, length
    */
    typedef boost::function<void(Sequence::Ptr, size_t, size_t)> AnchorHandler;

    /** Default anchor size */
    static const size_t ANCHOR_SIZE = 20;

    /** Default constructor */
    AnchorFinder();

    /** Add sequence */
    void add_sequnce(Sequence::Ptr sequence);

    /** Find possible anchors in added sequence.
    Each found anchor candidate is passed to anchor_handler.
    Same anchor can be passed more than once.

    \note If no anchor handler has been set, this method does nothing.
    */
    void run();

    /** Get anchor handler function */
    const AnchorHandler& anchor_handler() const {
        return anchor_handler_;
    }

    /** Set anchor handler function */
    void set_anchor_handler(const AnchorHandler& anchor_handler) {
        anchor_handler_ = anchor_handler;
    }

    /** Get anchor size */
    size_t anchor_size() const {
        return anchor_size_;
    }

    /** Set anchor size.
    Defaults to ANCHOR_SIZE.
    */
    void set_anchor_size(size_t anchor_size) {
        anchor_size_ = anchor_size;
    }

private:
    AnchorHandler anchor_handler_;
    std::vector<Sequence::Ptr> seqs_;
    size_t anchor_size_;

    void test_and_add(Sequence::Ptr sequence, BloomFilter& filter);
};

}

#endif

