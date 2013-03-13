/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SEQ_STORAGE_HPP_
#define BR_SEQ_STORAGE_HPP_

#include "OptionsPrefix.hpp"

namespace bloomrepeats {

/** Utility class storing type of sequence storage */
class SeqStorage : public OptionsPrefix {
public:
    /** Constructor */
    SeqStorage(SequenceType seq_type = COMPACT_SEQUENCE);

    /** Get storage mode */
    SequenceType seq_type() const {
        return seq_type_;
    }

    /** Set storage mode */
    void set_seq_type(SequenceType seq_type) {
        seq_type_ = seq_type;
    }

    /** Create new sequence */
    SequencePtr create_sequence() const;

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

private:
    SequenceType seq_type_;
};

}

#endif

