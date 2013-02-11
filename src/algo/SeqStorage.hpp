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
    SeqStorage(const std::string& storage = "asis");

    /** Get storage mode.
     - "asis": InMemorySequence
     - "compact": CompactSequence
    */
    const std::string& storage() const {
        return storage_;
    }

    /** Set storage mode */
    void set_storage(const std::string& storage) {
        storage_ = storage;
    }

    /** Create new sequence */
    SequencePtr create_sequence() const;

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

private:
    std::string storage_;
};

}

#endif

