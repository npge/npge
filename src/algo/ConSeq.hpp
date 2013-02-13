/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONSEQ_HPP_
#define BR_CONSEQ_HPP_

#include "Processor.hpp"
#include "SeqStorage.hpp"

namespace bloomrepeats {

/** Add consensus sequences, produced from blocks of source block set */
class ConSeq : public Processor, public SeqStorage {
public:
    /** Constructor */
    ConSeq(const BlockSetPtr& source = BlockSetPtr());

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

