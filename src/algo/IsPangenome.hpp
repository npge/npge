/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_IS_PANGENOME_HPP_
#define BR_IS_PANGENOME_HPP_

#include "Processor.hpp"
#include "FileWriter.hpp"
#include "SizeLimits.hpp"

namespace bloomrepeats {

/** Print if block set is a good pangenome.
Requirements of a good pangenome:
 - no overlapping blocks.
 - length of any fragment (not from 1-fragment blocks) >= limit.
 - identity of any fragment >= limit.
 - sequences are covered entirely by blocks (including 1-fragment blocks).
 - blast run on consensuses finds no blocks with length >= 2 * limit and
    identity >= limit.
*/
class IsPangenome : public Processor, public FileWriter, public SizeLimits {
protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

