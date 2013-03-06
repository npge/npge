/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FILTER_HPP_
#define BR_FILTER_HPP_

#include "Processor.hpp"
#include "SizeLimits.hpp"

namespace bloomrepeats {

/** Filter out short and invalid fragments.
Fragments are removed (and disconnected).
If block contains too few fragments, it is remved as well
with all its fragments.
*/
class Filter : public Processor, public SizeLimits {
public:
    /** Constructor */
    Filter(int min_fragment_length = 100, int min_block_size = 2);

    /** Process the block (utility method).
    Return if the block was changed.
    */
    bool filter_block(Block* block) const;

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    /** Make filter */
    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

