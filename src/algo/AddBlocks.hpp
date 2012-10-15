/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ADD_BLOCKS_HPP_
#define BR_ADD_BLOCKS_HPP_

#include <vector>

#include "Processor.hpp"

namespace bloomrepeats {

/** Add input blocks to the block set.
\note This processor depends on AddSequences.

Wrapper for stream >> block_set.
*/
class AddBlocks : public Processor {
public:
    /** Files list */
    typedef std::vector<std::string> Files;

    /** Get files list */
    const std::vector<std::string>& files() const {
        return files_;
    }

    /** Set files list */
    void set_files(const std::vector<std::string>& files) {
        files_ = files;
    }

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

private:
    std::vector<std::string> files_;
};

}

#endif

