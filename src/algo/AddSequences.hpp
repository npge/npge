/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ADD_SEQUENCES_HPP_
#define BR_ADD_SEQUENCES_HPP_

#include <vector>

#include "Processor.hpp"

namespace bloomrepeats {

/** Add input sequences to the block set */
class AddSequences : public Processor {
public:
    /** Files list */
    typedef std::vector<std::string> Files;

    /** Constructor */
    AddSequences(const std::string& storage = "asis");

    /** Get files list */
    const std::vector<std::string>& files() const {
        return files_;
    }

    /** Set files list */
    void set_files(const std::vector<std::string>& files) {
        files_ = files;
    }

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

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

private:
    std::vector<std::string> files_;
    std::string storage_;
};

}

#endif

