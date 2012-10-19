/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OUTPUT_HPP_
#define BR_OUTPUT_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Print blocks to file or to stdout */
class Output : public Processor {
public:
    /** Constructor */
    Output();

    /** Get output file with all blocks */
    const std::string& file() const {
        return file_;
    }

    /** Set output file with all blocks */
    void set_file(const std::string& file) {
        file_ = file;
    }

    /** Get mask of output files (${block} is replaced with block name) */
    const std::string& mask() const {
        return mask_;
    }

    /** Set mask of output files */
    void set_mask(const std::string& mask) {
        mask_ = mask;
    }

    /** Get if alignment will be used if it is available.
    Defaults to true.
    */
    bool export_alignment() const {
        return export_alignment_;
    }

    /** Set if alignment will be used if it is available */
    void set_export_alignment(bool export_alignment) {
        export_alignment_ = export_alignment;
    }

protected:
    /** Add options of all added processors */
    void add_options_impl(po::options_description& desc) const;

    /** Add options to all added processors */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

    const char* name_impl() const;

private:
    std::string file_;
    std::string mask_;
    bool export_alignment_;
};

}

#endif

