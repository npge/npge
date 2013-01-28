/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OUTPUT_HPP_
#define BR_OUTPUT_HPP_

#include "AbstractOutput.hpp"

namespace bloomrepeats {

/** Print blocks in fasta format to file or to stdout */
class Output : public AbstractOutput {
public:
    /** Constructor */
    Output();

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

    const char* name_impl() const;

    void print_block(std::ostream& o, Block* block) const;

private:
    bool export_alignment_;
};

}

#endif

