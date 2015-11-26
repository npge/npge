/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_INFO_HPP_
#define NPGE_INFO_HPP_

#include "Processor.hpp"

namespace npge {

class Stats;

/** Print human readable summary and statistics about blockset (+ subsets) */
class Info : public Processor {
public:
    /** Constructor */
    Info();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    Stats* stats_;

    BlockSetPtr filter_blocks() const;
    void print_seq() const;
    void print_all() const;
    void print_rest() const;
    void print_minor() const;
    void print_hemi() const;
    void print_repeats() const;
    void print_stem() const;
    void print_global() const;
};

}

#endif

