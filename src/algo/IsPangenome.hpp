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
 - CutGaps and MoveGaps applied to each block.
 - length of any fragment (not from 1-fragment blocks) >= limit.
 - alignment is defined for each block of >= 2 fragments.
 - identity of any block >= limit.
 - sequences are covered entirely by blocks (including 1-fragment blocks).
 - blast run on consensuses finds no blocks with length >= limit and
    identity of consensuses >= limit and
    identity of hits mapped to original blocks >= limit.
    Blast hits are passed through Align, then through Filter.
    Good hits found are saved to blockset "blast-hits".
*/
class IsPangenome : public Processor, public FileWriter, public SizeLimits {
public:
    /** Constructor */
    IsPangenome();

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;

private:
    MoveGaps* move_gaps_;
    CutGaps* cut_gaps_;
};

}

#endif

