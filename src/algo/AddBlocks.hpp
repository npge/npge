/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ADD_BLOCKS_HPP_
#define BR_ADD_BLOCKS_HPP_

#include "Processor.hpp"
#include "OptionsPrefix.hpp"
#include "FileReader.hpp"
#include "SeqStorage.hpp"
#include "RowStorage.hpp"

namespace bloomrepeats {

/** Add blocks and sequences to the block set.
If fragment is located on new sequences, then sequence contents is
revealed from fragment. In this case blocks must cover sequences entirely.

Input sequences can be passed too, before fragments.

Blocks are specified by "block=block_name".

Block set can be specified as "set=123" in fasta description of fragment
or sequence. Default block set is named "target". "set=all" means adding
this block or sequence to all block sets. Use "set=s1,s2,s3" to
specify multiple sets.

See stream >> block_set, stream >> alignment_row.
*/
class AddBlocks : public Processor, public FileReader,
    public OptionsPrefix {
    // FIXME prefix is ugly
public:
    SeqStorage seq_storage_;

    /** Constructor */
    AddBlocks();

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

