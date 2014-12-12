/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ADD_BLOCKS_HPP_
#define NPGE_ADD_BLOCKS_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"

namespace npge {

/** Add blocks and sequences to the blockset.
If fragment is located on new sequences, then sequence contents is
revealed from fragment. In this case blocks must cover sequences entirely.

Input sequences can be passed too, before fragments.

Blocks are specified by "block=block_name".

Blockset can be specified as "set=123" in fasta description of fragment
or sequence. Default blockset is named "target". "set=all" means adding
this block or sequence to all blocksets. Use "set=s1,s2,s3" to
specify multiple sets.

See stream >> block_set, stream >> alignment_row.
*/
class Read : public Processor {
public:
    /** Constructor */
    Read();

protected:
    /** Apply the action */
    void run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
};

}

#endif

