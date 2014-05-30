/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONSENSUS_TREE_HPP_
#define BR_CONSENSUS_TREE_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "FileWriter.hpp"

namespace npge {

class BranchGenerator;

/** Print consensus tree.

PrintTree is used to build tree.
*/
class ConsensusTree : public Processor {
public:
    /** Constructor */
    ConsensusTree();

protected:
    void run_impl() const;
    const char* name_impl() const;

private:
    BranchGenerator* branch_generator_;
    FileWriter file_writer_;
};

}

#endif

