/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_CONSENSUS_TREE_HPP_
#define NPGE_CONSENSUS_TREE_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "FileWriter.hpp"
#include "tree.hpp"

namespace npge {

class BranchGenerator;

bool check_bootstrap_values(std::string& message,
                            Processor* p);

bool check_bootstrap_print(std::string& message,
                           Processor* p);

enum BootstrapValues {
    BLOCKS,
    LENGTH,
    DIAGNOSTIC_POSITIONS
};

BootstrapValues parse_bv(const std::string& bv);

TreeNode::ShowBootstrap parse_bp(const std::string& bp);

void add_diagnostic(TreeNode* tree, BlockSetPtr bs,
                    int min_block, int workers);

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
    FileWriter branch_;
    FileWriter tre_;
};

}

#endif

