/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>

#include "Read.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "SeqStorage.hpp"
#include "RowStorage.hpp"
#include "name_to_stream.hpp"
#include "read_block_set.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"
#include "global.hpp"

namespace npge {

Read::Read():
    file_reader_(this, "in-blocks",
                 "input fasta or blockset file(s)") {
    add_seq_storage_options(this);
    add_row_storage_options(this);
    declare_bs("target",
               "Default blockset where blocks are added");
}

void Read::run_impl() const {
    Strings block_sets;
    get_block_sets(block_sets);
    typedef boost::shared_ptr<std::istream> IStreamPtr;
    std::vector<IStreamPtr> files;
    BOOST_FOREACH (std::string f, file_reader_.input_files()) {
        files.push_back(name_to_istream(f));
    }
    ASSERT_GTE(files.size(), 1);
    BlockSetFastaReader reader(*block_set(), *(files[0]),
                               row_type(this), seq_type(this));
    // add remaining files
    for (int i = 1; i < files.size(); i++) {
        reader.add_input(*(files[i]));
    }
    BOOST_FOREACH (const std::string& bs_name, block_sets) {
        reader.set_block_set(bs_name, get_bs(bs_name).get());
    }
    reader.set_workers(workers());
    reader.run();
    BOOST_FOREACH (const std::string& bs_name, block_sets) {
        BOOST_FOREACH (const Block* block, *get_bs(bs_name)) {
            test_block(block);
        }
    }
}

const char* Read::name_impl() const {
    return "Reads blocks from blockset file, "
           "of sequences from fasta file "
           "(in this case no blocks are added)";
}

}

