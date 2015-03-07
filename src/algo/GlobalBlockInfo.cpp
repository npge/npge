/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <map>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "GlobalBlockInfo.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "block_stat.hpp"
#include "FragmentCollection.hpp"
#include "boundaries.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

GlobalBlockInfo::GlobalBlockInfo(const std::string& prefix) {
    set_opt_prefix(prefix);
    declare_bs("global", "Global blocks");
    declare_bs("normal", "Normal block");
    set_block_set_name("global");
}

void GlobalBlockInfo::print_header(std::ostream& o) const {
    o << "block";
    o << '\t' << "fragments";
    o << '\t' << "min_length";
    o << '\t' << "median_length";
    o << '\t' << "avg_length";
    o << '\t' << "max_length";
    o << '\t' << "s_blocks";
    o << '\t' << "s_blocks_length";
    o << '\t' << "GC";
    o << "\n";
}

void GlobalBlockInfo::prepare() const {
    BlockSetPtr global_bs = get_bs("global");
    BlockSetPtr normal_bs = get_bs("normal");
    global2normal_.clear();
    VectorFc fc;
    fc.add_bs(*global_bs);
    fc.prepare();
    BOOST_FOREACH (Block* b, *normal_bs) {
        BOOST_FOREACH (Fragment* f, *b) {
            Fragments global_ff;
            fc.find_overlap_fragments(global_ff, f);
            BOOST_FOREACH (Fragment* f1, global_ff) {
                global2normal_[f1->block()].insert(b);
            }
        }
    }
}

void GlobalBlockInfo::print_block(std::ostream& o,
                                  Block* block) const {
    AlignmentStat stat;
    make_stat(stat, block);
    typedef Boundaries Ints;
    Ints lengths;
    BOOST_FOREACH (Fragment* f, *block) {
        lengths.push_back(f->length());
    }
    double min = 0, median = 0, avg = 0, max = 0;
    if (!lengths.empty()) {
        min = *std::min_element(lengths.begin(), lengths.end());
        max = *std::max_element(lengths.begin(), lengths.end());
        median = median_element(lengths);
        avg = avg_element_double(lengths);
    }
    //
    G2N::const_iterator it = global2normal_.find(block);
    ASSERT_TRUE(it != global2normal_.end());
    const BlocksSet& normal_blocks = it->second;
    int s_blocks = 0, s_blocks_length = 0;
    BOOST_FOREACH (Block* nb, normal_blocks) {
        using namespace boost::algorithm;
        if (starts_with(nb->name(), "s")) {
            s_blocks += 1;
            s_blocks_length += nb->alignment_length();
        }
    }
    //
    o << block->name();
    o << '\t' << block->size();
    o << '\t' << min;
    o << '\t' << median;
    o << '\t' << avg;
    o << '\t' << max;
    o << '\t' << s_blocks;
    o << '\t' << s_blocks_length;
    o << '\t' << stat.gc();
    o << "\n";
}

const char* GlobalBlockInfo::name_impl() const {
    return "Print information about global blocks";
}

}
