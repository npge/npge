/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "Filter.hpp"
#include "SizeLimits.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "goodSlices.hpp"
#include "block_stat.hpp"
#include "boundaries.hpp"
#include "char_to_size.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"

namespace npge {

LiteFilter::LiteFilter() {
    add_lite_size_limits_options(this);
    add_opt("remove-fragments", "Delete individual fragments "
            "instead of whole block", true);
}

typedef std::pair<Block*, Fragment*> BF;

struct LFData : public ThreadData {
    std::vector<BF> fragments_;
    Blocks blocks_;
    int min_fragment_;
    int min_block_;

    LFData(const Processor* p) {
        min_fragment_ = p->opt_value("min-fragment").as<int>();
        min_block_ = p->opt_value("min-block").as<int>();
    }
};

ThreadData* LiteFilter::before_thread_impl() const {
    return new LFData(this);
}

void LiteFilter::process_block_impl(Block* block,
                                    ThreadData* data) const {
    LFData* ld = D_CAST<LFData*>(data);
    if (block->size() < ld->min_block_) {
        ld->blocks_.push_back(block);
        return;
    }
    if (block->alignment_length() < ld->min_fragment_) {
        ld->blocks_.push_back(block);
        return;
    }
}

void LiteFilter::after_thread_impl(ThreadData* data) const {
    LFData* ld = D_CAST<LFData*>(data);
    // remove fragments refore blocks
    // because removing block can cause deletion of fragment
    BOOST_FOREACH (const BF& bf, ld->fragments_) {
        Block* block = bf.first;
        Fragment* fragment = bf.second;
        block->erase(fragment);
    }
    BlockSet& bs = *block_set();
    BOOST_FOREACH (Block* block, ld->blocks_) {
        bs.erase(block);
    }
}

const char* LiteFilter::name_impl() const {
    return "Filter blocks (checks only "
           "frangment length and block size)";
}

struct LengthRequirements {
    int min_fragment_length;
    Decimal min_identity;
    int min_end;
    int frame_length;

    LengthRequirements(const Processor* p) {
        min_fragment_length = p->opt_value("min-fragment").as<int>();
        frame_length = p->opt_value("frame-length").as<int>();
        min_identity = p->opt_value("min-identity").as<Decimal>();
        min_end = p->opt_value("min-end").as<int>();
    }
};

Filter::Filter() {
    add_size_limits_options(this);
    add_opt("find-subblocks", "Find and add good subblocks of bad blocks",
            true);
    add_opt("good-to-other", "Do not remove bad blocks, "
            "but copy good blocks to other blockset",
            false);
    declare_bs("target", "Filtered blockset");
    declare_bs("other", "Target blockset for good blocks "
               "(if --good-to-other)");
}

bool Filter::is_good_fragment(const Fragment* fragment) const {
    return fragment->valid();
}

bool Filter::filter_block(Block* block) const {
    TimeIncrementer ti(this);
    std::vector<Fragment*> block_copy(block->begin(), block->end());
    bool result = false;
    BOOST_FOREACH (Fragment* fragment, block_copy) {
        if (!is_good_fragment(fragment)) {
            block->erase(fragment);
            result = true;
        }
    }
    return result;
}

static int minIdentCount(Decimal min_identity) {
    Decimal min_gc = min_identity * MAX_COLUMN_SCORE;
    int min_good_count = min_gc.to_i();
    if (min_gc.fraction()) {
        min_good_count += 1;
    }
    return min_good_count;
}

static Coordinates goodSubblocks(const Block* block,
        const LengthRequirements& lr) {
    int nrows = block->size();
    Fragments ff(block->begin(), block->end());
    std::vector<std::string> rows(nrows);
    std::vector<const char*> crows(nrows);
    for (int i = 0; i < nrows; i++) {
        Fragment* f = ff[i];
        rows[i] = f->str();
        crows[i] = rows[i].c_str();
    }
    int length = block->alignment_length();
    int min_length = lr.min_fragment_length;
    int frame_length = lr.frame_length;
    int min_identity = minIdentCount(lr.min_identity);
    Scores scores = goodColumns(&crows[0], nrows, length,
            min_identity, min_length);
    return goodSlices(scores,
        frame_length, lr.min_end,
        min_identity, min_length);
}

static bool checkAlignment(const Block* block,
                           const LengthRequirements& lr) {
    int length = block->alignment_length();
    Coordinates slices = goodSubblocks(block, lr);
    return slices.size() == 1 &&
        slices.front() == StartStop(0, length - 1);
}

bool Filter::is_good_block(const Block* block) const {
    TimeIncrementer ti(this);
    int min_length = opt_value("min-fragment").as<int>();
    if (block->alignment_length() < min_length) {
        return false;
    }
    BOOST_FOREACH (Fragment* f, *block) {
        if (!is_good_fragment(f)) {
            return false;
        }
    }
    int min_block_size = opt_value("min-block").as<int>();
    int max_block_size = opt_value("max-block").as<int>();
    if (block->size() < min_block_size) {
        return false;
    }
    if (block->size() > max_block_size && max_block_size != -1) {
        return false;
    }
    AlignmentStat al_stat;
    make_stat(al_stat, block);
    Decimal min_identity = opt_value("min-identity").as<Decimal>();
    if (al_stat.alignment_rows() == block->size()) {
        Decimal identity = block_identity(al_stat);
        if (min_identity > 0.05) {
            LengthRequirements lr(this);
            if (!checkAlignment(block, lr)) {
                return false;
            }
        }
    }
    return true;
}

void Filter::find_good_subblocks(const Block* block,
                                 Blocks& good_subblocks) const {
    TimeIncrementer ti(this);
    int min_block_size = opt_value("min-block").as<int>();
    if (block->size() < min_block_size) {
        return;
    }
    const int length = block->alignment_length();
    BOOST_FOREACH (Fragment* fragment, *block) {
        if (!fragment->row()) {
            return;
        }
    }
    LengthRequirements lr(this);
    int min_length = lr.min_fragment_length;
    if (length < min_length) {
        return;
    }
    Coordinates slices = goodSubblocks(block, lr);
    BOOST_FOREACH (const StartStop& slice, slices) {
        Block* gb = block->slice(slice.first, slice.second);
        ASSERT_TRUE(is_good_block(gb));
        good_subblocks.push_back(gb);
    }
}

class FilterData : public ThreadData {
public:
    std::vector<Block*> blocks_to_erase;
    std::vector<Block*> blocks_to_insert;
};

ThreadData* Filter::before_thread_impl() const {
    return new FilterData;
}

void Filter::process_block_impl(Block* block, ThreadData* d) const {
    FilterData* data = boost::polymorphic_downcast<FilterData*>(d);
    bool g_t_o = opt_value("good-to-other").as<bool>();
    bool good = is_good_block(block);
    if (g_t_o && good) {
        data->blocks_to_insert.push_back(block->clone());
    }
    if (!g_t_o && !good) {
        bool find_subblocks = opt_value("find-subblocks").as<bool>();
        std::vector<Block*> subblocks;
        if (find_subblocks) {
            find_good_subblocks(block, subblocks);
        }
        if (!subblocks.empty()) {
            data->blocks_to_erase.push_back(block);
            BOOST_FOREACH (Block* subblock, subblocks) {
                ASSERT_TRUE(is_good_block(subblock));
                data->blocks_to_insert.push_back(subblock);
            }
            return;
        }
        if (filter_block(block)) {
            // some fragments were removed
            if (is_good_block(block)) {
                return;
            }
            subblocks.clear(); // useless
            if (find_subblocks) {
                find_good_subblocks(block, subblocks);
            }
            data->blocks_to_erase.push_back(block);
            BOOST_FOREACH (Block* subblock, subblocks) {
                ASSERT_TRUE(is_good_block(subblock));
                data->blocks_to_insert.push_back(subblock);
            }
            return;
        }
        data->blocks_to_erase.push_back(block);
    }
}

void Filter::after_thread_impl(ThreadData* d) const {
    FilterData* data = boost::polymorphic_downcast<FilterData*>(d);
    BlockSet& target = *block_set();
    BlockSet& o = *other();
    bool g_t_o = opt_value("good-to-other").as<bool>();
    BlockSet& bs_to_insert = g_t_o ? o : target;
    BOOST_FOREACH (Block* block, data->blocks_to_erase) {
        // blocks_to_erase is empty if g_t_o
        target.erase(block);
    }
    BOOST_FOREACH (Block* block, data->blocks_to_insert) {
        bs_to_insert.insert(block);
    }
}

const char* Filter::name_impl() const {
    return "Filter blocks";
}

}

