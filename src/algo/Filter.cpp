/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
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
    bool rf_;

    LFData(const Processor* p) {
        min_fragment_ = p->opt_value("min-fragment").as<int>();
        min_block_ = p->opt_value("min-block").as<int>();
        rf_ = p->opt_value("remove-fragments").as<bool>();
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
    int bad_fragments = 0;
    BOOST_FOREACH (Fragment* fragment, *block) {
        if (fragment->length() < ld->min_fragment_) {
            if (ld->rf_) {
                ld->fragments_.push_back(BF(block, fragment));
                bad_fragments += 1;
            } else {
                ld->blocks_.push_back(block);
                return;
            }
        }
    }
    if (block->size() - bad_fragments < ld->min_block_) {
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
    int max_fragment_length;
    Decimal min_identity;
    Decimal max_identity;

    LengthRequirements(const Processor* p) {
        min_fragment_length = p->opt_value("min-fragment").as<int>();
        max_fragment_length = p->opt_value("max-fragment").as<int>();
        min_identity = p->opt_value("min-identity").as<Decimal>();
        max_identity = p->opt_value("max-identity").as<Decimal>();
    }
};

struct IdentGapStat {
    int ident_nogap;
    int ident_gap;
    int noident_nogap;
    int noident_gap;

    IdentGapStat():
        ident_nogap(0), ident_gap(0), noident_nogap(0), noident_gap(0) {
    }

    Decimal strict_identity() const {
        return strict_block_identity(ident_nogap, ident_gap,
                                     noident_nogap, noident_gap);
    }
};

static void add_column(bool gap, bool ident, IdentGapStat& stat) {
    if (gap) {
        if (ident) {
            stat.ident_gap += 1;
        } else {
            stat.noident_gap += 1;
        }
    } else {
        if (ident) {
            stat.ident_nogap += 1;
        } else {
            stat.noident_nogap += 1;
        }
    }
}

static void del_column(bool gap, bool ident, IdentGapStat& stat) {
    if (gap) {
        if (ident) {
            stat.ident_gap -= 1;
        } else {
            stat.noident_gap -= 1;
        }
    } else {
        if (ident) {
            stat.ident_nogap -= 1;
        } else {
            stat.noident_nogap -= 1;
        }
    }
}

static bool strict_good_contents(const IdentGapStat& stat,
                                 const LengthRequirements& lr) {
    return stat.strict_identity() >= lr.min_identity;
}

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
    int min_fragment_length = opt_value("min-fragment").as<int>();
    int max_fragment_length = opt_value("max-fragment").as<int>();
    return fragment->valid() && fragment->length() >= min_fragment_length &&
           (fragment->length() <= max_fragment_length ||
            max_fragment_length == -1);
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

static bool checkAlignment(const Block* block,
                           const LengthRequirements& lr) {
    int length = block->alignment_length();
    int frame = std::min(lr.min_fragment_length, length);
    // get properties of all positions
    std::vector<bool> idents(length), gaps(length);
    for (int pos = 0; pos < length; pos++) {
        bool ident1, gap1;
        test_column(block, pos, ident1, gap1);
        idents[pos] = ident1;
        gaps[pos] = gap1;
    }
    // test first and last columns
    if (!idents[0] || gaps[0]) {
        return false;
    }
    int last = length - 1;
    if (!idents[last] || gaps[last]) {
        return false;
    }
    // test identity of all subblock of MIN_LENGTH
    IdentGapStat stat;
    for (int pos = 0; pos < frame; pos++) {
        add_column(gaps[pos], idents[pos], stat);
    }
    if (!strict_good_contents(stat, lr)) {
        return false;
    }
    for (int new_pos = frame; new_pos < length; new_pos++) {
        int old_pos = new_pos - frame;
        add_column(gaps[new_pos], idents[new_pos], stat);
        del_column(gaps[old_pos], idents[old_pos], stat);
        if (!strict_good_contents(stat, lr)) {
            return false;
        }
    }
    return true;
}

bool Filter::is_good_block(const Block* block) const {
    TimeIncrementer ti(this);
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
    Decimal max_identity = opt_value("max-identity").as<Decimal>();
    if (al_stat.alignment_rows() == block->size()) {
        Decimal identity = block_identity(al_stat);
        if (identity < min_identity || identity > max_identity) {
            return false;
        }
        if (min_identity > 0.05) {
            LengthRequirements lr(this);
            if (!checkAlignment(block, lr)) {
                return false;
            }
        }
    }
    return true;
}

static void findGoodColumns(std::vector<bool>& good_col,
                            const Block* block) {
    int length = block->alignment_length();
    for (int i = 0; i < length; i++) {
        bool ident1, gap1;
        test_column(block, i, ident1, gap1);
        good_col[i] = ident1 && !gap1;
    }
}

struct GoodIdentity {
    int min_good_count_;

    GoodIdentity(int length, Decimal min_identity) {
        Decimal min_gc = Decimal(length) * min_identity;
        min_good_count_ = min_gc.to_i();
        if (min_gc.fraction()) {
            min_good_count_ += 1;
        }
    }

    bool operator()(int good_count) const {
        return good_count >= min_good_count_;
    }
};

static void findGoodFrames(std::vector<bool>& good_frame,
                           const std::vector<bool>& good_col,
                           int frame, int length,
                           Decimal min_identity) {
    GoodIdentity good_identity((frame), min_identity);
    int good_count = 0;
    for (int i = 0; i < frame; i++) {
        good_count += good_col[i];
    }
    good_frame[0] = good_identity(good_count);
    for (int new_pos = frame; new_pos < length; new_pos++) {
        int old_pos = new_pos - frame;
        good_count += good_col[new_pos];
        good_count -= good_col[old_pos];
        int start = old_pos + 1;
        good_frame[start] = good_identity(good_count);
    }
}

struct Frame {
    int length, start;

    Frame(int s, int l):
        length(l), start(s) {
    }

    int stop() const {
        return start + length - 1;
    }

    bool operator<(const Frame& other) const {
        return length > other.length ||
            (length == other.length && start > other.start);
    }

    bool overlaps(const Frame& other) const {
        if (other.start <= start && start <= other.stop()) {
            return true;
        }
        if (other.start <= stop() && stop() <= other.stop()) {
            return true;
        }
        return false;
    }

    Frame exclude(const Frame& other) const {
        int start1 = start, stop1 = stop();
        if (other.start <= start && start <= other.stop()) {
            start1 = other.stop() + 1;
        }
        if (other.start <= stop() && stop() <= other.stop()) {
            stop1 = other.start - 1;
        }
        int length1 = stop1 - start1 + 1;
        return Frame(start1, length1);
    }

    bool valid(int block_length, int min_length) const {
        return length >= min_length && start >= 0 &&
            stop() < block_length;
    }
};

typedef std::vector<Frame> Frames;
typedef std::set<Frame> FramesSet;

static void joinFrames(Frames& frames,
                       const std::vector<bool>& good_frame,
                       int length, int frame) {
    for (int i = 0; i < length - frame + 1; i++) {
        if (good_frame[i]) {
            if (i > 0 && good_frame[i - 1]) {
                // increase previous frame
                ASSERT_GT(frames.size(), 0);
                frames.back().length += 1;
            } else {
                // add new frame
                frames.push_back(Frame(i, frame));
            }
        }
    }
}

static void excludeFrame(FramesSet& fs, const Frame& f,
                         int block_length, int min_length) {
    std::vector<FramesSet::iterator> to_remove;
    Frames to_insert;
    for (FramesSet::iterator it = fs.begin();
            it != fs.end(); ++it) {
        const Frame& f1 = *it;
        if (f1.overlaps(f)) {
            to_remove.push_back(it);
            Frame f2 = f1.exclude(f);
            if (f2.valid(block_length, min_length)) {
                to_insert.push_back(f2);
            }
        }
    }
    BOOST_FOREACH (const FramesSet::iterator& it, to_remove) {
        fs.erase(it);
    }
    BOOST_FOREACH (const Frame& f2, to_insert) {
        fs.insert(f2);
    }
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
    int frame = lr.min_fragment_length;
    if (length < frame) {
        return;
    }
    std::vector<bool> good_col(length);
    findGoodColumns(good_col, block);
    std::vector<bool> good_frame(length - frame + 1);
    // fill frame with first min_fragment_length
    findGoodFrames(good_frame, good_col, frame, length,
                   lr.min_identity);
    Frames frames;
    joinFrames(frames, good_frame, length, frame);
    FramesSet fs(frames.begin(), frames.end());
    while (!fs.empty() && fs.begin()->length >= frame) {
        Frame f = *fs.begin();
        fs.erase(fs.begin());
        Block* gb = block->slice(f.start, f.stop());
        if (is_good_block(gb)) {
            good_subblocks.push_back(gb);
            // exclude this frame from other frames
            excludeFrame(fs, f, length, frame);
        } else {
            // max-length? max-identity?
            delete gb;
        }
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

