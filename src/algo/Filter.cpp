/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
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
#include "to_s.hpp"

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
    Decimal min_spreading;
    Decimal max_spreading;
    Decimal min_identity;
    Decimal max_identity;
    Decimal min_gaps;
    Decimal max_gaps;

    LengthRequirements(const Processor* p) {
        min_fragment_length = p->opt_value("min-fragment").as<int>();
        max_fragment_length = p->opt_value("max-fragment").as<int>();
        min_spreading = p->opt_value("min-spreading").as<Decimal>();
        max_spreading = p->opt_value("max-spreading").as<Decimal>();
        min_identity = p->opt_value("min-identity").as<Decimal>();
        max_identity = p->opt_value("max-identity").as<Decimal>();
        min_gaps = p->opt_value("min-gaps").as<Decimal>();
        max_gaps = p->opt_value("max-gaps").as<Decimal>();
    }
};

// TODO rename Boundaries to smth
typedef Boundaries Integers;

static bool good_lengths(const Block* block, int start, int stop,
                         const LengthRequirements& lr) {
    if (block->empty()) {
        return false;
    }
    Integers lengths;
    BOOST_FOREACH (Fragment* fragment, *block) {
        AlignmentRow* row = fragment->row();
        ASSERT_TRUE(row);
        int f_start = row->nearest_in_fragment(start);
        int f_stop = row->nearest_in_fragment(stop);
        int row_start = row->map_to_alignment(f_start);
        if (row_start < start) {
            f_start += 1;
        }
        int row_stop = row->map_to_alignment(f_stop);
        if (row_stop > stop) {
            f_stop -= 1;
        }
        int f_length = f_stop - f_start + 1;
        if ((lr.max_fragment_length != -1 &&
                f_length > lr.max_fragment_length) ||
                f_length < lr.min_fragment_length) {
            return false;
        }
        lengths.push_back(f_length);
    }
    int max_length = *std::max_element(lengths.begin(), lengths.end());
    int min_length = *std::min_element(lengths.begin(), lengths.end());
    int avg_length = avg_element(lengths);
    Decimal spreading;
    if (avg_length == 0) {
        spreading = 0;
    } else {
        spreading = Decimal(max_length - min_length) / Decimal(avg_length);
    }
    if (spreading > lr.max_spreading || spreading < lr.min_spreading) {
        return false;
    }
    return true;
}

struct IdentGapStat {
    int ident_nogap;
    int ident_gap;
    int noident_nogap;
    int noident_gap;

    IdentGapStat():
        ident_nogap(0), ident_gap(0), noident_nogap(0), noident_gap(0) {
    }

    Decimal identity() const {
        return block_identity(ident_nogap, ident_gap,
                              noident_nogap, noident_gap);
    }

    Decimal strict_identity() const {
        return strict_block_identity(ident_nogap, ident_gap,
                                     noident_nogap, noident_gap);
    }

    Decimal gaps() const {
        int gaps = ident_gap + noident_gap;
        int nogaps = ident_nogap + noident_nogap;
        if (gaps + nogaps > 0) {
            return Decimal(gaps) / Decimal(gaps + nogaps);
        } else {
            return 0.0;
        }
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

static void add_column(int col,
                       const std::vector<char>& gap,
                       const std::vector<char>& ident,
                       IdentGapStat& stat) {
    ASSERT_GTE(col, 0);
    ASSERT_LT(col, gap.size());
    ASSERT_LT(col, ident.size());
    ASSERT_EQ(gap.size(), ident.size());
    add_column(gap[col], ident[col], stat);
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

static void del_column(int col,
                       const std::vector<char>& gap,
                       const std::vector<char>& ident,
                       IdentGapStat& stat) {
    ASSERT_GTE(col, 0);
    ASSERT_LT(col, gap.size());
    ASSERT_LT(col, ident.size());
    ASSERT_EQ(gap.size(), ident.size());
    del_column(gap[col], ident[col], stat);
}

static bool good_contents(const IdentGapStat& stat,
                          const LengthRequirements& lr) {
    Decimal identity = stat.identity();
    Decimal gaps = stat.gaps();
    return identity <= lr.max_identity && identity >= lr.min_identity &&
           gaps <= lr.max_gaps && gaps >= lr.min_gaps;
}

static bool strict_good_contents(const IdentGapStat& stat,
                                 const LengthRequirements& lr) {
    return stat.strict_identity() >= lr.min_identity;
}

static bool good_block(const Block* block, int start, int stop,
                       const IdentGapStat& stat,
                       const LengthRequirements& lr) {
    return good_contents(stat, lr) && good_lengths(block, start, stop, lr);
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
    Decimal min_spreading = opt_value("min-spreading").as<Decimal>();
    Decimal max_spreading = opt_value("max-spreading").as<Decimal>();
    if (al_stat.spreading() < min_spreading) {
        return false;
    }
    if (al_stat.spreading() > max_spreading) {
        return false;
    }
    Decimal min_identity = opt_value("min-identity").as<Decimal>();
    Decimal max_identity = opt_value("max-identity").as<Decimal>();
    Decimal min_gaps = opt_value("min-gaps").as<Decimal>();
    Decimal max_gaps = opt_value("max-gaps").as<Decimal>();
    if (al_stat.alignment_rows() == block->size()) {
        Decimal identity = block_identity(al_stat);
        int gaps = al_stat.ident_gap() + al_stat.noident_gap();
        Decimal gaps_p = Decimal(gaps) / al_stat.total();
        if (identity < min_identity || identity > max_identity) {
            return false;
        }
        if (gaps_p < min_gaps || gaps_p > max_gaps) {
            return false;
        }
        if (min_identity > 0.05) {
            LengthRequirements lr(this);
            int alignment_length = block->alignment_length();
            int frame = std::min(lr.min_fragment_length,
                                 alignment_length);
            bool ident1, gap1;
            IdentGapStat stat_start, stat_stop;
            for (int pos = 0; pos < frame; pos++) {
                test_column(block, pos, ident1, gap1);
                add_column(gap1, ident1, stat_start);
                if (pos == 0) {
                    if (!ident1 || gap1) {
                        return false;
                    }
                }
            }
            if (!strict_good_contents(stat_start, lr)) {
                return false;
            }
            for (int pos = alignment_length - frame;
                    pos < alignment_length; pos++) {
                test_column(block, pos, ident1, gap1);
                add_column(gap1, ident1, stat_stop);
                if (pos == alignment_length - 1) {
                    if (!ident1 || gap1) {
                        return false;
                    }
                }
            }
            if (!strict_good_contents(stat_stop, lr)) {
                return false;
            }
        }
    }
    return true;
}

static void cut_end(const Block* block, int start, int& stop,
                    const std::vector<char>& gap,
                    const std::vector<char>& ident,
                    IdentGapStat& stat,
                    const LengthRequirements& lr) {
    const int alignment_length = block->alignment_length();
    IdentGapStat local_stat;
    int frame = std::min(lr.min_fragment_length, alignment_length);
    int local_start = stop - frame + 1;
    for (int pos = local_start; pos <= stop; pos++) {
        add_column(gap[pos], ident[pos], local_stat);
    }
    while (local_start > start &&
            (gap[stop] || !ident[stop] ||
             !strict_good_contents(local_stat, lr))) {
        del_column(gap[stop], ident[stop], local_stat);
        del_column(gap[stop], ident[stop], stat);
        stop -= 1;
        local_start -= 1;
        add_column(gap[local_start], ident[local_start],
                   local_stat);
    }
    int best_stop = stop;
    int best_score = local_stat.ident_nogap;
    int sub_frame = (Decimal(lr.min_fragment_length) *
                     (D(1.0) - lr.min_identity) * 2).to_i();
    int sub_stop = stop;
    int sub_start = local_start;
    for (int i = 0; i < sub_frame; i++) {
        del_column(gap[sub_stop], ident[sub_stop], local_stat);
        sub_stop -= 1;
        sub_start -= 1;
        if (sub_start < 0) {
            break;
        }
        add_column(gap[sub_start], ident[sub_start], local_stat);
        if (!gap[sub_stop] && ident[sub_stop] &&
                sub_start > start &&
                local_stat.ident_nogap > best_score) {
            best_score = local_stat.ident_nogap;
            best_stop = sub_stop;
        }
    }
    while (stop > best_stop) {
        del_column(gap[stop], ident[stop], stat);
        stop -= 1;
    }
}

static void expand_end(const Block* block, int start, int& stop,
                       const std::vector<char>& gap,
                       const std::vector<char>& ident,
                       IdentGapStat& stat,
                       const std::vector<bool>& used,
                       const LengthRequirements& lr) {
    int step = 1;
    const int alignment_length = block->alignment_length();
    ASSERT_EQ(alignment_length, gap.size());
    ASSERT_EQ(alignment_length, ident.size());
    while (stop < alignment_length - 1) {
        step = std::min(step, alignment_length - stop - 1);
        if (step == 0) {
            break;
        }
        bool was_used = false;
        for (int i = 0; i < step; i++) {
            stop += 1;
            add_column(stop, gap, ident, stat);
            if (used[stop]) {
                was_used = true;
            }
        }
        bool good = good_block(block, start, stop, stat, lr);
        if (good && !was_used) {
            step *= 2;
        } else {
            for (int i = 0; i < step; i++) {
                del_column(stop, gap, ident, stat);
                stop -= 1;
            }
            if (step == 1) {
                break;
            } else {
                step = 1;
            }
        }
    }
    cut_end(block, start, stop, gap, ident, stat, lr);
}

static void cut_begin(const Block* block, int& start, int stop,
                      const std::vector<char>& gap,
                      const std::vector<char>& ident,
                      IdentGapStat& stat,
                      const LengthRequirements& lr) {
    const int alignment_length = block->alignment_length();
    IdentGapStat local_stat;
    int frame = std::min(lr.min_fragment_length, alignment_length);
    int local_stop = start + frame - 1;
    for (int pos = start; pos <= local_stop; pos++) {
        add_column(gap[pos], ident[pos], local_stat);
    }
    while (local_stop < stop &&
            (gap[start] || !ident[start] ||
             !strict_good_contents(local_stat, lr))) {
        del_column(gap[start], ident[start], local_stat);
        del_column(gap[start], ident[start], stat);
        start += 1;
        local_stop += 1;
        add_column(gap[local_stop], ident[local_stop], local_stat);
    }
    int best_start = start;
    int best_score = local_stat.ident_nogap;
    int sub_frame = (Decimal(lr.min_fragment_length) *
                     (D(1.0) - lr.min_identity) * 2).to_i();
    int sub_stop = local_stop;
    int sub_start = start;
    for (int i = 0; i < sub_frame; i++) {
        del_column(gap[sub_start], ident[sub_start], local_stat);
        sub_start += 1;
        sub_stop += 1;
        if (sub_stop >= alignment_length) {
            break;
        }
        add_column(gap[sub_stop], ident[sub_stop], local_stat);
        if (!gap[sub_start] && ident[sub_start] &&
                sub_stop < stop &&
                local_stat.ident_nogap > best_score) {
            best_score = local_stat.ident_nogap;
            best_start = sub_start;
        }
    }
    while (start < best_start) {
        del_column(gap[start], ident[start], stat);
        start += 1;
    }
}

static void expand_begin(const Block* block, int& start, int stop,
                         const std::vector<char>& gap,
                         const std::vector<char>& ident,
                         IdentGapStat& stat,
                         const std::vector<bool>& used,
                         const LengthRequirements& lr) {
    const int alignment_length = block->alignment_length();
    ASSERT_EQ(alignment_length, gap.size());
    ASSERT_EQ(alignment_length, ident.size());
    int step = 1;
    while (start > 0) {
        step = std::min(step, start);
        if (step == 0) {
            break;
        }
        bool was_used = false;
        for (int i = 0; i < step; i++) {
            start -= 1;
            add_column(start, gap, ident, stat);
            if (used[start]) {
                was_used = true;
            }
        }
        bool good = good_block(block, start, stop, stat, lr);
        if (good && !was_used) {
            step *= 2;
        } else {
            for (int i = 0; i < step; i++) {
                del_column(start, gap, ident, stat);
                start += 1;
            }
            if (step == 1) {
                break;
            } else {
                step = 1;
            }
        }
    }
    cut_begin(block, start, stop, gap, ident, stat, lr);
}

void Filter::find_good_subblocks(const Block* block,
                                 Blocks& good_subblocks) const {
    TimeIncrementer ti(this);
    int min_block_size = opt_value("min-block").as<int>();
    if (block->size() < min_block_size) {
        return;
    }
    const int alignment_length = block->alignment_length();
    BOOST_FOREACH (Fragment* fragment, *block) {
        if (!fragment->row()) {
            return;
        }
    }
    LengthRequirements lr(this);
    if (alignment_length < lr.min_fragment_length) {
        return;
    }
    std::vector<char> gap(alignment_length), ident(alignment_length);
    for (int i = 0; i < alignment_length; i++) {
        bool ident1, gap1;
        test_column(block, i, ident1, gap1);
        ident[i] = ident1;
        gap[i] = gap1;
    }
    int min_test = lr.min_fragment_length;
    int max_test = (Decimal(min_test) /
                    (D(1.0) - lr.max_gaps)).to_i() + 1;
    if (max_test > alignment_length) {
        max_test = alignment_length;
    }
    typedef std::pair<int, int> Candidate;
    typedef std::vector<Candidate> Candidates;
    Candidates candidates;
    for (int test = max_test; test >= min_test; test--) {
        int start = 0;
        int stop = start + test - 1;
        IdentGapStat stat;
        for (int pos = start; pos <= stop; pos++) {
            add_column(pos, gap, ident, stat);
        }
        int steps = alignment_length - stop - 1;
        int last_stop = -1;
        for (int i = 0; i < steps; i++) {
            if (start > last_stop &&
                    good_block(block, start, stop, stat, lr)) {
                candidates.push_back(Candidate(start, stop));
                last_stop = stop;
            }
            stop += 1;
            add_column(stop, gap, ident, stat);
            del_column(start, gap, ident, stat);
            start += 1;
        }
    }
    std::vector<bool> used(alignment_length, false);
    BOOST_FOREACH (const Candidate& candidate, candidates) {
        int start = candidate.first;
        int stop = candidate.second;
        if (used[start] || used[stop]) {
            continue;
        }
        IdentGapStat stat;
        bool bad = false;
        for (int pos = start; pos <= stop; pos++) {
            if (used[pos]) {
                bad = true;
                break;
            }
            add_column(pos, gap, ident, stat);
        }
        if (bad) {
            break;
        }
        ASSERT_TRUE(good_block(block, start, stop, stat, lr));
        // expand
        expand_end(block, start, stop, gap, ident, stat, used, lr);
        expand_begin(block, start, stop, gap, ident, stat, used, lr);
        if (good_block(block, start, stop, stat, lr)) {
            Block* gb = block->slice(start, stop);
            if (is_good_block(gb)) {
                good_subblocks.push_back(gb);
                for (int pos = start; pos <= stop; pos++) {
                    ASSERT_FALSE(used[pos]);
                    used[pos] = true;
                }
            } else {
                // max-length? max-identity?
                delete gb;
            }
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

void Filter::change_blocks_impl(std::vector<Block*>& blocks) const {
    BOOST_FOREACH (Block* block, blocks) {
        BOOST_FOREACH (Fragment* f, *block) {
            f->disconnect();
        }
    }
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

