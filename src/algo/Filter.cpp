/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "Filter.hpp"
#include "Union.hpp"
#include "SizeLimits.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "block_stat.hpp"
#include "boundaries.hpp"
#include "char_to_size.hpp"

namespace bloomrepeats {

Filter::Filter(int min_fragment_length, int min_block_size) {
    add_size_limits_options(this);
    set_opt_value("min-fragment", min_fragment_length);
    set_opt_value("min-block", min_block_size);
    add_opt("find-subblocks", "Find and add good subblocks of bad blocks",
            true);
    add_opt("good-to-other", "Do not remove bad blocks, "
            "but copy good blocks to other blockset",
            false);
}

bool Filter::is_good_fragment(const Fragment* fragment) const {
    int min_fragment_length = opt_value("min-fragment").as<int>();
    int max_fragment_length = opt_value("max-fragment").as<int>();
    return fragment->valid() && fragment->length() >= min_fragment_length &&
           (fragment->length() <= max_fragment_length ||
            max_fragment_length == -1);
}

bool Filter::filter_block(Block* block) const {
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
    double min_spreading = opt_value("min-spreading").as<double>();
    double max_spreading = opt_value("max-spreading").as<double>();
    if (al_stat.spreading() < min_spreading) {
        return false;
    }
    if (al_stat.spreading() > max_spreading) {
        return false;
    }
    double min_identity = opt_value("min-identity").as<double>();
    double max_identity = opt_value("max-identity").as<double>();
    double min_gaps = opt_value("min-gaps").as<double>();
    double max_gaps = opt_value("max-gaps").as<double>();
    if (al_stat.alignment_rows() == block->size()) {
        double identity = block_identity(al_stat);
        int gaps = al_stat.ident_gap() + al_stat.noident_gap();
        double gaps_p = double(gaps) / al_stat.total();
        if (identity < min_identity || identity > max_identity) {
            return false;
        }
        if (gaps_p < min_gaps || gaps_p > max_gaps) {
            return false;
        }
    }
    return true;
}

struct LengthRequirements {
    int min_fragment_length;
    int max_fragment_length;
    double min_spreading;
    double max_spreading;
    double min_identity;
    double max_identity;
    double min_gaps;
    double max_gaps;

    LengthRequirements(const Processor* p) {
        min_fragment_length = p->opt_value("min-fragment").as<int>();
        max_fragment_length = p->opt_value("max-fragment").as<int>();
        min_spreading = p->opt_value("min-spreading").as<double>();
        max_spreading = p->opt_value("max-spreading").as<double>();
        min_identity = p->opt_value("min-identity").as<double>();
        max_identity = p->opt_value("max-identity").as<double>();
        min_gaps = p->opt_value("min-gaps").as<double>();
        max_gaps = p->opt_value("max-gaps").as<double>();
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
        BOOST_ASSERT(row);
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
    double spreading;
    if (avg_length == 0) {
        spreading = 0;
    } else {
        spreading = double(max_length - min_length) / double(avg_length);
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
        ident_nogap(0), ident_gap(0), noident_nogap(0), noident_gap(0)
    { }

    double identity() const {
        return block_identity(ident_nogap, ident_gap,
                              noident_nogap, noident_gap);
    }

    double gaps() const {
        int gaps = ident_gap + noident_gap;
        int nogaps = ident_nogap + noident_nogap;
        return double(gaps) / double(gaps + nogaps);
    }
};

static void add_column(int col,
                       const std::vector<char>& gap,
                       const std::vector<char>& ident,
                       IdentGapStat& stat) {
    if (gap[col]) {
        if (ident[col]) {
            stat.ident_gap += 1;
        } else {
            stat.noident_gap += 1;
        }
    } else {
        if (ident[col]) {
            stat.ident_nogap += 1;
        } else {
            stat.noident_nogap += 1;
        }
    }
}

static void del_column(int col,
                       const std::vector<char>& gap,
                       const std::vector<char>& ident,
                       IdentGapStat& stat) {
    if (gap[col]) {
        if (ident[col]) {
            stat.ident_gap -= 1;
        } else {
            stat.noident_gap -= 1;
        }
    } else {
        if (ident[col]) {
            stat.ident_nogap -= 1;
        } else {
            stat.noident_nogap -= 1;
        }
    }
}

static bool good_contents(const IdentGapStat& stat,
                          const LengthRequirements& lr) {
    double identity = stat.identity();
    double gaps = stat.gaps();
    return identity <= lr.max_identity && identity >= lr.min_identity &&
           gaps <= lr.max_gaps && gaps >= lr.min_gaps;
}

static bool good_block(const Block* block, int start, int stop,
                       const IdentGapStat& stat,
                       const LengthRequirements& lr) {
    return good_contents(stat, lr) && good_lengths(block, start, stop, lr);
}

void expand_end(const Block* block, int start, int& stop,
                const std::vector<char>& gap,
                const std::vector<char>& ident,
                IdentGapStat& stat,
                std::vector<bool>& used,
                const LengthRequirements& lr) {
    int step = 1;
    int alignment_length = block->alignment_length();
    while (stop < alignment_length - 1) {
        for (int i = 0; i < step; i++) {
            stop += 1;
            add_column(stop, gap, ident, stat);
        }
        bool good = good_block(block, start, stop, stat, lr);
        if (good && !used[stop]) {
            step *= 2;
            step = std::min(step, alignment_length - stop - 1);
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
}

void expand_begin(const Block* block, int& start, int stop,
                  const std::vector<char>& gap,
                  const std::vector<char>& ident,
                  IdentGapStat& stat,
                  std::vector<bool>& used,
                  const LengthRequirements& lr) {
    int step = 1;
    while (start > 0) {
        for (int i = 0; i < step; i++) {
            start -= 1;
            add_column(start, gap, ident, stat);
        }
        bool good = good_block(block, start, stop, stat, lr);
        if (good && !used[start]) {
            step *= 2;
            step = std::min(step, start);
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
}

void Filter::find_good_subblocks(const Block* block,
                                 std::vector<Block*>& good_subblocks) const {
    int min_block_size = opt_value("min-block").as<int>();
    if (block->size() < min_block_size) {
        return;
    }
    int alignment_length = block->alignment_length();
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
        bool ident1, gap1, pure_gap;
        int atgc[LETTERS_NUMBER];
        test_column(block, i, ident1, gap1, pure_gap, atgc);
        ident[i] = ident1;
        gap[i] = gap1;
    }
    int min_test = lr.min_fragment_length;
    int max_test = int(double(min_test) / (1.0 - lr.max_gaps)) + 1;
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
        int steps = alignment_length - stop;
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
        for (int pos = start; pos <= stop; pos++) {
            add_column(pos, gap, ident, stat);
        }
        BOOST_ASSERT(good_block(block, start, stop, stat, lr));
        // expand
        expand_end(block, start, stop, gap, ident, stat, used, lr);
        expand_begin(block, start, stop, gap, ident, stat, used, lr);
        BOOST_ASSERT(good_block(block, start, stop, stat, lr));
        Block* gb = block->slice(start, stop);
        if (is_good_block(gb)) {
            good_subblocks.push_back(block->slice(start, stop));
            for (int pos = start; pos <= stop; pos++) {
                BOOST_ASSERT(!used[pos]);
                used[pos] = true;
            }
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

bool Filter::change_blocks_impl(std::vector<Block*>& blocks) const {
    BOOST_FOREACH (Block* block, blocks) {
        BOOST_FOREACH (Fragment* f, *block) {
            f->disconnect();
        }
    }
}

bool Filter::process_block_impl(Block* block, ThreadData* d) const {
    FilterData* data = boost::polymorphic_downcast<FilterData*>(d);
    bool g_t_o = opt_value("good-to-other").as<bool>();
    bool good = is_good_block(block);
    if (g_t_o && good) {
        data->blocks_to_insert.push_back(Union::clone_block(block));
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
                BOOST_ASSERT(is_good_block(subblock));
                data->blocks_to_insert.push_back(subblock);
            }
            return true; // TODO
        }
        if (filter_block(block)) {
            // some fragments were removed
            if (is_good_block(block)) {
                return true; // TODO
            }
            subblocks.clear(); // useless
            if (find_subblocks) {
                find_good_subblocks(block, subblocks);
            }
            data->blocks_to_erase.push_back(block);
            BOOST_FOREACH (Block* subblock, subblocks) {
                BOOST_ASSERT(is_good_block(subblock));
                data->blocks_to_insert.push_back(subblock);
            }
            return true; // TODO
        }
        data->blocks_to_erase.push_back(block);
    }
    return true; // TODO
}

bool Filter::after_thread_impl(ThreadData* d) const {
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
    return true; // TODO
}

const char* Filter::name_impl() const {
    return "Filter blocks";
}

}

