/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "FixEnds.hpp"
#include "global.hpp"
#include "SizeLimits.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "block_stat.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"

namespace npge {

FixEnds::FixEnds() {
    add_size_limits_options(this);
    declare_bs("target", "Target blockset");
}

struct FEData : public ThreadData {
    Blocks removed_;
    Blocks inserted_;
};

ThreadData* FixEnds::before_thread_impl() const {
    return new FEData;
}

struct GoodAlnFinder {
    typedef std::vector<char> Bools;

    Block* block;
    Bools good_col;
    int length;
    int min_fragment;
    Decimal min_identity;
    int good, min_good;
    int sub_frame;
    int start, stop;

    void init_frame() {
        good_col.clear();
        good_col.resize(min_fragment, false);
        min_good = (min_identity * min_fragment).to_i();
        sub_frame = ((D(1.0) - min_identity) *
                     min_fragment).to_i();
        good = 0;
        for (int i = 0; i < min_fragment; i++) {
            bool g = is_ident_nogap(block, i);
            good += g;
            good_col[i] = g;
        }
        start = 0;
        stop = min_fragment - 1;
    }

    void shift() {
        bool g = is_ident_nogap(block, stop);
        char& c = good_col[stop % min_fragment];
        good -= c;
        good += g;
        c = g;
    }

    bool start_is_good() const {
        return good_col[start % min_fragment];
    }

    bool find_first_good_frame() {
        while (true) {
            if (good >= min_good && start_is_good()) {
                return true;
            }
            start += 1;
            stop += 1;
            if (stop >= length) {
                return false;
            }
            shift();
        }
    }

    int find_start() {
        if (length < min_fragment) {
            return length;
        }
        init_frame();
        bool ok = find_first_good_frame();
        if (!ok) {
            return length;
        }
        int best_score = good;
        int best_start = start;
        while (true) {
            start += 1;
            stop += 1;
            if (stop >= length) {
                break;
            }
            shift();
            bool ok = find_first_good_frame();
            if (!ok || start - best_start > sub_frame) {
                break;
            }
            if (good > best_score) {
                best_score = good;
                best_start = start;
            }
        }
        return best_start;
    }
};

void FixEnds::process_block_impl(Block* b,
                                 ThreadData* d) const {
    ASSERT_TRUE(has_alignment(b));
    GoodAlnFinder gaf;
    gaf.block = b;
    gaf.length = b->alignment_length();
    gaf.min_fragment = opt_value("min-fragment").as<int>();
    gaf.min_identity = opt_value("min-identity").as<Decimal>();
    int start_direct = gaf.find_start();
    b->inverse();
    int start_reverse = gaf.find_start();
    b->inverse();
    if (start_direct == 0 && start_reverse == 0) {
        // block is already good
    } else {
        FEData* data = boost::polymorphic_cast<FEData*>(d);
        data->removed_.push_back(b);
        int stop_direct = gaf.length - start_reverse - 1;
        int slice_length = stop_direct - start_direct + 1;
        if (slice_length >= gaf.min_fragment) {
            Block* s = b->slice(start_direct, stop_direct);
            data->inserted_.push_back(s);
        }
    }
}

void FixEnds::after_thread_impl(ThreadData* d) const {
    FEData* data = boost::polymorphic_cast<FEData*>(d);
    BlockSet& t = *block_set();
    BOOST_FOREACH (Block* b, data->removed_) {
        t.erase(b);
    }
    BOOST_FOREACH (Block* b, data->inserted_) {
        t.insert(b);
    }
}

const char* FixEnds::name_impl() const {
    return "Cut bad aligned ends";
}

}

