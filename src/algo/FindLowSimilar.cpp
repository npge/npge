/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <memory>
#include <algorithm>
#include <boost/cast.hpp>
#include <boost/foreach.hpp>

#include "FindLowSimilar.hpp"
#include "SizeLimits.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "BlockSet.hpp"
#include "block_stat.hpp"
#include "block_hash.hpp"
#include "char_to_size.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

FindLowSimilar::FindLowSimilar() {
    add_size_limits_options(this);
    set_block_set_name("other");
    declare_bs("target", "Destination for low similarity sub-blocks");
    declare_bs("other", "Source blocks");
}

struct FindLowSimilarData : public ThreadData {
    Blocks subblocks_;
};

ThreadData* FindLowSimilar::before_thread_impl() const {
    return new FindLowSimilarData;
}

typedef FindLowSimilar::Region Region;
typedef std::vector<Region> Regions;

int Region::length() const {
    return stop_ - start_ + 1;
}

void Region::set_weight(int weight_factor) {
    if (good_) {
        weight_ = length();
    } else {
        weight_ = length() * weight_factor;
    }
}

int FindLowSimilar::get_weight_factor(Decimal min_identity) {
    return (D(1.0) / (D(1.0) - min_identity)).round();
}

Regions FindLowSimilar::make_regions(const std::vector<bool>& good_col,
                                     int weight_factor) {
    Regions result;
    Region* last = 0;
    for (int i = 0; i < good_col.size(); i++) {
        if (last && last->good_ == good_col[i]) {
            last->stop_ = i;
            last->set_weight(weight_factor);
        } else {
            result.push_back(Region());
            last = &result.back();
            last->start_ = i;
            last->stop_ = i;
            last->good_ = good_col[i];
            last->set_weight(weight_factor);
        }
    }
    return result;
}

int FindLowSimilar::find_min_region(const Regions& regions) {
    int min_index = 0;
    int min_weight = regions[min_index].weight_;
    for (int i = 0; i < regions.size(); i++) {
        if (regions[i].weight_ < min_weight) {
            min_index = i;
            min_weight = regions[min_index].weight_;
        }
    }
    return min_index;
}

Regions FindLowSimilar::merge_region(Regions& regions, int index) {
    Region region = regions[index];
    Region new_region = region;
    if (index > 0) {
        const Region& prev = regions[index - 1];
        new_region.start_ = prev.start_;
        new_region.weight_ += prev.weight_;
        ASSERT_NE(prev.good_, region.good_);
    }
    if (index < regions.size() - 1) {
        const Region& next = regions[index + 1];
        new_region.stop_ = next.stop_;
        new_region.weight_ += next.weight_;
        ASSERT_NE(next.good_, region.good_);
    }
    new_region.good_ = (region.good_ == 0) ? 1 : 0;
    Regions result;
    for (int i = 0; i < regions.size(); i++) {
        if (i == index) {
            result.push_back(new_region);
        } else if (i != index - 1 && i != index + 1) {
            result.push_back(regions[i]);
        }
    }
    return result;
}

void FindLowSimilar::reduce_regions(Regions& regions, int min_length) {
    while (regions.size() >= 2) {
        int min_region_index = find_min_region(regions);
        const Region& min_region = regions[min_region_index];
        if (min_region.weight_ >= min_length) {
            break;
        }
        regions = merge_region(regions, min_region_index);
    }
}

void FindLowSimilar::process_block_impl(Block* block,
                                        ThreadData* data) const {
    int L = block->alignment_length();
    std::vector<bool> good_col((L));
    for (int col = 0; col < L; col++) {
        bool ident, gap, pure_gap;
        int atgc[LETTERS_NUMBER];
        test_column(block, col, ident, gap, pure_gap, atgc);
        good_col[col] = ident && !gap;
    }
    int min_length = opt_value("min-fragment").as<int>();
    Decimal min_identity = opt_value("min-identity").as<Decimal>();
    int weight_factor = get_weight_factor(min_identity);
    Regions regions = make_regions(good_col, weight_factor);
    reduce_regions(regions, min_length);
    FindLowSimilarData* d;
    d = boost::polymorphic_downcast<FindLowSimilarData*>(data);
    Blocks& subblocks = d->subblocks_;
    BOOST_FOREACH (const Region& r, regions) {
        if (!r.good_) {
            typedef std::auto_ptr<Block> BPtr;
            BPtr subblock(block->slice(r.start_, r.stop_));
            Fragments fragments((subblock->begin()),
                                subblock->end());
            BOOST_FOREACH (Fragment* f, fragments) {
                if (f->length() <= 2) {
                    // only gaps
                    subblock->erase(f);
                }
            }
            if (!subblock->empty()) {
                std::string name = block_id(subblock.get());
                subblock->set_name("l" + name);
                subblocks.push_back(subblock.release());
            }
        }
    }
}

void FindLowSimilar::after_thread_impl(ThreadData* data) const {
    FindLowSimilarData* d;
    d = boost::polymorphic_downcast<FindLowSimilarData*>(data);
    BlockSet& target = *block_set();
    BOOST_FOREACH (Block* b, d->subblocks_) {
        target.insert(b);
    }
}

const char* FindLowSimilar::name_impl() const {
    return "Find regions of low similarity in blocks";
}

}

