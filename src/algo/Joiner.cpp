/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Joiner.hpp"
#include "MetaAligner.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Connector.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"

namespace npge {

Joiner::Joiner(int max_dist,
               Decimal ratio_to_fragment,
               Decimal gap_ratio) {
    aligner_ = new MetaAligner;
    aligner_->set_parent(this);
    add_opt("join-max-dist",
            "Max allowed distance when joining fragments", max_dist);
    add_opt("join-to-fragment",
            "Max allowed gap length to fragment length ratio "
            "when joining fragments", ratio_to_fragment);
    add_opt("join-to-gap",
            "Max allowed ratio of gaps' lengths (inside a block) "
            "when joining", gap_ratio);
    declare_bs("target", "Target blockset");
}

struct BlockGreater {
    bool operator()(const Block* b1, const Block* b2) const {
        return b1->size() > b2->size();
    }
};

static Block* neighbor_block(Block* b, int ori) {
    Block* result = 0;
    Fragment* f = b->front();
    if (f) {
        Fragment* neighbor_f = (ori == 1) ? f->next() : f->prev();
        if (neighbor_f) {
            result = neighbor_f->block();
        }
    }
    return result;
}

bool Joiner::can_join(Fragment* one, Fragment* another) {
    return one->seq() == another->seq() && one->ori() == another->ori() &&
           one->is_neighbor(*another);
}

int Joiner::can_join(Block* one, Block* another) {
    if (one->weak() || another->weak()) {
        return false;
    }
    if (one->size() != another->size()) {
        return false;
    }
    bool all[3] = {true, false, true};
    for (int ori = 1; ori >= -1; ori -= 2) {
        BOOST_FOREACH (Fragment* f, *one) {
            Fragment* f1 = f->logical_neighbor(ori);
            if (!f1 || f1->block() != another || !Joiner::can_join(f, f1)) {
                all[ori + 1] = false;
                break;
            }
        }
        if (all[ori + 1]) {
            break;
        }
    }
    int result = all[1 + 1] ? 1 : all[-1 + 1] ? -1 : 0;
    ASSERT_FALSE(result && !Block::match(one, another));
    return result;
}

void Joiner::build_alignment(Strings& rows,
                             const Fragments& fragments,
                             const Block* another,
                             int logical_ori) const {
    Strings middle;
    int size = fragments.size();
    middle.resize(size);
    for (int i = 0; i < size; i++) {
        Fragment* f = fragments[i];
        Fragment* f1 = f->logical_neighbor(logical_ori);
        ASSERT_TRUE(f1);
        ASSERT_EQ(f1->block(), another);
        ASSERT_EQ(f1->ori(), f->ori());
        std::string& seq = middle[i];
        int min_pos, max_pos;
        if (f->next() == f1) {
            min_pos = f->max_pos() + 1;
            max_pos = f1->min_pos() - 1;
        } else {
            min_pos = f1->max_pos() + 1;
            max_pos = f->min_pos() - 1;
        }
        if (max_pos >= min_pos) {
            Fragment between(f->seq(), min_pos, max_pos, f->ori());
            seq = between.str(0);
        }
    }
    aligner_->align_seqs(middle);
    rows.resize(size);
    for (int i = 0; i < size; i++) {
        Fragment* f = fragments[i];
        Fragment* f1 = f->logical_neighbor(logical_ori);
        std::string& row = rows[i];
        if (logical_ori == 1) {
            row = f->str() + middle[i] + f1->str();
        } else {
            row = f1->str() + middle[i] + f->str();
        }
    }
}

Block* Joiner::join_blocks(Block* one, Block* another,
                           int logical_ori) const {
    TimeIncrementer ti(this);
    ASSERT_FALSE(one->weak());
    ASSERT_FALSE(another->weak());
    ASSERT_EQ(Joiner::can_join(one, another), logical_ori);
    Block* result = new Block();
    std::set<Fragment*> to_delete;
    Fragments fragments((one->begin()), one->end());
    int size = fragments.size();
    ASSERT_GT(size, 0);
    ASSERT_EQ(another->size(), size);
    Strings rows;
    RowType type;
    bool aln = has_alignment(one) && has_alignment(another);
    if (aln) {
        build_alignment(rows, fragments, another, logical_ori);
        type = one->front()->row()->type();
    }
    Fragments new_fragments;
    BOOST_FOREACH (Fragment* f, fragments) {
        Fragment* f1 = f->logical_neighbor(logical_ori);
        ASSERT_TRUE(f1);
        ASSERT_EQ(f1->block(), another);
        Fragment* new_fragment = join(f, f1);
        result->insert(new_fragment);
        new_fragments.push_back(new_fragment);
        to_delete.insert(f);
        to_delete.insert(f1);
    }
    BOOST_FOREACH (Fragment* f, to_delete) {
        delete f;
    }
    ASSERT_EQ(new_fragments.size(), size);
    if (aln) {
        ASSERT_EQ(rows.size(), size);
        for (int i = 0; i < size; i++) {
            Fragment* new_fragment = new_fragments[i];
            AlignmentRow* new_row = AlignmentRow::new_row(type);
            new_fragment->set_row(new_row);
            new_row->grow(rows[i]);
        }
    }
    return result;
}

Fragment* Joiner::join(Fragment* one, Fragment* another) {
    ASSERT_TRUE(Joiner::can_join(one, another));
    if (another->next() == one) {
        std::swap(one, another);
    }
    ASSERT_EQ(one->next(), another);
    Fragment* new_fragment = new Fragment(one->seq());
    new_fragment->set_min_pos(std::min(one->min_pos(), another->min_pos()));
    new_fragment->set_max_pos(std::max(one->max_pos(), another->max_pos()));
    new_fragment->set_ori(one->ori());
    if (one->prev()) {
        Fragment::connect(one->prev(), new_fragment);
    }
    if (another->next()) {
        Fragment::connect(new_fragment, another->next());
    }
    return new_fragment;
}

bool Joiner::can_join_fragments(Fragment* f1, Fragment* f2) const {
    TimeIncrementer ti(this);
    if (!Joiner::can_join(f1, f2)) {
        return false;
    }
    int dist = f1->dist_to(*f2);
    int min_length = std::min(f1->length(), f2->length());
    ASSERT_GT(min_length, 0);
    Decimal ratio = Decimal(dist) / Decimal(min_length);
    int max_dist = opt_value("join-max-dist").as<int>();
    Decimal to_fragment = opt_value("join-to-fragment").as<Decimal>();
    return (max_dist == -1 || dist <= max_dist) &&
           (to_fragment < 0 || ratio <= to_fragment);
}

bool Joiner::can_join_blocks(Block* b1, Block* b2) const {
    TimeIncrementer ti(this);
    int ori = Joiner::can_join(b1, b2);
    if (ori == 0) {
        return false;
    }
    ASSERT_TRUE(ori);
    ASSERT_FALSE(b1->empty());
    ASSERT_FALSE(b2->empty());
    int min_gap = -1, max_gap = -1;
    BOOST_FOREACH (Fragment* f1, *b1) {
        Fragment* f2 = f1->logical_neighbor(ori);
        ASSERT_TRUE(f2);
        ASSERT_EQ(f2->block(), b2);
        if (!can_join_fragments(f1, f2)) {
            return false;
        }
        int dist = f1->dist_to(*f2);
        min_gap = (min_gap == -1 || dist < min_gap) ? dist : min_gap;
        max_gap = (max_gap == -1 || dist > max_gap) ? dist : max_gap;
    }
    Decimal gap_ratio = opt_value("join-to-gap").as<Decimal>();
    if (gap_ratio >= 0 &&
            Decimal(max_gap) / min_gap > gap_ratio) {
        return false;
    }
    return true;
}

Block* Joiner::try_join(Block* one, Block* another) const {
    TimeIncrementer ti(this);
    Block* result = 0;
    int match_ori = Block::match(one, another);
    if (match_ori == -1) {
        another->inverse();
    }
    if (match_ori) {
        int logical_ori = Joiner::can_join(one, another);
        if (logical_ori && can_join_blocks(one, another)) {
            result = join_blocks(one, another, logical_ori);
        }
    }
    return result;
}

void Joiner::run_impl() const {
    Connector c;
    c.apply(block_set());
    std::vector<Block*> bs(block_set()->begin(), block_set()->end());
    std::sort(bs.begin(), bs.end(), BlockGreater());
    BOOST_FOREACH (Block* block, bs) {
        if (block_set()->has(block)) {
            for (int ori = -1; ori <= 1; ori += 2) {
                while (Block* other_block = neighbor_block(block, ori)) {
                    Block* new_block = try_join(block, other_block);
                    if (new_block) {
                        block_set()->erase(block);
                        block_set()->erase(other_block);
                        block_set()->insert(new_block);
                        block = new_block;
                    } else {
                        break;
                    }
                }
            }
        }
    }
    c.apply(block_set());
}

const char* Joiner::name_impl() const {
    return "Join blocks";
}

}

