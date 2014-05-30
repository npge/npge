/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>

#include "stick_impl.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "SortedVector.hpp"
#include "Graph.hpp"
#include "boundaries.hpp"
#include "throw_assert.hpp"

namespace npge {

void bs_to_sb(Seq2Boundaries& sb, const BlockSet& bs) {
    BOOST_FOREACH (Block* block, bs) {
        BOOST_FOREACH (Fragment* fragment, *block) {
            sb[fragment->seq()].push_back(fragment->min_pos());
            size_t max_pos = fragment->max_pos();
            if (max_pos > 0) {
                max_pos += 1;
            }
            sb[fragment->seq()].push_back(max_pos);
        }
    }
}

bool sb_equal(const Seq2Boundaries& x, const Seq2Boundaries& y) {
    if (x.size() != y.size()) {
        return false;
    }
    BOOST_FOREACH (const Seq2Boundaries::value_type& s_and_b, y) {
        const Boundaries& y_b = s_and_b.second;
        Sequence* seq = s_and_b.first;
        Seq2Boundaries::const_iterator it = x.find(seq);
        if (it == x.end()) {
            return false;
        }
        const Boundaries& x_b = it->second;
        if (x_b.size() != y_b.size()) {
            return false;
        }
        for (int i = 0; i < x_b.size(); ++i) {
            if (x_b[i] != y_b[i]) {
                return false;
            }
        }
    }
    return true;
}

bool sb_match_bs(const Seq2Boundaries& sb, const BlockSet& bs) {
    Seq2Boundaries used_sb;
    bs_to_sb(used_sb, bs);
    BOOST_FOREACH (Seq2Boundaries::value_type& s_and_b, used_sb) {
        Boundaries& b = s_and_b.second;
        b.sort_unique();
    }
    return sb_equal(used_sb, sb);
}

bool stick_fragments(BlockSet& bs, const Seq2Boundaries& sb, int min_distance) {
    bool result = false;
    BOOST_FOREACH (Block* block, bs) {
        bool block_changed = false;
        BOOST_FOREACH (Fragment* f, *block) {
            Seq2Boundaries::const_iterator it = sb.find(f->seq());
            ASSERT_TRUE(it != sb.end());
            const Boundaries& boundaries = it->second;
            int min_pos = nearest_element(boundaries, f->min_pos());
            const int max_dist = 2 * min_distance;
            ASSERT_TRUE(std::abs(int(min_pos - f->min_pos())) < max_dist
                        || min_pos == 0 || min_pos == f->seq()->size() - 1);
            int max_pos = nearest_element(boundaries, f->max_pos() + 1);
            if (max_pos > 0) {
                max_pos -= 1;
            }
            ASSERT_TRUE(std::abs(int(max_pos - f->max_pos())) < max_dist
                        || max_pos == 0 || max_pos == f->seq()->size() - 1);
            if (min_pos != f->min_pos() || max_pos != f->max_pos()) {
                f->set_min_pos(min_pos);
                f->set_max_pos(max_pos);
                block_changed = true;
            }
        }
        if (block_changed) {
            result = true;
            BOOST_FOREACH (Fragment* f, *block) {
                f->set_row(0);
            }
        }
    }
    return result;
}

void stick_boundaries(Seq2Boundaries& sb, int min_distance) {
    BOOST_FOREACH (Seq2Boundaries::value_type& s_and_b, sb) {
        const Sequence* seq = s_and_b.first;
        Boundaries& b = s_and_b.second;
        select_boundaries(b, min_distance, seq->size());
    }
}

void remove_extra_boundaries(Boundaries& x, const Boundaries& y) {
    x.erase(std::remove_if(x.begin(), x.end(),
                           !boost::bind(&Boundaries::has_elem, &y, _1)),
            x.end());
}

void remove_extra_sb(Seq2Boundaries& x, const Seq2Boundaries& y) {
    std::vector<Sequence*> to_remove;
    BOOST_FOREACH (Seq2Boundaries::value_type& s_and_b, x) {
        Sequence* seq = s_and_b.first;
        Boundaries& x_b = s_and_b.second;
        Seq2Boundaries::const_iterator it = y.find(seq);
        if (it == y.end()) {
            to_remove.push_back(seq);
        } else {
            const Boundaries& y_b = it->second;
            remove_extra_boundaries(x_b, y_b);
        }
    }
    BOOST_FOREACH (Sequence* seq, to_remove) {
        x.erase(seq);
    }
}

void remove_extra_sb(Seq2Boundaries& sb, const BlockSet& bs) {
    Seq2Boundaries used_sb;
    bs_to_sb(used_sb, bs);
    BOOST_FOREACH (Seq2Boundaries::value_type& s_and_b, used_sb) {
        Boundaries& b = s_and_b.second;
        b.sort_unique();
    }
    remove_extra_sb(sb, used_sb);
}

std::ostream& operator<<(std::ostream& o, const Seq2Boundaries& sb) {
    BOOST_FOREACH (const Seq2Boundaries::value_type& s_and_b, sb) {
        Sequence* seq = s_and_b.first;
        const Boundaries& b = s_and_b.second;
        o << seq->name() << ": " << b << std::endl;
    }
    return o;
}

}

