/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <stdint.h>
#include <sstream>
#include <ostream>
#include <algorithm>
#include "boost-xtime.hpp"
#include <boost/lexical_cast.hpp>

#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Block.hpp"
#include "Sequence.hpp"
#include "complement.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"

namespace npge {

const Fragment Fragment::INVALID = Fragment(SequencePtr(), 1, 0);

Fragment::Fragment(Sequence* seq, size_t min_pos, size_t max_pos, int ori):
    seq_(seq), min_pos_(min_pos), max_pos_(max_pos),
    prev_(0), next_(0), block_and_ori_(0), row_(0) {
    set_ori(ori);
}

Fragment::Fragment(SequencePtr seq, size_t min_pos, size_t max_pos, int ori):
    seq_(seq.get()), min_pos_(min_pos), max_pos_(max_pos),
    prev_(0), next_(0), block_and_ori_(0), row_(0) {
    set_ori(ori);
}

Fragment::Fragment(const Fragment& other):
    prev_(0), next_(0), block_and_ori_(0), row_(0) {
    apply_coords(other);
}

Fragment::~Fragment() {
    disconnect();
    if (block_raw_ptr()) {
        Block* b = block_raw_ptr();
        set_block(0);
        b->erase(this);
    }
    set_row(0);
}

Block* Fragment::block() const {
    return block_raw_ptr();
}

Fragment* Fragment::prev() const {
    ASSERT_TRUE(!prev_ || prev_->next_ == this);
    return prev_;
}

Fragment* Fragment::next() const {
    ASSERT_TRUE(!next_ || next_->prev_ == this);
    return next_;
}

Fragment* Fragment::neighbor(int ori) const {
    return (ori == 1) ? next() : prev();
}

Fragment* Fragment::logical_neighbor(int ori) const {
    return neighbor(this->ori() * ori);
}

bool Fragment::is_neighbor(const Fragment& other) const {
    return prev() == &other || next() == &other;
}

Fragment* Fragment::another_neighbor(const Fragment& other) const {
    ASSERT_TRUE(is_neighbor(other));
    return prev() == &other ? next() : prev();
}

const uintptr_t LAST_BIT = 1;

int Fragment::ori() const {
    return (uintptr_t(block_and_ori_) & LAST_BIT) ? 1 : -1;
}

size_t Fragment::length() const {
    return max_pos() - min_pos() + 1;
}

size_t Fragment::alignment_length() const {
    size_t result = row() ? row()->length() : length();
    ASSERT_MSG(result >= length(),
               ("result=" + TO_S(result) +
                " length=" + TO_S(length())).c_str());
    return result;
}

void Fragment::set_ori(int ori, bool inverse_row) {
    ASSERT_TRUE(ori == 1 || ori == -1);
    if (inverse_row && ori == this->ori() * -1 && row()) {
        InversedRow* r = dynamic_cast<InversedRow*>(row());
        if (r) {
            AlignmentRow* s = r->source();
            r->detach_source();
            set_row(s);
        } else {
            set_row(new InversedRow(row()));
        }
    }
    uintptr_t block_and_ori = uintptr_t(block_and_ori_);
    block_and_ori &= ~LAST_BIT;
    block_and_ori |= (ori == 1) ? LAST_BIT : 0;
    block_and_ori_ = (Block*)block_and_ori;
}

size_t Fragment::begin_pos() const {
    return ori() == 1 ? min_pos() : max_pos();
}

void Fragment::set_begin_pos(size_t begin_pos) {
    if (ori() == 1) {
        set_min_pos(begin_pos);
    } else {
        set_max_pos(begin_pos);
    }
}

size_t Fragment::last_pos() const {
    return ori() == 1 ? max_pos() : min_pos();
}

void Fragment::set_last_pos(size_t last_pos) {
    if (ori() == 1) {
        set_max_pos(last_pos);
    } else {
        set_min_pos(last_pos);
    }
}

void Fragment::set_begin_last(size_t begin_pos, size_t last_pos) {
    if (begin_pos <= last_pos) {
        set_min_pos(begin_pos);
        set_max_pos(last_pos);
        set_ori(1);
    } else {
        set_max_pos(begin_pos);
        set_min_pos(last_pos);
        set_ori(-1);
    }
}

size_t Fragment::end_pos() const {
    return ori() == 1 ? max_pos() + 1 : min_pos() - 1;
}

void Fragment::inverse(bool inverse_row) {
    set_ori(ori() == 1 ? -1 : 1, inverse_row);
}

std::string Fragment::str(char gap) const {
    std::stringstream result;
    print_contents(result, gap, /* line */ 0);
    return result.str();
}

std::string Fragment::substr(int min, int max) const {
    min = min < 0 ? length() + min : min;
    max = max < 0 ? length() + max : max;
    int l = max - min + 1;
    size_t seq_pos = frag_to_seq(this, min);
    return seq()->substr(seq_pos, l, ori());
}

Fragment* Fragment::subfragment(size_t from, size_t to) const {
    Fragment* result = new Fragment(*this);
    bool inverse_needed = from > to;
    if (from > to) {
        std::swap(from, to);
    }
    result->set_begin_pos(begin_pos() + from * ori());
    result->set_last_pos(begin_pos() + to * ori());
    if (inverse_needed) {
        result->inverse();
    }
    return result;
}

Fragment* Fragment::clone() const {
    Fragment* f1 = new Fragment(*this);
    if (row()) {
        f1->set_row(row()->clone());
    }
    return f1;
}

std::string Fragment::id() const {
    return seq()->name() + "_" +
           boost::lexical_cast<std::string>(begin_pos())
           + "_" +
           boost::lexical_cast<std::string>(last_pos());
}

hash_t Fragment::hash() const {
    return seq()->hash(begin_pos(), length(), ori());
}

std::string Fragment::seq_name_from_id(const std::string& id) {
    if (id.empty()) {
        return "";
    }
    size_t u1 = id.find('_');
    if (u1 == std::string::npos) {
        return "";
    }
    return id.substr(0, u1);
}

void Fragment::shift_end(int shift) {
    if (ori() == 1) {
        set_max_pos(max_pos() + shift);
    } else { /* if (ori() == -1) */
        set_min_pos(min_pos() - shift);
    }
}

int Fragment::max_shift_end(int max_overlap) const {
    int result = (ori() == 1) ? (seq()->size() - max_pos() - 1) : min_pos();
    if (max_overlap != -1) {
        Fragment* neighbor = logical_neighbor(1);
        if (neighbor) {
            int n_shift = (ori() == 1) ? (neighbor->min_pos() - max_pos() - 1) :
                          (min_pos() - neighbor->max_pos() - 1);
            n_shift += max_overlap;
            result = std::min(result, n_shift);
        }
    }
    return result;
}

bool Fragment::valid() const {
    return min_pos() <= max_pos() && max_pos() < seq()->size();
}

bool Fragment::operator==(const Fragment& other) const {
    return min_pos() == other.min_pos() && max_pos() == other.max_pos() &&
           ori() == other.ori() && seq() == other.seq();
}

bool Fragment::operator!=(const Fragment& other) const {
    return !(*this == other);
}

bool Fragment::operator<(const Fragment& other) const {
    return min_pos() < other.min_pos() ||
           (min_pos() == other.min_pos() &&
            (max_pos() < other.max_pos() ||
             (max_pos() == other.max_pos() &&
              (ori() < other.ori() ||
               (ori() == other.ori() &&
                seq() < other.seq())))));
}

bool Fragment::has(size_t pos) const {
    return min_pos_ <= pos && pos <= max_pos_;
}

char Fragment::raw_at(int pos) const {
    ASSERT_TRUE(seq_);
    char raw = seq_->char_at(begin_pos() + ori() * pos);
    return ori() == 1 ? raw : complement(raw);
}

char Fragment::at(int pos) const {
    return raw_at(pos >= 0 ? pos : length() + pos);
}

char Fragment::alignment_at(int pos) const {
    if (row()) {
        pos = row()->map_to_fragment(pos);
    }
    return (pos >= 0 && pos < length()) ? raw_at(pos) : 0;
}

void Fragment::connect(Fragment* first, Fragment* second) {
    ASSERT_TRUE(first);
    ASSERT_TRUE(second);
    if (first->next_ != second) {
        if (first->next_) {
            first->next_->prev_ = 0;
        }
        if (second->prev_) {
            second->prev_->next_ = 0;
        }
        first->next_ = second;
        second->prev_ = first;
    } else {
        ASSERT_EQ(second->prev_, first);
    }
#ifndef NDEBUG
    first->next();
    first->prev();
    second->next();
    second->prev();
#endif
}

void Fragment::connect(Fragment* first, Fragment* second, int ori) {
    if (ori == 1) {
        connect(first, second);
    } else {
        connect(second, first);
    }
}

void Fragment::rearrange_with(Fragment* other) {
    Fragment* this_prev = prev();
    Fragment* this_next = next();
    Fragment* other_prev = other->prev();
    Fragment* other_next = other->next();
    Fragment* this_ptr = this;
    this->disconnect(/* connect_neighbors */ false);
    other->disconnect(/* connect_neighbors */ false);
    if (this_prev && this_prev != other) {
        connect(this_prev, other);
    }
    if (this_next && this_next != other) {
        connect(other, this_next);
    }
    if (other_prev && other_prev != this_ptr) {
        connect(other_prev, this_ptr);
    }
    if (other_next && other_next != this_ptr) {
        connect(this_ptr, other_next);
    }
    if (this_next == other) {
        connect(other, this_ptr);
    }
    if (other_next == this_ptr) {
        connect(this_ptr, other);
    }
}

void Fragment::find_place() {
    for (int ori = -1; ori <= 1; ori += 2) {
        while (Fragment* n = neighbor(ori)) {
            if ((ori == 1 && *n < *this) || (ori == -1 && *this < *n)) {
                rearrange_with(n);
            } else {
                break;
            }
        }
    }
}

void Fragment::find_place(Fragment* start_from) {
    disconnect();
    if (start_from->next()) {
        Fragment::connect(this, start_from->next());
    }
    Fragment::connect(start_from, this);
    find_place();
}

void Fragment::disconnect(bool connect_neighbors) {
    if (connect_neighbors && next_ && next_ != this &&
            prev_ && prev_ != this) {
        connect(prev(), next());
    } else {
        if (next_) {
            next_->prev_ = 0;
        }
        if (prev_) {
            prev_->next_ = 0;
        }
    }
    next_ = 0;
    prev_ = 0;
}

size_t Fragment::common_positions(const Fragment& other) const {
    size_t result = 0;
    if (seq() == other.seq()) {
        size_t max_min = std::max(min_pos(), other.min_pos());
        size_t min_max = std::min(max_pos(), other.max_pos());
        if (max_min <= min_max) {
            result = min_max - max_min + 1;
        }
    }
    return result;
}

size_t Fragment::dist_to(const Fragment& other) const {
    ASSERT_EQ(seq(), other.seq());
    if (common_positions(other)) {
        return 0;
    } else if (*this < other) {
        return other.min_pos() - max_pos() - 1;
    } else {
        return min_pos() - other.max_pos() - 1;
    }
}

Fragment Fragment::common_fragment(const Fragment& other) const {
    if (seq() == other.seq()) {
        size_t max_min = std::max(min_pos(), other.min_pos());
        size_t min_max = std::min(max_pos(), other.max_pos());
        if (max_min <= min_max) {
            Fragment res(seq(), max_min, min_max, ori());
            ASSERT_EQ(res.length(), common_positions(other));
            return res;
        }
    }
    return INVALID;
}

bool Fragment::is_subfragment_of(const Fragment& other) const {
    bool result = seq() == other.seq() &&
                  min_pos() >= other.min_pos() && max_pos() <= other.max_pos();
    ASSERT_EQ(result, (common_positions(other) == length()));
    return result;
}

bool Fragment::is_internal_subfragment_of(const Fragment& other) const {
    bool result = seq() == other.seq() &&
                  min_pos() > other.min_pos() && max_pos() < other.max_pos();
    ASSERT_TRUE(!result || is_subfragment_of(other));
    return result;
}

Fragment::Diff Fragment::diff_to(const Fragment& other) const {
    ASSERT_EQ(seq(), other.seq());
    Diff diff;
    diff.begin = ori() * (int(other.begin_pos()) - int(begin_pos()));
    diff.last = ori() * (int(other.last_pos()) - int(last_pos()));
    diff.ori = other.ori() == ori() ? 1 : -1;
    return diff;
}

void Fragment::patch(const Fragment::Diff& diff) {
    size_t new_begin = begin_pos() + ori() * diff.begin;
    size_t new_last = last_pos() + ori() * diff.last;
    set_ori(ori() * diff.ori);
    set_begin_pos(new_begin);
    set_last_pos(new_last);
}

void Fragment::apply_coords(const Fragment& other) {
    seq_ = other.seq();
    set_min_pos(other.min_pos());
    set_max_pos(other.max_pos());
    set_ori(other.ori());
}

Fragment& Fragment::operator=(const Fragment& other) {
    if (this != &other) {
        apply_coords(other);
    }
    return *this;
}

void Fragment::exclude(const Fragment& other) {
    ASSERT_EQ(seq(), other.seq());
    size_t max_min = std::max(min_pos(), other.min_pos());
    size_t min_max = std::min(max_pos(), other.max_pos());
    if (max_min <= min_max) {
        if (min_pos() < other.min_pos()) {
            set_max_pos(other.min_pos() - 1);
        } else if (max_pos() > other.max_pos()) {
            set_min_pos(other.max_pos() + 1);
        } else {
            size_t old_min = min_pos();
            set_min_pos(max_pos() + 1); // +1 for fragments of length=1
            set_max_pos(old_min);
            ASSERT_FALSE(valid());
        }
    }
}

Fragment::Diff Fragment::exclusion_diff(const Fragment& other) const {
    Fragment fr;
    fr.apply_coords(*this);
    fr.exclude(other);
    return diff_to(fr);
}

Fragment* Fragment::split(size_t new_length) {
    Fragment* result = 0;
    if (length() > new_length) {
        result = new Fragment();
        result->apply_coords(*this);
        result->set_begin_pos(begin_pos() + ori() * new_length);
        ASSERT_EQ(result->length() + new_length, length());
        set_last_pos(begin_pos() + ori() * (new_length - 1));
        ASSERT_EQ(length(), new_length);
        ASSERT_FALSE(common_positions(*result));
        find_place();
        result->find_place(this);
        ASSERT_TRUE(result->valid());
    }
    ASSERT_TRUE(valid());
    return result;
}

AlignmentRow* Fragment::detach_row() {
    AlignmentRow* result = row_;
    if (row_) {
        row_->set_fragment(0);
        row_ = 0;
    }
    return result;
}

void Fragment::set_row(AlignmentRow* row) {
    if (row_ && row_->fragment() && row != row_) {
        row_->set_fragment(0);
        delete row_;
    }
    row_ = row;
    if (row_) {
        row_->set_fragment(this);
    }
}

void Fragment::print_header(std::ostream& o, const Block* b) const {
    const Block* bb = b ? : block();
    o << id();
    if (bb) {
        if (bb->name().find(' ') == std::string::npos) {
            o << " block=" << bb->name();
        } else {
            o << ' ' << '"' << "block=" << bb->name() << '"';
        }
    }
    if (prev()) {
        o << " prev=" << prev()->id();
    }
    if (next()) {
        o << " next=" << next()->id();
    }
    if (!row()) {
        o << " norow";
    }
}

void Fragment::print_contents(std::ostream& o, char gap, int line) const {
    int l = 0;
    if (row_ && gap) {
        int row_length = row_->length();
        ASSERT_GTE(row_length, length());
        for (int align_pos = 0; align_pos < row_length; align_pos++) {
            if (l >= line && line != 0) {
                o << std::endl;
                l = 0;
            }
            int fragment_pos = row_->map_to_fragment(align_pos);
            if (fragment_pos == -1) {
                o << gap;
            } else {
                o << raw_at(fragment_pos);
            }
            l += 1;
        }
    } else {
        for (size_t i = 0; i < length(); i++) {
            if (l >= line && line != 0) {
                o << std::endl;
                l = 0;
            }
            o << raw_at(i);
            l += 1;
        }
    }
}

void Fragment::set_block(Block* block) {
    ASSERT_FALSE(uintptr_t(block) & LAST_BIT);
    uintptr_t block_and_ori = uintptr_t(block_and_ori_);
    block_and_ori &= LAST_BIT;
    block_and_ori |= uintptr_t(block);
    block_and_ori_ = (Block*)block_and_ori;
}

Block* Fragment::block_raw_ptr() const {
    uintptr_t result = uintptr_t(block_and_ori_);
    result &= ~LAST_BIT;
    return (Block*)result;
}

std::ostream& operator<<(std::ostream& o, const Fragment& f) {
    o << '>';
    f.print_header(o);
    o << std::endl;
    f.print_contents(o);
    o << std::endl;
    return o;
}

}

