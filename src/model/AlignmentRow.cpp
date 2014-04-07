/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <algorithm>
#include <boost/lexical_cast.hpp>

#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

AlignmentRow::AlignmentRow(Fragment* fragment):
    length_(0), fragment_(0) {
    if (fragment) {
        fragment->set_row(this);
    }
}

AlignmentRow::~AlignmentRow() {
    if (fragment_) {
        Fragment* f = fragment_;
        fragment_ = 0;
        f->set_row(0);
    }
}

void AlignmentRow::grow(const std::string& alignment_string) {
    int align_pos = length();
    int fragment_pos = nearest_in_fragment(align_pos) + 1; // -1 -> 0
    for (int i = 0; i < alignment_string.size(); i++) {
        if (isalpha(alignment_string[i])) {
            BOOST_ASSERT_MSG(!fragment() || fragment()->length() == 1 ||
                             fragment_pos < fragment()->length(),
                             ("Fragment: " + fragment()->id() + "\n" +
                              "Fragment string: " + fragment()->str(/*gap*/ 0) +
                              "\n" +
                              "Alignment string: " + alignment_string + "\n" +
                              "Fragment index: " +
                              boost::lexical_cast<std::string>(fragment_pos) +
                              "\n" +
                              "Alignment index: " +
                              boost::lexical_cast<std::string>(i) + "\n"
                             ).c_str());
            BOOST_ASSERT_MSG(!fragment() ||
                             fragment()->length() == 1 ||
                             tolower(fragment()->raw_at(fragment_pos)) ==
                             tolower(alignment_string[i]),
                             ("Fragment: " + fragment()->id() + "\n" +
                              "Fragment string: " + fragment()->str(/*gap*/ 0) +
                              "\n" +
                              "Alignment string: " + alignment_string + "\n" +
                              "Fragment index: " +
                              boost::lexical_cast<std::string>(fragment_pos) +
                              "\n" +
                              "Alignment index: " +
                              boost::lexical_cast<std::string>(i) + "\n" +
                              "Fragment char: " +
                              fragment()->raw_at(fragment_pos) + "\n" +
                              "Alignment char: " + alignment_string[i]
                             ).c_str());
            bind(fragment_pos, align_pos);
            fragment_pos += 1;
        }
        align_pos += 1;
    }
    set_length(length() + alignment_string.length());
}

int AlignmentRow::nearest_in_fragment(int align_pos) const {
    // FIXME do smth with this
    for (int distance = 0; distance <= length(); distance++) {
        for (int ori = -1; ori <= 1; ori += 2) {
            int new_align_pos = align_pos + ori * distance;
            if (map_to_fragment(new_align_pos) != -1) {
                return map_to_fragment(new_align_pos);
            }
        }
    }
    return -1;
}

void AlignmentRow::assign(const AlignmentRow& other, int start, int stop) {
    clear();
    int length = (stop == -1) ? (other.length() - start) : (stop - start + 1);
    int align_pos = start;
    for (int align_pos = start; align_pos < start + length; align_pos++) {
        int fragment_pos = other.map_to_fragment(align_pos);
        if (fragment_pos != -1) {
            bind(fragment_pos, align_pos);
        }
    }
    set_length(length);
}

AlignmentRow* AlignmentRow::new_row(RowType type) {
    if (type == COMPACT_ROW) {
        return new CompactAlignmentRow;
    } else {
        // default = MAP_ROW
        return new MapAlignmentRow;
    }
}

RowType AlignmentRow::type() const {
    return type_impl();
}

AlignmentRow* AlignmentRow::clone() const {
    AlignmentRow* result = new_row(type());
    result->assign(*this);
    return result;
}

AlignmentRow* AlignmentRow::slice(int start, int stop) const {
    BOOST_ASSERT(stop < length());
    BOOST_ASSERT(start < length());
    int min = std::min(start, stop);
    int max = std::max(start, stop);
    int ori = (min == start) ? 1 : -1;
    AlignmentRow* new_row = AlignmentRow::new_row(type());
    int l = max - min + 1;
    new_row->set_length(l);
    int fragment_pos = 0;
    for (int new_row_pos = 0; new_row_pos < l; new_row_pos++) {
        int old_row_pos = start + new_row_pos * ori;
        if (map_to_fragment(old_row_pos) != -1) {
            new_row->bind(fragment_pos, new_row_pos);
            fragment_pos += 1;
        }
    }
    return new_row;
}

MapAlignmentRow::MapAlignmentRow(const std::string& alignment_string,
                                 Fragment* fragment):
    AlignmentRow(fragment) {
    grow(alignment_string);
}

void MapAlignmentRow::clear() {
    fragment_to_alignment_.clear();
    alignment_to_fragment_.clear();
    set_length(0);
}

void MapAlignmentRow::bind(int fragment_pos, int align_pos) {
    fragment_to_alignment_[fragment_pos] = align_pos;
    alignment_to_fragment_[align_pos] = fragment_pos;
}

int MapAlignmentRow::map_to_alignment(int fragment_pos) const {
    if (fragment_pos >= length() || fragment_pos < 0) {
        return -1;
    }
    if (fragment() && fragment_pos >= fragment()->length()) {
        return -1;
    }
    Pos2Pos::const_iterator it2 = fragment_to_alignment_.find(fragment_pos);
    if (it2 == fragment_to_alignment_.end()) {
        return -1;
    } else {
        return it2->second;
    }
}

int MapAlignmentRow::map_to_fragment(int align_pos) const {
    if (align_pos >= length() || align_pos < 0) {
        return -1;
    }
    Pos2Pos::const_iterator it2 = alignment_to_fragment_.find(align_pos);
    if (it2 == alignment_to_fragment_.end()) {
        return -1;
    } else {
        return it2->second;
    }
}

RowType MapAlignmentRow::type_impl() const {
    return MAP_ROW;
}

CompactAlignmentRow::CompactAlignmentRow(const std::string& alignment_string,
        Fragment* fragment):
    AlignmentRow(fragment) {
    grow(alignment_string);
}

void CompactAlignmentRow::clear() {
    data_.clear();
    set_length(0);
}

void CompactAlignmentRow::bind(int /* fragment_pos */, int align_pos) {
    Chunk& c = chunk(chunk_index(align_pos));
    int internal_pos = pos_in_chunk(align_pos);
    c.set(internal_pos);
}

static struct ChunkCompare {
    typedef CompactAlignmentRow::Chunk Chunk;
    bool operator()(const Chunk& c1, int pos) const {
        return c1.pos_in_fragment > pos;
    }
} cc;

int CompactAlignmentRow::map_to_alignment(int fragment_pos) const {
    if (fragment_pos >= length() || fragment_pos < 0) {
        return -1;
    }
    if (fragment() && fragment_pos >= fragment()->length()) {
        return -1;
    }
    Data::const_reverse_iterator it = std::lower_bound(data_.rbegin(),
                                      data_.rend(), fragment_pos, cc);
    if (it == data_.rend()) {
        return -1;
    } else {
        const Chunk& c = *it;
        int internal_pos = fragment_pos - c.pos_in_fragment;
        BOOST_ASSERT(0 <= internal_pos && internal_pos < BITS_IN_CHUNK);
        return to_align_pos(&c) + c.map_to_alignment(internal_pos);
    }
}

int CompactAlignmentRow::map_to_fragment(int align_pos) const {
    if (align_pos >= length() || align_pos < 0) {
        return -1;
    }
    int index = chunk_index(align_pos);
    if (index >= data_.size()) {
        return -1;
    }
    int internal_pos = pos_in_chunk(align_pos);
    const Chunk& chunk = data_[chunk_index(align_pos)];
    int shift = chunk.map_to_fragment(internal_pos);
    return shift == -1 ? -1 : chunk.pos_in_fragment + shift;
}

RowType CompactAlignmentRow::type_impl() const {
    return COMPACT_ROW;
}

CompactAlignmentRow::Chunk::Chunk():
    pos_in_fragment(0), bitset(0)
{ }

int CompactAlignmentRow::Chunk::size() const {
    int s = 0;
    for (int i = 0; i < BITS_IN_CHUNK; i++) {
        if (get(i)) {
            s += 1;
        }
    }
    return s;
}

int CompactAlignmentRow::Chunk::map_to_alignment(int fragment_pos) const {
    BOOST_ASSERT(fragment_pos < BITS_IN_CHUNK);
    int fragment = 0;
    for (int i = 0; i < BITS_IN_CHUNK; i++) {
        if (get(i)) {
            if (fragment == fragment_pos) {
                return i;
            }
            fragment += 1;
        }
    }
    BOOST_ASSERT(fragment_pos == 0);
    return 0;
}

int CompactAlignmentRow::Chunk::map_to_fragment(int align_pos) const {
    BOOST_ASSERT(align_pos < BITS_IN_CHUNK);
    if (!get(align_pos)) {
        return -1;
    }
    int result = 0;
    for (int i = 0; i <= align_pos; i++) {
        if (get(i)) {
            result += 1;
        }
    }
    return result - 1;
}

bool CompactAlignmentRow::Chunk::get(int align_pos) const {
    return (bitset >> align_pos) & 0x01;
}

void CompactAlignmentRow::Chunk::set(int align_pos) {
    bitset |= (0x01 << align_pos);
}

int CompactAlignmentRow::chunk_index(int align_pos) {
    return align_pos / BITS_IN_CHUNK;
}

int CompactAlignmentRow::pos_in_chunk(int align_pos) {
    return align_pos % BITS_IN_CHUNK;
}

CompactAlignmentRow::Chunk& CompactAlignmentRow::chunk(int index) {
    if (index >= data_.size()) {
        Chunk c;
        if (!data_.empty()) {
            const Chunk&  back = data_.back();
            c.pos_in_fragment = back.pos_in_fragment + back.size();
        }
        data_.resize(index + 1, c);
    }
    return data_[index];
}

int CompactAlignmentRow::to_align_pos(const Chunk* chunk) const {
    return (reinterpret_cast<const char*>(chunk) -
            reinterpret_cast<const char*>(&data_[0]))
           / sizeof(Chunk) * BITS_IN_CHUNK;
}

InversedRow::InversedRow(AlignmentRow* source):
    source_(0), fragment_length_(0) {
    set_source(source);
}

InversedRow::~InversedRow() {
    delete source_;
    source_ = 0;
}

void InversedRow::clear() {
    throw "Tried to clear InversedRow";
}

void InversedRow::bind(int fragment_pos, int align_pos) {
    throw "Tried to bind InversedRow";
}

int InversedRow::map_to_alignment(int fragment_pos) const {
    if (fragment_pos >= length() || fragment_pos < 0) {
        return -1;
    }
    if (fragment_pos >= fragment_length_) {
        return -1;
    }
    BOOST_ASSERT(fragment_pos >= 0);
    BOOST_ASSERT(fragment_pos < fragment_length_);
    fragment_pos = fragment_length_ - fragment_pos - 1;
    int align_pos = source()->map_to_alignment(fragment_pos);
    if (align_pos == -1) {
        return -1;
    } else {
        BOOST_ASSERT(align_pos >= 0);
        BOOST_ASSERT(align_pos < length());
        align_pos = length() - align_pos - 1;
        return align_pos;
    }
}

int InversedRow::map_to_fragment(int align_pos) const {
    if (align_pos >= length() || align_pos < 0) {
        return -1;
    }
    BOOST_ASSERT(align_pos >= 0);
    BOOST_ASSERT(align_pos < length());
    align_pos = length() - align_pos - 1;
    int fragment_pos = source()->map_to_fragment(align_pos);
    if (fragment_pos == -1) {
        return -1;
    } else {
        BOOST_ASSERT(fragment_pos >= 0);
        BOOST_ASSERT(fragment_pos < fragment_length_);
        fragment_pos = fragment_length_ - fragment_pos - 1;
        return fragment_pos;
    }
}

AlignmentRow* InversedRow::source() const {
    return source_;
}

void InversedRow::set_source(AlignmentRow* source) {
    BOOST_ASSERT(source->fragment());
    fragment_length_ = source->fragment()->length();
    source->fragment()->detach_row();
    delete source_;
    source_ = source;
    set_length(source_->length());
}

void InversedRow::detach_source() {
    source_ = 0;
}

RowType InversedRow::type_impl() const {
    return source()->type();
}

}

