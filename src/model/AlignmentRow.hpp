/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ALIGNMENT_ROW_HPP_
#define NPGE_ALIGNMENT_ROW_HPP_

#include <map>
#include <vector>
#include <string>

#include "global.hpp"

namespace npge {

class AlignmentRow {
public:
    AlignmentRow(Fragment* fragment = 0);

    virtual ~AlignmentRow();

    virtual void clear() = 0;

    /** Grow alignment row with string representing a part of alignment.
    If debug and fragment(), correspondece with Fragment.raw_at() is tested.
    */
    virtual void grow(const std::string& alignment_string);

    virtual void bind(int fragment_pos, int align_pos) = 0;

    /** Return position in alignment, corresponding to position in fragment.
    In case of non-existing position in fragment return -1.
    */
    virtual int map_to_alignment(int fragment_pos) const = 0;

    virtual int map_to_fragment(int align_pos) const = 0;

    int length() const {
        return length_;
    }

    void set_length(int length) {
        length_ = length;
    }

    Fragment* fragment() const {
        return fragment_;
    }

    virtual int nearest_in_fragment(int align_pos) const;

    virtual void assign(const AlignmentRow& other,
                        int start = 0, int stop = -1);

    static AlignmentRow* new_row(RowType type);

    AlignmentRow* clone() const;

    AlignmentRow* slice(int min, int max) const;

    RowType type() const;

protected:
    virtual RowType type_impl() const = 0;

private:
    int length_;
    Fragment* fragment_;

    void set_fragment(Fragment* fragment) {
        fragment_ = fragment;
    }

    friend class Fragment;
};

class MapAlignmentRow : public AlignmentRow {
public:
    MapAlignmentRow(const std::string& alignment_string = "",
                    Fragment* fragment = 0);

    void clear();

    void bind(int fragment_pos, int align_pos);

    int map_to_alignment(int fragment_pos) const;

    int map_to_fragment(int align_pos) const;

protected:
    RowType type_impl() const;

private:
    typedef std::map<int, int> Pos2Pos;

    Pos2Pos fragment_to_alignment_;
    Pos2Pos alignment_to_fragment_;
};

typedef unsigned int CAR_Bitset;
const int BITS_IN_CHUNK = sizeof(CAR_Bitset) * 8;

class CompactAlignmentRow : public AlignmentRow {
public:
    CompactAlignmentRow(const std::string& alignment_string = "",
                        Fragment* fragment = 0);

    void clear();

    // TODO Currently works only forward
    void bind(int fragment_pos, int align_pos);

    int map_to_alignment(int fragment_pos) const;

    int map_to_fragment(int align_pos) const;

protected:
    RowType type_impl() const;

private:
    typedef CAR_Bitset Bitset;
    typedef unsigned int Index;
    struct Chunk {
        Index pos_in_fragment;
        Bitset bitset;

        Chunk();

        int size() const;
        int map_to_alignment(int fragment_pos) const;
        int map_to_fragment(int align_pos) const;

        bool get(int align_pos) const;
        void set(int align_pos); // TODO value = true|false
    };
    typedef std::vector<Chunk> Data;

    Data data_;

    static int chunk_index(int align_pos);
    static int pos_in_chunk(int align_pos);
    Chunk& chunk(int index);
    int to_align_pos(const Chunk* chunk) const;

    friend struct ChunkCompare;
};

/** Proxy class for inversed row.
Read-only.
*/
class InversedRow : public AlignmentRow {
public:
    /** Constructor.
    Owns source row.
    Detach source from fragment.
    Source row must have fragment.
    */
    InversedRow(AlignmentRow* source);

    /** Destructor */
    ~InversedRow();

    /** throws */
    void clear();

    /** throws */
    void bind(int fragment_pos, int align_pos);

    int map_to_alignment(int fragment_pos) const;

    int map_to_fragment(int align_pos) const;

    AlignmentRow* source() const;

    /** Set new source.
    Delete previous source.
    Detach source from fragment.
    Source row must have fragment.
    */
    void set_source(AlignmentRow* source);

    void detach_source();

protected:
    RowType type_impl() const;

private:
    AlignmentRow* source_;
    int fragment_length_;
};

}

#endif

