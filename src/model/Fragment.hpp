/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FRAGMENT_HPP_
#define NPGE_FRAGMENT_HPP_

#include <iosfwd>
#include <new>
#include <string>

#include "global.hpp"

namespace npge {

/** Part of sequence.
  - min_pos and max_pos
  - seq
  - ori
  - block
  - row
*/
class Fragment {
public:
    /** Invalid fragment */
    static const Fragment INVALID;

    /** Constructor.
    \param seq Sequence
    \param min_pos Minimum position of sequence occupied by the fragment
    \param max_pos Maximum position of sequence occupied by the fragment
    \param ori 1 for forward orientation, -1 for reverse
    */
    Fragment(Sequence* seq = 0,
             size_t min_pos = 0, size_t max_pos = 0, int ori = 1);

    /** Constructor.
    \deprecated Use Constructor(Sequence*)
    */
    Fragment(SequencePtr seq,
             size_t min_pos = 0, size_t max_pos = 0, int ori = 1);

    /** Copy constructor.
    Same as default constructor, followed by apply_coords().
    */
    Fragment(const Fragment& other);

    /** Destructor.
    Call disconnect(), \ref Block::erase "erase itself" from block.
    */
    ~Fragment();

    /** Get sequence */
    Sequence* seq() const {
        return seq_;
    }

    /** Get block, if any.
    \see Block::insert
    */
    Block* block() const;

    /** Get minimum position of sequence occupied by the fragment */
    size_t min_pos() const {
        return min_pos_;
    }

    /** Set minimum position of sequence occupied by the fragment */
    void set_min_pos(size_t min_pos) {
        min_pos_ = min_pos;
    }

    /** Get maximum position of sequence occupied by the fragment */
    size_t max_pos() const {
        return max_pos_;
    }

    /** Set maximum position of sequence occupied by the fragment */
    void set_max_pos(size_t max_pos) {
        max_pos_ = max_pos;
    }

    /** Get orientation (1 for forward, -1 for reverse) */
    int ori() const;

    /** Get the number of sequence positions occupied by the fragment */
    size_t length() const;

    /** Get the length of alignment row, fallback to length() */
    size_t alignment_length() const;

    /** Get orientation (1 for forward, -1 for reverse) */
    void set_ori(int ori, bool inverse_row = true);

    /** Change orientation to the opposite */
    void inverse(bool inverse_row = true);

    /** Get beginning position of sequence occupied by the fragment */
    size_t begin_pos() const;

    /** Set beginning position of sequence occupied by the fragment */
    void set_begin_pos(size_t begin_pos);

    /** Get last position of sequence occupied by the fragment.
    Last position (last inside).
    */
    size_t last_pos() const;

    /** Set last position of sequence occupied by the fragment */
    void set_last_pos(size_t last_pos);

    /** Set begin and last position of sequence occupied by the fragment.
    Orientation can be changed if needed.
    */
    void set_begin_last(size_t begin_pos, size_t last_pos);

    /** Get end of of sequence occupied by the fragment.
    STL-like end (first outside).
    */
    size_t end_pos() const;

    /** Return string representation of the fragment.
    \param gap Gap character
    If gap != 0 and row() != 0, then output is gapped.
    */
    std::string str(char gap = '-') const;

    /** Return string representation of fragment part.
    \param from Beginning position in fragment
    \param to end position in fragment
    */
    std::string substr(int from, int to) const;

    /** Create new slice of this fragment.
    \param from Beginning position in fragment (>= 0).
    \param to end position in fragment (< length()).
    If from > to, then the resulting fragment will be inversed.
    \warning This allocates new Fragment. Make sure it is not leaked.
    */
    Fragment* subfragment(size_t from, size_t to) const;

    /** Copy the fragment and the row */
    Fragment* clone() const;

    /** Return fragment identifier.
    if is formed from sequence name, begin_pos() and f.last_pos(),
    separated by '_'.
    */
    std::string id() const;

    /** Return hash of this fragment */
    hash_t hash() const;

    /** Return sequence name built from fragment id.
    On error, returns empty string.
    */
    static std::string seq_name_from_id(const std::string& id);

    /** Return if fragment is valid.
    Fragment is valid if and only if:
     - min_pos() <= max_pos() and
     - max_pos() < seq().size()
    */
    bool valid() const;

    /** Comparison operator.
    Fragments are considered equal, if they share the same
    sequence, ori, min_pos and max_pos.
    */
    bool operator==(const Fragment& other) const;

    /** Comparison operator */
    bool operator!=(const Fragment& other) const;

    /** Comparison operator.
     - by min_pos,
     - by max_pos (if min_pos is equal),
     - by ori (if max_pos is equal),
     - by sequence pointer.
    */
    bool operator<(const Fragment& other) const;

    /** Return if the fragment occupies the sequence position */
    bool has(size_t pos) const;

    /** Return fragment letter by index in fragment.
    Negative indexes are interpreted as is.
    */
    char raw_at(int pos) const;

    /** Return fragment letter by index in fragment.
    Negative indexes are interpreted as length - abs(pos).
    So -1 means last.
    */
    char at(int pos) const;

    /** Return fragment letter by index in fragment row.
    If row is not set, then fragment is used as is.
    Negative indexes are interpreted as length - abs(pos).
    On error return 0.
    */
    char alignment_at(int pos) const;

    /** Return number of positions, occupied by both fragments */
    size_t common_positions(const Fragment& other) const;

    /** Return number of positions between two fragments */
    size_t dist_to(const Fragment& other) const;

    /** Return fragment of positions, occupied by both fragments.
    New fragment inherits ori from this fragment.

    If input fragments do not have common positions,
    invalid fragment is returned.
    */
    Fragment common_fragment(const Fragment& other) const;

    /** Return if this fragment belongs to other fragment.
    This means, common_positions(other) == this->length().
    */
    bool is_subfragment_of(const Fragment& other) const;

    /** Return if this belongs to other and does not share boundaries with it.
    If this method returns true, then is_subfragment_of must also be true.
    */
    bool is_internal_subfragment_of(const Fragment& other) const;

    /** Copy seq, min_pos, max_pos and ori from other fragment */
    void apply_coords(const Fragment& other);

    /** Assignment operator (same as apply_coords) */
    Fragment& operator=(const Fragment& other);

    /** Return alignemnt row of this fragment */
    AlignmentRow* row() const {
        return row_;
    }

    /** Detach a row from the fragment and return it.
    Detached row is not owned by a fragment.
    */
    AlignmentRow* detach_row();

    /** Set alignemnt row of this fragment.
    Ownership is transferred.
    Previous alignemnt row is deleted if set.
    Does nothing if row == row().

    \note Alignment row is not changed by other methods of Fragment
    */
    void set_row(AlignmentRow* row);

    /** Output id() and description.
    Description includes "block=... prev=... next=...".
    If block is passed it is used, otherwise block() is used.
    \warning Leading '>' is not printed.
    */
    void print_header(std::ostream& o, const Block* block = 0) const;

    /** Print contents of fragment.
    \param o Output stream
    \param gap Gap character
    \param line Line length. 0 means infinite
    If gap != 0 and row() != 0, then output is gapped.
    */
    void print_contents(std::ostream& o, char gap = '-', int line = 60) const;

private:
    Sequence* seq_;
    Block* block_and_ori_; // pointer XOR (ori == 1 ? 0x01 : 0x00)
    size_t min_pos_;
    size_t max_pos_;
    AlignmentRow* row_;

    void set_block(Block* block);

    Block* block_raw_ptr() const;

    friend class Block;
};

/** Streaming operator.
\see Fragment::print_header()
\see Fragment::print_contents(o, '-')
*/
std::ostream& operator<<(std::ostream& o, const Fragment& fragment);

}

#endif

