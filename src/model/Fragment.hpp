/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FRAGMENT_HPP_
#define BR_FRAGMENT_HPP_

#include <iosfwd>
#include <new>
#include <string>

#include "global.hpp"

namespace bloomrepeats {

/** Difference between two fragments.
Logical first and last positions are taken into account.
\see diff_to(), patch()
*/
struct FragmentDiff {
#ifndef DOXYGEN_ONLY
    int begin; /**< Difference of begin */
    int last; /**< Difference of last */
    int ori; /**< 1 if same, -1 otherwise */
#endif
};

/** Part of sequence.
  - min_pos and max_pos
  - ori
  - block
  - prev and next
*/
class Fragment {
public:
    /** Difference between two fragments */
    typedef FragmentDiff Diff;

    /** Invalid fragment */
    static const Fragment INVALID;

    /** Constructor.
    \param seq Sequence
    \param min_pos Minimal position of sequence occupied by the fragment
    \param max_pos Maximal position of sequence occupied by the fragment
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

    /** Allocate storage */
    void* operator new(size_t x);

    /** Deallocate storage */
    void operator delete(void* ptr);

    /** Get sequence */
    Sequence* seq() const {
        return seq_;
    }

    /** Get block, if any.
    \see Block::insert
    */
    Block* block() const;

    /** Get previous fragment, if any.
    \see connect()
    */
    Fragment* prev() const;

    /** Get next fragment, if any.
    \see connect()
    */
    Fragment* next() const;

    /** Get next (ori=1) or previous (ori=-1) fragment */
    Fragment* neighbor(int ori) const;

    /** Get next or previous taking fragment ori into account.
    This is an equivalent to \code neighbor(ori() * ori) \endcode
    */
    Fragment* logical_neighbor(int ori) const;

    /** Return if a fragment is previous or next fragment of this fragment */
    bool is_neighbor(const Fragment& other) const;

    /** Return another neighbor of this fragment.
    Given other fragment must be a neighbor of this fragment.
    */
    Fragment* another_neighbor(const Fragment& other) const;

    /** Get minimal position of sequence occupied by the fragment */
    size_t min_pos() const {
        return min_pos_;
    }

    /** Set minimal position of sequence occupied by the fragment */
    void set_min_pos(size_t min_pos) {
        min_pos_ = min_pos;
    }

    /** Get maximal position of sequence occupied by the fragment */
    size_t max_pos() const {
        return max_pos_;
    }

    /** Set maximal position of sequence occupied by the fragment */
    void set_max_pos(size_t max_pos) {
        max_pos_ = max_pos;
    }

    /** Get orientation (1 for forward, -1 for reverse) */
    int ori() const;

    /** Get the number of sequence positions occupied by the fragment */
    size_t length() const;

    /** Get orientation (1 for forward, -1 for reverse) */
    void set_ori(int ori);

    /** Change orientation to the opposite */
    void inverse();

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

    /** Get end of of sequence occupied by the fragment.
    STL-like end (first outside).
    */
    size_t end_pos() const;

    /** Return string representation of the fragment */
    std::string str() const;

    /** Return string representation of fragment part.
    \param from Beginning position in fragment
    \param to end position in fragment
    */
    std::string substr(int from, int to) const;

    /** Create new slice of this fragment.
    \param from Beginning position in fragment (>= 0).
    \param to end position in fragment (>= 0).
    If from > to, then the resulting fragment will be inversed.
    \warning This allocates new Fragment. Make sure it is not leaked.
    */
    Fragment* subfragment(size_t from, size_t to) const;

    /** Return fragment identifier.
    if is formed from sequence name, begin_pos() and f.last_pos(),
    separated by '_'.
    */
    std::string id() const;

    /** Return hash of this fragment */
    size_t hash() const;

    /** Change fragment size.
    \param shift Difference of fragment length.
    Beginning position is constant.
    End position is shifted.
    */
    void shift_end(int shift = 1);

    /** Max valid shift of the fragment.
    \param max_overlap Max number of positions, that are allowed to be added
       to the block after first overlap occured.
       -1 means "overlaps of any length are allowed".
       Fragments must be \ref BlockSet::connect_fragments "connected"
       for this to work correctly.

    Return max value, that can be passed to shift_end(),
    keeping the fragment valid().
    May be negative, if the fragment is already invalid.
    */
    int max_shift_end(int max_overlap = 0) const;

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

    /** Make first fragment be previous of second and second -- next of first */
    static void connect(Fragment* first, Fragment* second);

    /** Behaves as connect() if ori == 1, else vice-versa */
    static void connect(Fragment* first, Fragment* second, int ori);

    /** Swap this and other positions (prev, next) */
    void rearrange_with(Fragment* other);

    /** Rearrange this fragment before or after its neighbors */
    void find_place();

    /** Disconnect this fragment, connect it near start_from and find_place() */
    void find_place(Fragment* start_from);

    /** Return if two fragments can be joined.
    \param one Fragment.
    \param another Fragment.

    Fragments can be joined if they share the same sequence and ori
    and \ref is_neighbor "are neighbors".

    \see JoinApprover
    */
    static bool can_join(Fragment* one, Fragment* another);

    /** Merge fragments and return new larger fragment.
    \warning Fragments must be \ref can_join "joinable".
    */
    static Fragment* join(Fragment* one, Fragment* another);

    /** Disconnect this fragment from its neighbors.
    \param connect_neighbors If has prev() and next(), they would be connected
    */
    void disconnect(bool connect_neighbors = true);

    /** Return number of positions, occupied by both fragments */
    size_t common_positions(const Fragment& other);

    /** Return number of positions between two fragments */
    size_t dist_to(const Fragment& other);

    /** Return fragment of positions, occupied by both fragments.
    New fragment inherits ori from this fragment.

    If input fragments do not have common positions, empty pointer is returned.
    */
    Fragment common_fragment(const Fragment& other);

    /** Return if this fragment belongs to other fragment.
    This means, common_positions(other) == this->length().
    */
    bool is_subfragment_of(const Fragment& other);

    /** Return if this belongs to other and does not share boundaries with it.
    If this method returns true, then is_subfragment_of must also be true.
    */
    bool is_internal_subfragment_of(const Fragment& other);

    /** Return difference, which can be applied to this to get other.
    \warning Fragments MUST be of same sequence.
    \see patch()
    */
    Diff diff_to(const Fragment& other) const;

    /** Apply difference to this fragment.
    \see diff_to()
    */
    void patch(const Diff& diff);

    /** Copy seq, min_pos, max_pos and ori from other fragment */
    void apply_coords(const Fragment& other);

    /** Assignment operator (same as apply_coords) */
    Fragment& operator=(const Fragment& other);

    /** Exclude positions of other fragment from this fragment.
    If other is strongly inside this, one of "flank" fragments is produced.

    If this is inside other, \ref valid() "invalid" fragment is produced.

    This method keeps ori unchanged.

    \verbatim
    This   : ---xxx---
    Other  : -----x---
    Result : ---xx----

    This   : ---xxx---
    Other  : ----x----
    Result : ---x-----
      or     -----x---

    This   : ---xxx---
    Other  : -------x-
    Result : ---xxx---

    This   : ---xxx---
    Other  : -xxxxxxx-
    Result :  invalid
    \endverbatim
    \warning Fragments MUST be of same sequence.
    */
    void exclude(const Fragment& other);

    /** Return diff, applying of which is same to exclude() */
    Diff exclusion_diff(const Fragment& other) const;

    /** Split this fragment into two fragments.
    End of this fragment is changed so that its new length is \p new_length.
    Another part of this fragment (if any) is returned (with ori of this).

    This method \ref find_place() "finds place" for this and new fragment.

    If \p new_length >= length(), nothing is done, empty pointer is returned.
    */
    Fragment* split(size_t new_length);

    /** Return if two fragments can be aligned */
    bool aligned(const Fragment& other, PairAligner* pa = 0, int batch = 100);

    /** Output id() and description.
    Description includes "block=... prev=... next=...".
    \warning Leading '>' is not printed.
    */
    void print_header(std::ostream& o) const;

private:
    Sequence* seq_;
    Block* block_and_ori_; // pointer XOR (ori == 1 ? 0x01 : 0x00)
    Fragment* prev_;
    Fragment* next_;
    size_t min_pos_;
    size_t max_pos_;

    void set_block(Block* block);

    Block* block_raw_ptr() const;

    friend class Block;
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Fragment& fragment);

}

#endif

