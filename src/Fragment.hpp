/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FRAGMENT_HPP_
#define BR_FRAGMENT_HPP_

#include <iosfwd>
#include <string>
#include <boost/weak_ptr.hpp>

#include "global.hpp"

namespace bloomrepeats {

/** Part of sequence.
  - min_pos and max_pos
  - ori
  - block
  - prev and next
*/
class Fragment {
public:
    /** Constructor.
    \param seq Sequence
    \param min_pos Minimal position of sequence occupied by the fragment
    \param max_pos Maximal position of sequence occupied by the fragment
    \param ori 1 for forward orientation, -1 for reverse
    */
    Fragment(SequencePtr seq = SequencePtr(),
             size_t min_pos = 0, size_t max_pos = 0, int ori = 1);

    /** Get sequence */
    SequencePtr seq() const {
        return seq_;
    }

    /** Get block, if any.
    \see Block::insert
    */
    BlockPtr block() const;

    /** Get previous fragment, if any.
    \see connect()
    */
    FragmentPtr prev() const;

    /** Get next fragment, if any.
    \see connect()
    */
    FragmentPtr next() const;

    /** Get next (ori=1) or previous (ori=-1) fragment */
    FragmentPtr neighbour(int ori) const;

    /** Get next or previous taking fragment ori into account.
    This is an equivalent to \code neighbour(ori() * ori) \endcode
    */
    FragmentPtr logical_neighbour(int ori) const;

    /** Return if a fragment is previous or next fragment of this fragment */
    bool is_neighbour(const Fragment& other) const;

    /** Return another neighbour of this fragment.
    Given other fragment must be a neighbour of this fragment.
    */
    FragmentPtr another_neighbour(const Fragment& other) const;

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
    int ori() const {
        return ori_;
    }

    /** Get the number of sequence positions occupied by the fragment */
    size_t length() const;

    /** Get orientation (1 for forward, -1 for reverse) */
    void set_ori(int ori) {
        ori_ = ori;
    }

    /** Change orientation to the opposite */
    void inverse();

    /** Get beginning position of sequence occupied by the fragment */
    size_t begin_pos() const;

    /** Get pointer to beginning of sequence occupied by the fragment */
    const char* begin() const;

    /** Get end position of sequence occupied by the fragment */
    size_t end_pos() const;

    /** Get pointer to end of sequence occupied by the fragment */
    const char* end() const;

    /** Return string representation of the fragment */
    std::string str() const;

    /** Return string representation of fragment part.
    \param from Beginning position in fragment
    \param to end position in fragment
    */
    std::string substr(int from, int to) const;

    /** Return hash of this fragment */
    size_t hash() const;

    /** Change fragment size.
    \param shift Difference of fragment length.
    Beginning position is constant.
    End position is shifted.
    */
    void shift_end(int shift = 1);

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
    static void connect(FragmentPtr first, FragmentPtr second);

    /** Return if two fragments can be merged.
    Fragments can be merged if they share the same sequence and ori
    and \ref is_neighbour "are neighbours".
    */
    static bool can_merge(FragmentPtr one, FragmentPtr another);

    /** Merge fragments and return new larger fragment.
    \warning Fragments must be \ref can_merge "mergeable".
    */
    static FragmentPtr merge(FragmentPtr one, FragmentPtr another);

    /** Disconnect this fragment from its neighbours */
    void disconnect();

    /** Return number of positions, occupied by both fragments */
    size_t common_positions(const Fragment& other);

private:
    SequencePtr seq_;
    size_t min_pos_;
    size_t max_pos_;
    int ori_;
    boost::weak_ptr<Block> block_;
    boost::weak_ptr<Fragment> prev_;
    boost::weak_ptr<Fragment> next_;

    friend class Block;
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Fragment& fragment);

}

#endif

