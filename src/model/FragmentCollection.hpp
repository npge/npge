/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FRAGMENT_COLLECTION_HPP_
#define NPGE_FRAGMENT_COLLECTION_HPP_

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>

#include "global.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace npge {

class FragmentCompare {
public:
    bool operator()(const Fragment& a, const Fragment& b) const {
        return a < b;
    }

    bool operator()(const Fragment* a, const Fragment* b) const {
        return *a < *b;
    }
};

static FragmentCompare fc_;

template<typename F>
struct AssignFragment {
    void operator()(F& dest, Fragment* source) const;

    Fragment* operator()(const F& source) const;
};

template<>
struct AssignFragment<Fragment> {
    typedef Fragment F;

    void operator()(F& dest, Fragment* source) const {
        dest = *source;
    }

    Fragment* operator()(const F& source) const {
        return const_cast<Fragment*>(&source);
    }
};

template<>
struct AssignFragment<Fragment*> {
    typedef Fragment* F;

    void operator()(F& dest, Fragment* source) const {
        dest = source;
    }

    Fragment* operator()(const F& source) const {
        return const_cast<Fragment*>(source);
    }
};

template<typename F, typename C>
struct InsertFragment {
    void operator()(C& col, const F& fragment) const;
};

template<typename F>
struct InsertFragment<F, std::vector<F> > {
    typedef std::vector<F> C;

    void operator()(C& col, const F& f) const {
        col.push_back(f);
    }
};

template<typename F>
struct InsertFragment<F, std::set<F, FragmentCompare> > {
    typedef std::set<F, FragmentCompare> C;

    void operator()(C& col, const F& f) const {
        col.insert(f);
    }
};

template<typename F, typename C>
struct RemoveFragment {
    void operator()(C& col, const F& fragment) const;
};

template<typename F>
struct RemoveFragment<F, std::vector<F> > {
    typedef std::vector<F> C;

    void operator()(C& c, const F& f) const {
        c.erase(std::remove(c.begin(), c.end(), f), c.end());
    }
};

template<typename F>
struct RemoveFragment<F, std::set<F, FragmentCompare> > {
    typedef std::set<F, FragmentCompare> C;

    void operator()(C& col, const F& f) const {
        col.erase(f);
    }
};

template<typename C>
struct SortFragments {
    void operator()(C& col) const;
};

template<typename F>
struct SortFragments<std::vector<F> > {
    typedef std::vector<F> C;

    void operator()(C& col) const {
        std::sort(col.begin(), col.end(), fc_);
    }
};

template<typename F>
struct SortFragments<std::set<F, FragmentCompare> > {
    typedef std::set<F, FragmentCompare> C;

    void operator()(C& col) const {
        // set is already sorted
    }
};

template<typename F, typename C>
struct LowerBound {
    typename C::const_iterator
    operator()(const C& col, const F& fragment) const;
};

template<typename F>
struct LowerBound<F, std::vector<F> > {
    typedef std::vector<F> C;

    typename C::const_iterator
    operator()(const C& col, const F& fragment) const {
        return std::lower_bound(col.begin(), col.end(), fragment, fc_);
    }
};

template<typename F>
struct LowerBound<F, std::set<F, FragmentCompare> > {
    typedef std::set<F, FragmentCompare> C;

    typename C::const_iterator
    operator()(const C& col, const F& fragment) const {
        return col.lower_bound(fragment);
    }
};

/** Collection of fragments.
Template class.
First template parameter is type used to store fragments:
 - Fragment
 - Fragment*

Second template parameter is type used to store several fragments:
 - std::vector<F>
 - std::set<F, FragmentCompare>

Add fragments using add_fragment(), add_block() and add_bs().
Then call prepare().
After this you can use has_overlap() and find_overlaps() methods.

If std::set is used to store fragments, then call of prepare() is not needed.
You can add fragments and check overlaps in any order in this case.
*/
template<typename F, typename C>
class FragmentCollection {
public:
    /** Constructor */
    FragmentCollection():
        cycles_allowed_(true) {
    }

    /** Return if cycles in neighborhood are allowed */
    bool cycles_allowed() const {
        return cycles_allowed_;
    }

    /** Set if cycles in neighborhood are allowed */
    void set_cycles_allowed(bool cycles_allowed) {
        cycles_allowed_ = cycles_allowed;
    }

    /** Add a fragment to the collection */
    void add_fragment(Fragment* fragment) {
        ASSERT_TRUE(fragment);
        F f;
        assigner_(f, fragment);
        Sequence* seq = fragment->seq();
        ASSERT_TRUE(seq);
        C& col = data_[seq];
        inserter_(col, f);
    }

    /** Remove a fragment to the collection.
    Does nothing if the fragment is not in collection.
    O(N) if vector is used as storage.
    prepare() is not needed after removing.
    */
    void remove_fragment(Fragment* fragment) {
        ASSERT_TRUE(fragment);
        F f;
        assigner_(f, fragment);
        Sequence* seq = fragment->seq();
        ASSERT_TRUE(seq);
        C& col = data_[seq];
        remover_(col, f);
    }

    /** Add fragments of block to the collection */
    void add_block(Block* b) {
        ASSERT_TRUE(b);
        BOOST_FOREACH (Fragment* f, *b) {
            add_fragment(f);
        }
    }

    /** Remove fragments of block from the collection */
    void remove_block(Block* b) {
        ASSERT_TRUE(b);
        BOOST_FOREACH (Fragment* f, *b) {
            remove_fragment(f);
        }
    }

    /** Add fragments of blockset to the collection */
    void add_bs(const BlockSet& bs) {
        BOOST_FOREACH (Block* b, bs) {
            add_block(b);
        }
    }

    /** Remove fragments of blockset from the collection */
    void remove_bs(const BlockSet& bs) {
        BOOST_FOREACH (Block* b, bs) {
            remove_block(b);
        }
    }

    /** Prepare collection for overlaps search.
    This step is not needed if std::set is used to store fragments.
    */
    void prepare() {
        BOOST_FOREACH (typename Seq2Fragments::value_type& seq_and_col, data_) {
            C& col = seq_and_col.second;
            sorter_(col);
        }
    }

    /** Clear collection */
    void clear() {
        data_.clear();
    }

    /** Return if the fragment overlaps any fragment from the collection */
    bool has_overlap(Fragment* fragment) const {
        F f;
        assigner_(f, fragment);
        Sequence* seq = fragment->seq();
        typename Seq2Fragments::const_iterator it = data_.find(seq);
        if (it == data_.end()) {
            return false;
        }
        const C& fragments = it->second;
        if (fragments.empty()) {
            return false;
        }
        typename C::const_iterator i2 = lower_bound_(fragments, f);
        if (i2 != fragments.end() &&
                assigner_(*i2)->common_positions(*fragment)) {
            return true;
        } else if (i2 != fragments.begin()) {
            i2--;
            if (assigner_(*i2)->common_positions(*fragment)) {
                return true;
            }
        }
        return false;
    }

    /** Return if the block overlaps any fragment from the collection */
    bool block_has_overlap(Block* block) const {
        BOOST_FOREACH (Fragment* f, *block) {
            if (has_overlap(f)) {
                return true;
            }
        }
        return false;
    }

    /** Return if the blockset overlaps any fragment from the collection */
    bool bs_has_overlap(const BlockSet& bs) const {
        BOOST_FOREACH (Block* block, bs) {
            if (block_has_overlap(block)) {
                return true;
            }
        }
        return false;
    }

    /** Find fragments overlapping with the fragment */
    void find_overlap_fragments(Fragments& overlap_fragments,
                                Fragment* fragment) const {
        F f;
        assigner_(f, fragment);
        Sequence* seq = fragment->seq();
        typename Seq2Fragments::const_iterator it = data_.find(seq);
        if (it == data_.end()) {
            return;
        }
        const C& fragments = it->second;
        ASSERT_FALSE(fragments.empty());
        typename C::const_iterator i2 = lower_bound_(fragments, f);
        typename C::const_iterator i2r = i2, i2l = i2;
        if (i2 != fragments.end() &&
                assigner_(*i2)->common_positions(*fragment)) {
            overlap_fragments.push_back(assigner_(*i2));
        }
        while (i2l != fragments.begin()) {
            i2l--;
            if (assigner_(*i2l)->common_positions(*fragment)) {
                overlap_fragments.push_back(assigner_(*i2l));
            }
        }
        while (i2r != fragments.end()) {
            i2r++;
            if (i2r != fragments.end() &&
                    assigner_(*i2r)->common_positions(*fragment)) {
                overlap_fragments.push_back(assigner_(*i2r));
            }
        }
    }

    /** Find overlaps between the fragment and fragments from the collection.
    Ori of fragment from collection is used.
    */
    void find_overlaps(std::vector<Fragment>& overlaps,
                       Fragment* fragment) const {
        Fragments overlap_fragments;
        find_overlap_fragments(overlap_fragments, fragment);
        BOOST_FOREACH (Fragment* f, overlap_fragments) {
            ASSERT_TRUE(f->common_positions(*fragment));
            overlaps.push_back(f->common_fragment(*fragment));
            ASSERT_EQ(overlaps.back().ori(), f->ori());
        }
    }

    typedef typename C::const_iterator CIt2;
    typedef std::pair<const C*, CIt2> FrIt;

    /** Find fragment */
    FrIt find_fragment(Fragment* fragment) const {
        typedef typename Seq2Fragments::const_iterator CIt;
        Sequence* seq = fragment->seq();
        CIt it = data_.find(seq);
        if (it == data_.end()) {
            return FrIt();
        }
        const C& fragments = it->second;
        if (fragments.empty()) {
            return FrIt();
        }
        F f;
        assigner_(f, fragment);
        CIt2 i2 = lower_bound_(fragments, f);
        if (i2 == fragments.end()) {
            return FrIt();
        }
        Fragment* found = *i2;
        if (*assigner_(*i2) != *fragment) {
            return FrIt();
        }
        return FrIt(&fragments, i2);
    }

    /** Return next fragment in collection.
    If a sequence is circular and cycles_allowed(),
    then any fragment has previous
    and next fragment, otherwise first fragment has no
    previous and last fragment has no next.
    If searched fragment is not part of the collection,
    then return 0.
    */
    Fragment* next(Fragment* fragment) const {
        return next_from_it(find_fragment(fragment));
    }

    /** Return next fragment in collection */
    Fragment* next_from_it(const FrIt& frit) const {
        if (frit.first == 0) {
            return 0;
        }
        const C& fragments = *frit.first;
        CIt2 i2 = frit.second;
        ASSERT_TRUE(i2 != fragments.end());
        Sequence* seq = assigner_(*i2)->seq();
        i2++;
        if (i2 != fragments.end()) {
            return assigner_(*i2);
        } else if (cycles_allowed() && seq->circular()) {
            return assigner_(*fragments.begin());
        } else {
            return 0;
        }
    }

    /** Return prev fragment in collection */
    Fragment* prev(Fragment* fragment) const {
        return prev_from_it(find_fragment(fragment));
    }

    /** Return prev fragment in collection */
    Fragment* prev_from_it(const FrIt& frit) const {
        if (frit.first == 0) {
            return 0;
        }
        const C& fragments = *frit.first;
        CIt2 i2 = frit.second;
        ASSERT_TRUE(i2 != fragments.end());
        Sequence* seq = assigner_(*i2)->seq();
        if (frit.second != fragments.begin()) {
            i2--;
            return assigner_(*i2);
        } else if (cycles_allowed() && seq->circular()) {
            CIt2 rbegin = fragments.end();
            rbegin--;
            return assigner_(*rbegin);
        } else {
            return 0;
        }
    }

    /** Get next (ori=1) or previous (ori=-1) fragment */
    Fragment* neighbor(Fragment* fragment, int ori) const {
        if (ori == 1) {
            return next(fragment);
        } else {
            return prev(fragment);
        }
    }

    /** Get next or prev taking fragment ori into account */
    Fragment* logical_neighbor(Fragment* f, int ori) const {
        return neighbor(f, f->ori() * ori);
    }

    /** Return if fragments are neighbors and ori.
    Returns 1, if f1's next is f2.
    Returns -1, if f2's next is f1.
    Returns 0, if they are not neighbors.
    */
    int are_neighbors(Fragment* f1, Fragment* f2) const {
        FrIt it = find_fragment(f1);
        if (next_from_it(it) == f2) {
            return 1;
        } else if (prev_from_it(it) == f2) {
            return -1;
        } else {
            return 0;
        }
    }

    /** Return another neighbor of f1 */
    Fragment* another_neighbor(Fragment* f1,
                               Fragment* f2) const {
        FrIt it = find_fragment(f1);
        if (next_from_it(it) == f2) {
            return prev_from_it(it);
        } else if (prev_from_it(it) == f2) {
            return next_from_it(it);
        } else {
            return 0;
        }
    }

    /** Return list of sequences */
    std::vector<Sequence*> seqs() const {
        std::vector<Sequence*> result;
        typedef typename Seq2Fragments::value_type V;
        BOOST_FOREACH (const V& v, data_) {
            result.push_back(v.first);
        }
        return result;
    }

    /** Return if it contains this sequence */
    bool has_seq(Sequence* seq) const {
        return data_.find(seq) != data_.end();
    }

    /** Return container with fragments of the sequence */
    const C& fragments_of(Sequence* seq) const {
        typedef typename Seq2Fragments::const_iterator CIt;
        CIt it = data_.find(seq);
        ASSERT_TRUE(it != data_.end());
        return it->second;
    }

private:
    typedef std::map<Sequence*, C> Seq2Fragments;
    Seq2Fragments data_;
    AssignFragment<F> assigner_;
    InsertFragment<F, C> inserter_;
    RemoveFragment<F, C> remover_;
    SortFragments<C> sorter_;
    LowerBound<F, C> lower_bound_;
    bool cycles_allowed_;
};

typedef std::set<Fragment*, FragmentCompare> FSet;
typedef FragmentCollection<Fragment*, FSet> SetFc;
typedef std::vector<Fragment*> FVec;
typedef FragmentCollection<Fragment*, FVec> VectorFc;

typedef std::set<Fragment, FragmentCompare> DirectFSet;
typedef FragmentCollection<Fragment, DirectFSet> DirectSetFc;
typedef std::vector<Fragment> DirectFVec;
typedef FragmentCollection<Fragment, DirectFVec> DirectVectorFc;

}

#endif

