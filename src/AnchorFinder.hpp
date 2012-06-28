/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ANCHOR_HINDER_HPP_
#define BR_ANCHOR_HINDER_HPP_

#include <vector>
#include <boost/function.hpp>

#include "global.hpp"

namespace bloomrepeats {

/** Finder of short anchors.

For large repeats one short part is selected and returned as a anchor.

\note All hairpins are considered anchors.
*/
class AnchorFinder {
public:
    /** Function, called with an anchor */
    typedef boost::function<void(BlockPtr)> AnchorHandler;

    /** Default anchor size */
    static const size_t ANCHOR_SIZE = 20;

    /** Default constructor */
    AnchorFinder();

    /** Add options to options description */
    void add_options(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options(po::variables_map& vm);

    /** Add sequence.
    If a block set was \ref set_block_set "set",
    the sequence is added to it too.
    */
    void add_sequence(SequencePtr sequence);

    /** Find anchors in added sequence.
    Each found anchor is passed to anchor_handler.

    \note If no anchor handler has been set, this method does nothing.
    */
    void run();

    /** Set block set to which blocks will be added.
    \note This method invalidates previous called set_anchor_handler().
    */
    void set_block_set(BlockSetPtr block_set);

    /** Set anchor handler function.
    Sequences are \ref BlockSet::add_sequence "added" to block set as well.
    \note This method invalidates previous called set_block_set().
    */
    void set_anchor_handler(const AnchorHandler& anchor_handler) {
        anchor_handler_ = anchor_handler;
    }

    /** Get min number of fragments in a block to accept this block */
    size_t min_fragments() const {
        return min_fragments_;
    }

    /** Set min number of fragments in a block to accept this block.
    Defaults to 2.
    */
    void set_min_fragments(size_t min_fragments) {
        min_fragments_ = min_fragments;
    }

    /** Get anchor size */
    size_t anchor_size() const {
        return anchor_size_;
    }

    /** Set anchor size.
    Defaults to ANCHOR_SIZE.
    */
    void set_anchor_size(size_t anchor_size) {
        anchor_size_ = anchor_size;
    }

    /** Return if palindromes in anchors are disabled */
    bool palindromes_elimination() const;

    /** Enable or disable palindromes in anchors.
    Palindromes (hairpins) are sequences like "ataggttaatattaacctat".
    Complementary sequence of a palindrome is equal to this palindrome.

    When palindrome elimination is enabled (by default),
    only one (ori=1 or ori=-1) fragment of palindrome sequence is added
    to a block.

    Note that palindromes are eliminated for odd anchor_size,
    regardless of this option.

    Defaults to true (palindrome
    */
    void set_palindromes_elimination(bool eliminate);

    /** Return the only ori considered, if any.
    \see set_only_ori()
    */
    int only_ori() const {
        return only_ori_;
    }

    /** Set the only ori considered.
    0 means "both".
    */
    void set_only_ori(int only_ori) {
        only_ori_ = only_ori;
    }

    /** Return number of threads used to find anchors */
    int workers() const {
        return workers_;
    }

    /** Set number of threads used to find anchors.
    Defaults to 1.
    \note Using >= 2 workers may (very unlikely) cause races,
        since bloom filter is not protected by a mutex.
        Such a races may cause some anchors not to be found.
    \note The smallest piece of work, passed to a worker, is one sequence.
        So it is useless to set workers > sequences.
    */
    void set_workers(int workers) {
        workers_ = workers;
    }

private:
    AnchorHandler anchor_handler_;
    BlockSetPtr block_set_;
    std::vector<SequencePtr> seqs_;
    size_t min_fragments_;
    size_t anchor_size_;
    int add_ori_;
    int only_ori_;
    int workers_;
};

}

#endif

