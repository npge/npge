/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ANCHOR_HINDER_HPP_
#define BR_ANCHOR_HINDER_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace bloomrepeats {

/** Finder of short anchors.

For large repeats one short part is selected and returned as a anchor.

\note All hairpins are considered anchors.

\note Using >= 2 workers may (very unlikely) cause races,
    since bloom filter is not protected by a mutex.
    Such a races may cause some anchors not to be found.
\note The smallest piece of work, passed to a worker, is one sequence.
    So it is useless to set workers > sequences.
*/
class AnchorFinder : public Processor {
public:
    /** Default anchor size */
    static const size_t ANCHOR_SIZE = 20;

    /** Default constructor */
    AnchorFinder();

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

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    /** Find anchors in added sequence.
    Each found anchor is inserted into block_set().
    */
    bool run_impl() const;

    const char* name_impl() const;

private:
    size_t anchor_size_;
    int add_ori_;
    int only_ori_;
};

}

#endif

