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

Default anchor size is ANCHOR_SIZE.

<h5>Enable or disable palindromes in anchors</h5>
Option --no-palindromes.

Palindromes (hairpins) are sequences like "ATAGGTTAATATTAACCTAT".
Complementary sequence of a palindrome is equal to this palindrome.

When palindrome elimination is enabled (by default),
only one (ori=1 or ori=-1) fragment of palindrome sequence is added
to a block.

Note that palindromes are eliminated for odd anchor_size,
regardless of this option.

Defaults to true (palindromes are eliminated).
*/
class AnchorFinder : public Processor {
public:
    /** Default anchor size */
    static const size_t ANCHOR_SIZE = 20;

    /** Default constructor */
    AnchorFinder();

protected:
    /** Find anchors in added sequence.
    Each found anchor is inserted into block_set().
    */
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

