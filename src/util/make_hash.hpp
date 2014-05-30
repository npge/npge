/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MAKE_HASH_HPP_
#define NPGE_MAKE_HASH_HPP_

#include "complement.hpp"

namespace npge {

/** Make hash value from fragment of sequence.
\param start Beginning of the fragment
\param length Length of the fragment
\param ori Orientation of the fragment (1 or -1)
*/
size_t make_hash(const char* start, size_t length, int ori = 1);

/** Make hash value from previous hash value (optimization).
\param old_hash Previous hash value
\param length Length of the fragment
\param forward If moved in according with the fragment's direction
\param remove_char Nucleotide, removed from the hash value.
\param add_char New nucleotide, added to the hash value.
Nucleotides add_char and remove_char should be pre-complement'ed, if needed.
*/
size_t reuse_hash(size_t old_hash, size_t length,
                  char remove_char, char add_char, bool forward = true);

}

#endif

