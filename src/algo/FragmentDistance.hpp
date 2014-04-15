/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FRAGMENT_DISTANCE_HPP_
#define BR_FRAGMENT_DISTANCE_HPP_

#include "AbstractOutput.hpp"
#include "global.hpp"

namespace bloomrepeats {

/** Calculate distance between fragments in block.
\warning Fragmnets must be aligned! Otherwise Exception is thrown.
*/
class FragmentDistance : public AbstractOutput {
public:
    /** Constructor */
    FragmentDistance();

    /** Structure describing distance between fragments */
    struct Distance {
        int penalty; /**< Penalty, "numerator" */
        int total; /**< Number of all columns, "denominator" */
        double ratio() const; /**< Distance as number in [0, 1] */
    };

    /** Distance between fragments in block.
    Long gaps are counted as one mutation.
    \warning Fragments must be aligned! Otherwise Exception is thrown.
    */
    Distance fragment_distance(const Fragment* a, const Fragment* b) const;

    /** Print table block - fr1 - fr2 - distance */
    void print_block(std::ostream& o, Block* block) const;

    void print_header(std::ostream& o) const;

protected:
    const char* name_impl() const;
};

}

#endif

