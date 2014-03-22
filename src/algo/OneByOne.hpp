/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ONE_BY_ONE_HPP_
#define BR_ONE_BY_ONE_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Split "target" blocks according other blocks, processed one by one.

Blocksets: "target", "other".

New blocks are aligned and checked using Filter.

Target blockset is modified.
*/
class OneByOne : public Processor {
public:
    /** Constructor */
    OneByOne();

    /** Destructor */
    ~OneByOne();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    class Impl;
    Impl* impl_;
};

}

#endif

