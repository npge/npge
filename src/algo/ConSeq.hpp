/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONSEQ_HPP_
#define BR_CONSEQ_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Add consensus sequences, produced from blocks of source block set.
Depends on UniqueNames and ExternalAligner.
*/
class ConSeq : public Processor {
public:
    /** Constructor */
    ConSeq(const BlockSetPtr& source = BlockSetPtr());

protected:
    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

