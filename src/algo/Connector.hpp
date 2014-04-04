/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONNECTOR_HPP_
#define BR_CONNECTOR_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Connect all the fragments (prev-next).
Returns "false", because this operation
is considered to be constant.
*/
class Connector : public Processor {
public:
    /** Constructor */
    Connector();

protected:
    /** Apply the action */
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

