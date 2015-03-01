/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REPLACE_NAMES_HPP_
#define NPGE_REPLACE_NAMES_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"

namespace npge {

/** Replace names in downloaded sequences according to table.

Example.

Table:
AE017224 BRUAB chr2 c

In:
>ENA|AE017224|AE017224.1 Brucella abortus ...

Out:
>BRUAB&chr1&c ac=AE017223 Brucella abortus ...

Replace in "target" blockset.
*/
class ReplaceNames : public Processor {
public:
    /** Constructor */
    ReplaceNames();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileReader table_;
};

}

#endif

