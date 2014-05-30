/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_META_PROCESSOR_HPP_
#define BR_META_PROCESSOR_HPP_

#include "Processor.hpp"

namespace npge {

/** Any processor.
Processor is specified through an option --[prefix]processor.
Options (Processor::set_options()) for internal processor are specified
through an option --[prefix]opts.

run() will fail if no processor was set.
*/
class MetaProcessor : public Processor {
public:
    /** Constructor */
    MetaProcessor(const std::string& prefix = "",
                  const std::string& processor = "",
                  const std::string& opts = "");

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    mutable Processor* p_;
};

}

#endif

