/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Processor.hpp"
#include "FileWriter.hpp"

namespace npge {

/** Print HTML help about all processors */
class AllProcessors : public Processor {
public:
    /** Constructor */
    AllProcessors();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileWriter out_;
};

/** Print HTML help about all global configuration options */
class AllOptions : public Processor {
public:
    /** Constructor */
    AllOptions();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileWriter out_;
};

}

