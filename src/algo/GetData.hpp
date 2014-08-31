/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GET_DATA_HPP_
#define NPGE_GET_DATA_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"
#include "FileWriter.hpp"

namespace npge {

/** Download genomes from Web */
class GetData : public Processor {
public:
    /** Constructor */
    GetData();

protected:
    void run_impl() const;

    const char* name_impl() const;

    /** Process one line of table file */
    void process_line(const std::string& line) const;

private:
    FileReader table_;
    FileWriter out_;
};

}

#endif

