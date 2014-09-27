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

/** Parameters of line from table */
struct SequenceParams {
    std::string record_type_;
    std::string database_;
    std::string id_;
    std::string fname_;
    std::string id_in_file_;
    std::string genome_;
    std::string chromosome_;
    std::string circular_;

    /** Constructor */
    SequenceParams(const std::string& line);
};

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

