/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <cctype>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "ReplaceNames.hpp"
#include "GetData.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"

namespace npge {

ReplaceNames::ReplaceNames():
    table_(this, "table", "Table of genomes") {
    declare_bs("target", "Modified blockset");
}

typedef std::map<std::string, std::string> S2S;

static void read_table(
    const Processor* p,
    S2S& name2name,
    std::istream& input) {
    for (std::string line; std::getline(input, line);) {
        using namespace boost::algorithm;
        trim(line);
        if (!line.empty()) {
            SequenceParams par(line);
            if (par.fasta_id_.empty()) {
                p->write_log("Can't parse table row: " + line);
                continue;
            }
            std::string new_name = par.genome_ + "&" +
                                   par.chromosome_ + "&" +
                                   par.circular_;
            name2name[par.fasta_id_] = new_name;
        }
    }
}

static void replace_names(
    const Processor* p,
    const S2S& name2name) {
    BlockSet& bs = *(p->block_set());
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        std::string old_name = seq->name();
        BOOST_FOREACH (const S2S::value_type& n2n, name2name) {
            const std::string& ac = n2n.first;
            const std::string& new_name = n2n.second;
            if (old_name.find(ac) != std::string::npos) {
                seq->set_name(new_name);
                std::string d = seq->description();
                using namespace boost::algorithm;
                trim(d);
                d = "ac=" + ac + " " + d;
                seq->set_description(d);
            }
        }
    }
}

void ReplaceNames::run_impl() const {
    S2S name2name;
    std::istream& input = table_.input();
    read_table(this, name2name, input);
    replace_names(this, name2name);
}

const char* ReplaceNames::name_impl() const {
    return "Replace names in downloaded sequences "
           "according to table";
}

}

