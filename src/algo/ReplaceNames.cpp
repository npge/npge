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
#include "throw_assert.hpp"

namespace npge {

ReplaceNames::ReplaceNames():
    table_(this, "table", "Table of genomes") {
    declare_bs("target", "Modified blockset");
}

typedef std::map<std::string, std::string> S2S;

static void read_table(
    const Processor* p,
    S2S& old2new,
    S2S& new2ac,
    std::istream& input) {
    for (std::string line; std::getline(input, line);) {
        using namespace boost::algorithm;
        trim(line);
        if (!line.empty()) {
            SequenceParams par(line);
            if (par.id_in_file_.empty()) {
                p->write_log("Can't parse table row: " + line);
                continue;
            }
            std::string new_name = par.genome_ + "&" +
                                   par.chromosome_ + "&" +
                                   par.circular_;
            if (par.record_type_ == "fasta" ||
                    par.record_type_ == "all") {
                old2new[par.id_in_file_] = new_name;
            }
            if (par.record_type_ == "features" ||
                    par.record_type_ == "all") {
                new2ac[new_name] = par.id_in_file_;
            }
        }
    }
}

static void replace_names(
    const Processor* p,
    const S2S& old2new,
    S2S& new2ac) {
    BlockSet& bs = *(p->block_set());
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        std::string old_name = seq->name();
        BOOST_FOREACH (const S2S::value_type& n2n, old2new) {
            const std::string& old_name1 = n2n.first;
            const std::string& new_name = n2n.second;
            if (old_name.find(old_name1) != std::string::npos) {
                seq->set_name(new_name);
                std::string d = seq->description();
                using namespace boost::algorithm;
                trim(d);
                std::string ac = new2ac[new_name];
                if (!ac.empty()) {
                    d = "ac=" + ac + " " + d;
                }
                seq->set_description(d);
                break;
            }
        }
    }
}

void ReplaceNames::run_impl() const {
    S2S old2new, new2ac;
    std::istream& input = table_.input();
    read_table(this, old2new, new2ac, input);
    replace_names(this, old2new, new2ac);
}

const char* ReplaceNames::name_impl() const {
    return "Replace names in downloaded sequences "
           "according to table";
}

}

