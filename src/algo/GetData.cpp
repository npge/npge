/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "GetData.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "download_file.hpp"
#include "name_to_stream.hpp"
#include "read_file.hpp"

namespace npge {

SequenceParams::SequenceParams(const std::string& line) {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, line, isspace, token_compress_on);
    if (!starts_with(line, "#") && parts.size() >= 4) {
        std::string usa = parts[0];
        Strings usa_parts;
        split(usa_parts, usa, is_any_of(":"));
        if (usa_parts.size() == 3) {
            record_type_ = usa_parts[0];
            database_ = usa_parts[1];
            id_ = usa_parts[2];
            id_in_file_ = id_;
            if (ends_with(id_, "]")) {
                std::string i1 = id_;
                // cut ']'
                i1.resize(i1.size() - 1);
                Strings id_parts;
                split(id_parts, i1, is_any_of("["));
                if (id_parts.size() == 2) {
                    fname_ = id_parts[0];
                    id_in_file_ = id_parts[1];
                }
            }
            //
            genome_ = parts[1];
            chromosome_ = parts[2];
            circular_ = parts[3];
        }
    }
}

static bool check_type(Processor* p, std::string& m) {
    std::string t = p->opt_value("type").as<std::string>();
    if (t != "fasta" && t != "features") {
        m = "'type' must be 'fasta' or 'features'";
        return false;
    }
    return true;
}

GetData::GetData():
    table_(this, "table", "Table of genomes"),
    out_(this, "data", "Output file", true) {
    add_opt("type",
            "Type of content downloaded (fasta|features)",
            std::string("fasta"));
    add_opt_check(boost::bind(check_type, this, _1));
}

const char* DBFETCH_URL = "http://www.ebi.ac.uk/Tools/"
                          "dbfetch/dbfetch?db={db}&id={id}"
                          "&format={format}&style=raw";

void GetData::run_impl() const {
    std::istream& input = table_.input();
    // make sure output file is opened (and created)
    out_.output();
    for (std::string line; std::getline(input, line);) {
        using namespace boost::algorithm;
        trim(line);
        if (!line.empty()) {
            process_line(line);
        }
    }
}

static void read_fasta_from_file(
    std::ostream& out,
    const SequenceParams& par) {
    BlockSet bs;
    typedef boost::shared_ptr<std::istream> IPtr;
    IPtr ifile = name_to_istream(par.fname_);
    *ifile >> bs;
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        if (seq->name() == par.id_in_file_) {
            out << *seq;
            break;
        }
    }
}

static void read_features_from_file(
    std::ostream& out,
    const SequenceParams& par) {
    using namespace boost::algorithm;
    typedef boost::shared_ptr<std::istream> IPtr;
    IPtr ifile = name_to_istream(par.fname_);
    bool inside = false;
    for (std::string line; std::getline(*ifile, line);) {
        trim(line);
        if (starts_with(line, "ID ")) {
            if (!inside) {
                std::string line1 = line.substr(3);
                trim(line1);
                Strings line_parts;
                split(line_parts, line1, is_any_of(";"));
                ASSERT_GTE(line_parts.size(), 1);
                if (line_parts[0] == par.id_in_file_) {
                    inside = true;
                }
            } else {
                break;
            }
        }
        if (inside) {
            out << line << "\n";
        }
    }
}

static void read_from_file(std::ostream& out,
                           const SequenceParams& par) {
    if (!par.fname_.empty()) {
        if (par.record_type_ == "fasta") {
            read_fasta_from_file(out, par);
        } else if (par.record_type_ == "features") {
            read_features_from_file(out, par);
        }
    } else {
        out << read_file(par.id_);
    }
}

void GetData::process_line(const std::string& line) const {
    using namespace boost::algorithm;
    std::string type = opt_value("type").as<std::string>();
    SequenceParams par(line);
    if (par.id_.empty()) {
        if (!starts_with(line, "#")) {
            write_log("Can't parse table row: " + line);
        }
        return;
    }
    if (par.record_type_ != "fasta" &&
            par.record_type_ != "features" &&
            par.record_type_ != "all") {
        write_log("Unknown record type: " + par.record_type_);
        return;
    }
    if (par.record_type_ != type && par.record_type_ != "all") {
        return;
    }
    std::string format = "default";
    if (type == "fasta") {
        format = "fasta";
    }
    std::string db = par.database_;
    if (db == "refseq") {
        db = "refseqn";
    }
    if (db == "ena") {
        db = "embl";
    }
    if (db == "file") {
        read_from_file(out_.output(), par);
        return;
    }
    std::string url(DBFETCH_URL);
    replace_first(url, "{db}", db);
    replace_first(url, "{id}", par.id_);
    replace_first(url, "{format}", format);
    write_log("Downloading " + url);
    set_sstream(":downloaded");
    bool ok = download_file(url, ":downloaded");
    out_.output() << read_file(":downloaded");
    remove_stream(":downloaded");
    if (ok) {
        write_log(".. downloaded");
    } else {
        write_log(url + " - problems");
    }
}

const char* GetData::name_impl() const {
    return "Download genomes from Web";
}

}

