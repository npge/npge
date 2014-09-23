/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <boost/bind.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "GetData.hpp"
#include "download_file.hpp"
#include "name_to_stream.hpp"
#include "read_file.hpp"

namespace npge {

SequenceParams::SequenceParams(const std::string& line) {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, line, isspace, token_compress_on);
    if (parts.size() >= 4) {
        std::string usa = parts[0];
        Strings usa_parts;
        split(usa_parts, usa, is_any_of(":"));
        if (usa_parts.size() == 3) {
            record_type_ = usa_parts[0];
            database_ = usa_parts[1];
            id_ = usa_parts[2];
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
    for (std::string line; std::getline(input, line);) {
        using namespace boost::algorithm;
        trim(line);
        if (!line.empty()) {
            process_line(line);
        }
    }
}

void GetData::process_line(const std::string& line) const {
    using namespace boost::algorithm;
    std::string type = opt_value("type").as<std::string>();
    SequenceParams par(line);
    if (par.id_.empty()) {
        write_log("Can't parse table row: " + line);
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
        write_log(url + " downloaded");
    } else {
        write_log(url + " - problems");
    }
}

const char* GetData::name_impl() const {
    return "Download genomes from Web";
}

}

