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

#include "GetData.hpp"
#include "curl_download.hpp"

namespace npge {

SequenceParams::SequenceParams(const std::string& line) {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, line, isspace, token_compress_on);
    if (parts.size() >= 4) {
        fasta_id_ = parts[0];
        genome_ = parts[1];
        chromosome_ = parts[2];
        circular_ = parts[3];
    }
}

static bool check_type(Processor* p, std::string& m) {
    std::string t = p->opt_value("type").as<std::string>();
    if (t != "fasta" && t != "genes") {
        m = "'type' must be 'fasta' or 'genes'";
        return false;
    }
    return true;
}

GetData::GetData():
    table_(this, "table", "Table of genomes"),
    out_(this, "data", "Output file", true) {
    add_opt("type", "Type of content downloaded (fasta|genes)",
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
    std::string format = "default";
    if (type == "fasta") {
        format = "fasta";
    }
    SequenceParams par(line);
    if (par.fasta_id_.empty()) {
        write_log("Can't parse table row: " + line);
        return;
    }
    std::string db = "embl";
    if (par.fasta_id_.size() > 3 && par.fasta_id_[2] == '_') {
        db = "refseqn";
    }
    std::string url(DBFETCH_URL);
    replace_first(url, "{db}", db);
    replace_first(url, "{id}", par.fasta_id_);
    replace_first(url, "{format}", format);
    write_log("Downloading " + url);
    bool ok = download_file(url, out_.output_file());
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

