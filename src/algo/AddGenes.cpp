/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <istream>
#include <map>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "AddGenes.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "annotation.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"
#include "global.hpp"

namespace npge {

AddGenes::AddGenes():
    file_reader_(this, "in-genes", "input database files with genes") {
    add_opt("product", "Append product name to block name", true);
    declare_bs("target", "Blockset where genes are added");
}

typedef std::map<std::string, Sequence*> Ac2Seq;

static Sequence* find_seq(const Ac2Seq& ac2seq,
                          const std::string& ac) {
    BOOST_FOREACH (const Ac2Seq::value_type& ac_seq, ac2seq) {
        using namespace boost::algorithm;
        if (starts_with(ac_seq.first, ac)) {
            // example: CP000001.1 and CP000001
            return ac_seq.second;
        }
    }
    return 0;
}

static bool is_accession(const std::string& line) {
    using namespace boost::algorithm;
    return starts_with(line, "AC ") ||
        starts_with(line, "ACCESSION ");
}

static std::string get_accession(const std::string& line) {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, line, isspace, token_compress_on);
    ASSERT_GTE(parts.size(), 2);
    std::string ac = parts[1];
    if (ends_with(ac, ";")) {
        // cut ";"
        ac.resize(ac.size() - 1);
    }
    return ac;
}

static bool is_locus_tag(const std::string& line) {
    using namespace boost::algorithm;
    // locus_tag, refseq_locus_tag, alt_locus_tag
    size_t npos = std::string::npos;
    return line.find("/locus_tag=\"") != npos ||
           line.find("/refseq_locus_tag=\"") != npos ||
           line.find("/alt_locus_tag=\"") != npos;
}

static std::string get_locus_tag(const std::string& line) {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, line, is_any_of("\""));
    return parts[1];
}

static bool is_product(const std::string& line) {
    using namespace boost::algorithm;
    return line.length() >= 29 &&
           line.substr(21, 8) == "/product";
}

static std::string get_product(const std::string& line) {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, line, is_any_of("\""));
    return parts[1];
}

static bool is_gene(const std::string& line0) {
    using namespace boost::algorithm;
    if (line0.size() < 6) {
        return false;
    }
    // https://github.com/npge/npge/issues/19
    std::string prefix2 = line0.substr(0, 2);
    if (prefix2 != "  " && prefix2 != "FT") {
        return false;
    }
    std::string line = line0.substr(5);
    if (starts_with(line, "CDS")) {
        return true;
    }
    if (starts_with(line, "tRNA")) {
        return true;
    }
    if (starts_with(line, "rRNA")) {
        return true;
    }
    if (starts_with(line, "misc_RNA")) {
        return true;
    }
    if (starts_with(line, "ncRNA")) {
        return true;
    }
    if (starts_with(line, "tmRNA")) {
        return true;
    }
    return false;
}

static void get_gene(
    const std::string& line0,
    std::string& feature_type,
    std::string& coords) {
    using namespace boost::algorithm;
    std::string line = line0.substr(2); // remove FT
    trim(line);
    Strings parts;
    split(parts, line, isspace, token_compress_on);
    ASSERT_GTE(parts.size(), 2);
    feature_type = parts[0];
    coords = parts[1];
}

void AddGenes::run_impl() const {
    BlockSet& bs = *block_set();
    Ac2Seq ac2seq;
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        ac2seq[seq->ac()] = seq.get();
    }
    bool use_product = opt_value("product").as<bool>();
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
        Sequence* seq = 0;
        Block* b = 0;
        Block* locus_tag_block = 0;
        std::string feature_type;
        for (std::string line; std::getline(input_file, line);) {
            using namespace boost::algorithm;
            if (is_id(line)) {
                seq = 0;
            } else if (is_accession(line) && !seq) {
                std::string ac = get_accession(line);
                seq = find_seq(ac2seq, ac);
                ASSERT_MSG(seq, ("No sequence with ac=" +
                                 ac).c_str());
                b = 0;
                locus_tag_block = 0;
            } else if (b && is_locus_tag(line)) {
                std::string locus_tag = get_locus_tag(line);
                if (locus_tag_block) {
                    // append to name
                    std::string name = locus_tag_block->name();
                    if (!locus_tag.empty() &&
                            name.find(locus_tag) != std::string::npos) {
                        name += "_" + locus_tag;
                        locus_tag_block->set_name(name);
                    }
                } else {
                    b->set_name(feature_type + " " + locus_tag);
                    locus_tag_block = b;
                }
            } else if (use_product && locus_tag_block &&
                       is_product(line)) {
                std::string product = get_product(line);
                std::string locus_tag = locus_tag_block->name();
                std::string genome = seq->genome();
                locus_tag += " " + product + " (" + genome + ")";
                locus_tag_block->set_name(locus_tag);
                locus_tag_block = 0;
            } else if (is_gene(line)) {
                ASSERT_TRUE(seq);
                std::string coords;
                // feature_type declared above
                get_gene(line, feature_type, coords);
                int ori = 1;
                if (starts_with(coords, "complement(")) {
                    ori = -1;
                    int slice_begin = 11;
                    int slice_end = coords.size() - 1;
                    int slice_length = slice_end - slice_begin;
                    coords = coords.substr(slice_begin, slice_length);
                }
                if (starts_with(coords, "join(")) {
                    int slice_begin = 5;
                    int slice_end = coords.size() - 1;
                    int slice_length = slice_end - slice_begin;
                    coords = coords.substr(slice_begin, slice_length);
                }
                ASSERT_GT(coords.size(), 4);
                Strings boundaries;
                split(boundaries, coords, is_any_of(".,"),
                      token_compress_on);
                if (boundaries.size() == 1) {
                    // CDS is single number
                    // interpret it as both start and stop
                    boundaries.push_back(boundaries[0]);
                }
                ASSERT_GTE(boundaries.size(), 2);
                ASSERT_EQ(boundaries.size() % 2, 0);
                b = new Block;
                bs.insert(b);
                for (int i = 0; i < boundaries.size() / 2; i++) {
                    std::string& min_pos_str = boundaries[i * 2];
                    std::string& max_pos_str = boundaries[i * 2 + 1];
                    // <1375315..1375356
                    if (!isdigit(min_pos_str[0])) {
                        min_pos_str = min_pos_str.substr(1);
                    }
                    if (!isdigit(max_pos_str[0])) {
                        max_pos_str = max_pos_str.substr(1);
                    }
                    int min_pos = boost::lexical_cast<int>(min_pos_str) - 1;
                    int max_pos = boost::lexical_cast<int>(max_pos_str) - 1;
                    Fragment* f = new Fragment(seq, min_pos, max_pos, ori);
                    b->insert(f);
                }
            }
        }
    }
    int index = 1;
    BOOST_FOREACH (Block* block, bs) {
        std::string name = block->name();
        std::string prefix = TO_S(index);
        block->set_name(prefix + " " + name);
        index += 1;
    }
}

const char* AddGenes::name_impl() const {
    return "Add genes from EBI or GenBank genes description";
}

}

