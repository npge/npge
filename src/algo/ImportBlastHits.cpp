/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "ImportBlastHits.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "convert_position.hpp"
#include "Exception.hpp"
#include "to_s.hpp"
#include "po.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

ImportBlastHits::ImportBlastHits(const BlockSetPtr& block_set, int min_length,
                                 float min_ident, float max_evalue):
    file_reader_(this, "blast-hits", "results of blast -m 8"),
    min_length_(min_length), min_ident_(min_ident), max_evalue_(max_evalue) {
    set_other(block_set);
}

void ImportBlastHits::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("blast-min-length", po::value<int>()->default_value(min_length()),
     "min length of blast hit")
    ("blast-min-ident", po::value<float>()->default_value(min_ident()),
     "min ident of blast hit")
    ("blast-max-evalue", po::value<float>()->default_value(max_evalue()),
     "max e-value of blast hit")
   ;
}

void ImportBlastHits::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("blast-min-length")) {
        int min_length = vm["blast-min-length"].as<int>();
        if (min_length < 0) {
            throw Exception("'blast-min-length' must be >= 0");
        }
        set_min_length(min_length);
    }
    if (vm.count("blast-min-ident")) {
        float min_ident = vm["blast-min-ident"].as<float>();
        if (min_ident < 0 || min_ident > 1) {
            throw Exception("'blast-min-ident' must be in [0, 1]");
        }
        set_min_ident(min_ident);
    }
    if (vm.count("blast-max-evalue")) {
        float max_evalue = vm["blast-max-evalue"].as<float>();
        if (max_evalue < 0) {
            throw Exception("'blast-max-evalue' must be >= 0");
        }
        set_max_evalue(max_evalue);
    }
}

struct BlastItem {
    std::string id;
    int start;
    int stop;

    bool operator<(const BlastItem& o) const {
        typedef boost::tuple<const std::string&, int, int> Tie;
        return Tie(id, start, stop) < Tie(o.id, o.start, o.stop);
    }
};

struct BlastHit {
    BlastHit(std::string line) {
        using namespace boost::algorithm;
        trim(line);
        std::vector<std::string> parts;
        split(parts, line, is_any_of("\t"));
        if (parts.size() < 12) {
            throw Exception("Number of fields in blast hits fasta"
                            " (" + TO_S(parts.size()) + ") "
                            "must be >= 12");
        }
        items[0].id = parts[0];
        items[1].id = parts[1];
        ident = boost::lexical_cast<float>(parts[2]) / 100;
        length = boost::lexical_cast<int>(parts[3]);
        mismatches = boost::lexical_cast<int>(parts[4]);
        gap_openings = boost::lexical_cast<int>(parts[5]);
        items[0].start = boost::lexical_cast<int>(parts[6]);
        items[0].stop = boost::lexical_cast<int>(parts[7]);
        items[1].start = boost::lexical_cast<int>(parts[8]);
        items[1].stop = boost::lexical_cast<int>(parts[9]);
        evalue = boost::lexical_cast<double>(parts[10]);
        //bit_score = boost::lexical_cast<int>(parts[11]);
    }

    BlastItem items[2];
    float ident;
    int length;
    int mismatches;
    int gap_openings;
    double evalue;
    //int bit_score;
};

typedef std::map<std::string, Block*> NameToBlock;

static void add_blast_item(const BlockSet* bs, const NameToBlock& name2block,
                           Block* new_block, const BlastItem& item) {
    if (SequencePtr seq = bs->seq_from_name(item.id)) {
        Fragment* new_fragment = new Fragment(seq);
        new_fragment->set_begin_last(item.start - 1, item.stop - 1);
        new_block->insert(new_fragment);
    } else if (Fragment* f = bs->fragment_from_id(item.id)) {
        new_block->insert(f->subfragment(item.start - 1, item.stop - 1));
        delete f;
    } else {
        NameToBlock::const_iterator it = name2block.find(item.id);
        if (it == name2block.end()) {
            throw Exception("Bad block name: " + item.id);
        }
        const Block* block = it->second;
        BOOST_ASSERT(block);
        int block_length = block->alignment_length();
        BOOST_FOREACH (Fragment* fr, *block) {
            int start = fragment_pos(fr, item.start - 1, block_length);
            BOOST_ASSERT(start != -1);
            int stop = fragment_pos(fr, item.stop - 1, block_length);
            BOOST_ASSERT(stop != -1);
            new_block->insert(fr->subfragment(start, stop));
        }
    }
}

bool ImportBlastHits::run_impl() const {
    int size_before = block_set()->size();
    NameToBlock name2block;
    BOOST_FOREACH (Block* block, *other()) {
        name2block[block->name()] = block;
    }
    BlockSet* bs = other().get();
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
        for (std::string line; std::getline(input_file, line);) {
            BlastHit hit(line);
            if (hit.items[0] < hit.items[1] &&
                    hit.length >= min_length() &&
                    hit.ident >= min_ident() &&
                    hit.evalue <= max_evalue()) {
                Block* new_block = new Block;
                add_blast_item(bs, name2block, new_block, hit.items[0]);
                add_blast_item(bs, name2block, new_block, hit.items[1]);
                block_set()->insert(new_block);
            }
        }
    }
    return block_set()->size() > size_before;
}

const char* ImportBlastHits::name_impl() const {
    return "Import blast hits";
}

}

