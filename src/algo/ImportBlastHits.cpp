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
#include <boost/algorithm/string/classification.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "ImportBlastHits.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "convert_position.hpp"
#include "Exception.hpp"
#include "to_s.hpp"
#include "po.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace bloomrepeats {

ImportBlastHits::ImportBlastHits(const BlockSetPtr& block_set,
                                 int min_length,
                                 double min_ident, double max_evalue):
    file_reader_(this, "blast-hits", "results of blast -m 8") {
    add_opt("blast-min-length", "min length of blast hit", min_length);
    add_opt_rule("blast-min-length >= 0");
    set_other(block_set);
    declare_bs("target", "Where hits are added");
    declare_bs("other", "Consensuses (blocks or sequences) "
               "on which blast was run");
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
        Strings parts;
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
typedef std::map<std::string, Sequence*> NameToSeq;

static void add_blast_item(const BlockSet* bs,
                           const NameToSeq& name2seq,
                           const NameToBlock& name2block,
                           Block* new_block, const BlastItem& item) {
    Sequence* seq = 0;
    NameToSeq::const_iterator it = name2seq.find(item.id);
    if (it != name2seq.end()) {
        seq = it->second;
    }
    Fragment* f = 0;
    std::string f_seq_name = Fragment::seq_name_from_id(item.id);
    if (!f_seq_name.empty()) {
        NameToSeq::const_iterator it2 = name2seq.find(f_seq_name);
        if (it2 != name2seq.end()) {
            Sequence* s = it2->second;
            f = s->fragment_from_id(item.id);
        }
    }
    if (seq) {
        Fragment* new_fragment = new Fragment(seq);
        new_fragment->set_begin_last(item.start - 1, item.stop - 1);
        new_block->insert(new_fragment);
    } else if (f) {
        new_block->insert(f->subfragment(item.start - 1,
                                         item.stop - 1));
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

void ImportBlastHits::run_impl() const {
    NameToSeq name2seq;
    BOOST_FOREACH (SequencePtr seq, other()->seqs()) {
        name2seq[seq->name()] = seq.get();
    }
    NameToBlock name2block;
    BOOST_FOREACH (Block* block, *other()) {
        name2block[block->name()] = block;
    }
    BlockSet* bs = other().get();
    int min_length = opt_value("blast-min-length").as<int>();
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
        for (std::string line; std::getline(input_file, line);) {
            BlastHit hit(line);
            if (hit.items[0] < hit.items[1] &&
                    hit.length >= min_length) {
                Block* new_block = new Block;
                add_blast_item(bs, name2seq, name2block,
                               new_block, hit.items[0]);
                add_blast_item(bs, name2seq, name2block,
                               new_block, hit.items[1]);
                block_set()->insert(new_block);
            }
        }
    }
}

const char* ImportBlastHits::name_impl() const {
    return "Import blast hits";
}

}

