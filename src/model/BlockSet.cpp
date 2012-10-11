/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "BlockSet.hpp"
#include "FastaReader.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "Connector.hpp" // FIXME
#include "po.hpp"

namespace bloomrepeats {

BlockSet::~BlockSet() {
    clear();
}

BlockSetPtr BlockSet::clone() const {
    BlockSetPtr result = boost::make_shared<BlockSet>();
    BOOST_FOREACH (Block* b, *this) {
        result->insert(b->clone());
    }
    Connector connector;
    connector.apply(result);
    return result;
}

void BlockSet::add_sequence(SequencePtr seq) {
    seqs_.push_back(seq);
}

void BlockSet::add_sequences(const std::vector<SequencePtr>& sequences) {
    BOOST_FOREACH (const SequencePtr& seq, sequences) {
        add_sequence(seq);
    }
}

SequencePtr BlockSet::seq_from_name(const std::string& name) const {
    Name2Seq::const_iterator it = name2seq_.find(name);
    if (it == name2seq_.end()) {
        BOOST_FOREACH (SequencePtr seq, seqs_) {
            name2seq_[seq->name()] = seq;
        }
        it = name2seq_.find(name);
    }
    if (it != name2seq_.end()) {
        return it->second;
    } else {
        return SequencePtr();
    }
}

Fragment* BlockSet::fragment_from_id(const std::string& id) const {
    if (id.empty()) {
        return 0;
    }
    size_t u1 = id.find('_');
    if (u1 == std::string::npos) {
        return 0;
    }
    std::string seq_name = id.substr(0, u1);
    if (seq_name.empty()) {
        return 0;
    }
    SequencePtr seq = seq_from_name(seq_name);
    if (!seq) {
        return 0;
    }
    size_t u2 = id.find('_', u1 + 1);
    if (u2 == std::string::npos) {
        return 0;
    }
    std::string begin_pos_str = id.substr(u1 + 1, u2 - u1 - 1);
    size_t begin_pos = boost::lexical_cast<size_t>(begin_pos_str);
    std::string last_pos_str = id.substr(u2 + 1);
    size_t last_pos = boost::lexical_cast<size_t>(last_pos_str);
    Fragment* f = new Fragment(seq);
    f->set_ori(begin_pos <= last_pos ? 1 : -1);
    f->set_begin_pos(begin_pos);
    f->set_last_pos(last_pos);
    return f;
}

std::string BlockSet::block_from_description(const std::string& description) {
    size_t block_pos = description.find("block=");
    if (block_pos == std::string::npos) {
        return "";
    }
    size_t block_name_start = block_pos + std::string("block=").size();
    size_t space_pos = description.find(' ', block_name_start); // or npos
    return description.substr(block_name_start, space_pos - block_name_start);
}

void BlockSet::insert(Block* block) {
#ifndef NDEBUG
    BOOST_FOREACH (Block* b, *this) {
        BOOST_ASSERT(block != b);
    }
#endif
    blocks_.insert(block);
}

void BlockSet::erase(Block* block) {
    blocks_.erase(block);
    delete block;
}

size_t BlockSet::size() const {
    return blocks_.size();
}

bool BlockSet::empty() const {
    return blocks_.empty();
}

bool BlockSet::has(const Block* block) const {
    return blocks_.find(const_cast<Block*>(block)) != blocks_.end();
}

void BlockSet::clear() {
    BOOST_FOREACH (Block* block, *this) {
        delete block;
    }
    blocks_.clear();
}

Block* BlockSet::front() const {
    return empty() ? 0 : *(begin());
}

BlockSet::iterator BlockSet::begin() {
    return blocks_.begin();
}

BlockSet::const_iterator BlockSet::begin() const {
    return blocks_.begin();
}

BlockSet::iterator BlockSet::end() {
    return blocks_.end();
}

BlockSet::const_iterator BlockSet::end() const {
    return blocks_.end();
}

static struct BlockGreater {
    bool operator()(const Block* b1, const Block* b2) const {
        return b1->size() > b2->size();
    }
} block_greater;

static void try_new_block(BlockSet& set, const Fragment& f, int ori,
                          Fragment** prev) {
    Fragment* n = f.neighbor(ori);
    Fragment* new_f = new Fragment(f.seq());
    if (ori == -1) {
        new_f->set_min_pos(n ? n->max_pos() + 1 : 0);
        new_f->set_max_pos(f.min_pos() - 1);
    } else {
        new_f->set_min_pos(f.max_pos() + 1);
        new_f->set_max_pos(n ? n->min_pos() - 1 : f.seq()->size() - 1);
    }
    if (new_f->valid()) {
        if (*prev) {
            BOOST_ASSERT(!(*new_f < **prev));
            Fragment::connect(*prev, new_f);
        }
        *prev = new_f;
        Block* block = new Block();
        block->insert(new_f);
        set.insert(block);
    } else {
        delete new_f;
    }
}

BlockSetPtr BlockSet::rest() const {
    BlockSetPtr result = boost::make_shared<BlockSet>();
    result->seqs_ = seqs_;
    std::set<Sequence*> used;
    BOOST_FOREACH (Block* block, *this) {
        BOOST_FOREACH (Fragment* f, *block) {
            Sequence* seq = f->seq();
            if (used.find(seq) == used.end()) {
                used.insert(seq);
                Fragment* prev = 0;
                while (Fragment* fr = f->neighbor(-1)) {
                    f = fr;
                }
                try_new_block(*result, *f, -1, &prev);
                while (Fragment* fr = f->neighbor(1)) {
                    f = fr;
                    try_new_block(*result, *f, -1, &prev);
                }
                try_new_block(*result, *f, 1, &prev);
            }
        }
    }
    return result;
}

void BlockSet::add_output_options(po::options_description& desc) {
    add_unique_options(desc)
    ("out-file,o", po::value<std::string>(), "output file with all blocks")
    ("out-mask", po::value<std::string>(),
     "mask of output files (${block} is replaced with block name)");
}

void BlockSet::make_output(const po::variables_map& vm) {
    if (!vm["out-mask"].empty()) {
        std::string mask = vm["out-mask"].as<std::string>();
        BOOST_ASSERT(mask.find("${block}") != std::string::npos);
        BOOST_FOREACH (Block* b, *this) {
            using namespace boost::algorithm;
            std::string path = replace_all_copy(mask, "${block}", b->name());
            std::ofstream o(path.c_str());
            o << *b << std::endl;
        }
    }
    if (!vm["out-file"].empty()) {
        std::string path = vm["out-file"].as<std::string>();
        std::ofstream o(path.c_str());
        o << *this << std::endl;
    }
    if (vm["out-file"].empty() && vm["out-mask"].empty()) {
        std::cout << *this << std::endl;
    }
}

void BlockSet::set_unique_block_names() {
    std::set<std::string> names;
    std::string null_name = Block().name(); // 0000 0000
    BOOST_FOREACH (Block* b, *this) {
        if (b->name() == null_name) {
            b->set_name_from_fragments();
        }
        while (names.find(b->name()) != names.end()) {
            b->set_random_name();
        }
        names.insert(b->name());
    }
}

class BlockSetFastaReader : public FastaReader {
public:
    BlockSetFastaReader(BlockSet& block_set, std::istream& input):
        FastaReader(input), block_set_(block_set)
    { }

    void new_sequence(const std::string& name, const std::string& description) {
        Fragment* f = block_set_.fragment_from_id(name);
        BOOST_ASSERT(f);
        std::string block_name = block_set_.block_from_description(description);
        BOOST_ASSERT(!block_name.empty());
        Block* block = name2block_[block_name];
        if (!block) {
            block = new Block(block_name);
            name2block_[block_name] = block;
            block_set_.insert(block);
        }
        block->insert(f);
    }

    void grow_sequence(const std::string& data)
    { }

private:
    BlockSet& block_set_;
    std::map<std::string, Block*> name2block_;
};

std::istream& operator>>(std::istream& input, BlockSet& block_set) {
    BlockSetFastaReader reader(block_set, input);
    reader.read_all_sequences();
    return input;
}

static struct BlockCompareName {
    bool operator()(const Block* b1, const Block* b2) const {
        typedef boost::tuple<int, const std::string&> Tie;
        return Tie(-b1->size(), b1->name()) < Tie(-b2->size(), b2->name());
    }
} bcn;

std::ostream& operator<<(std::ostream& o, const BlockSet& block_set) {
    std::vector<Block*> blocks(block_set.begin(), block_set.end());
    std::sort(blocks.begin(), blocks.end(), bcn);
    BOOST_FOREACH (Block* block, blocks) {
        o << *block;
        o << std::endl; // empty line
    }
    return o;
}

}

