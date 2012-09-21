/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "BlockSet.hpp"
#include "FastaReader.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "PairAligner.hpp"
#include "JoinApprover.hpp"
#include "po.hpp"

namespace bloomrepeats {

BlockSet::~BlockSet() {
    clear();
}

BlockSetPtr BlockSet::clone() const {
    BlockSetPtr result = boost::make_shared<BlockSet>();
    BOOST_FOREACH (BlockPtr b, *this) {
        result->insert(b->clone());
    }
    result->connect_fragments();
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

void BlockSet::insert(BlockPtr block) {
#ifndef NDEBUG
    BOOST_FOREACH (BlockPtr b, *this) {
        BOOST_ASSERT(block != b);
    }
#endif
    blocks_.insert(block);
}

void BlockSet::erase(BlockPtr block) {
    blocks_.erase(block);
    delete block;
}

size_t BlockSet::size() const {
    return blocks_.size();
}

bool BlockSet::empty() const {
    return blocks_.empty();
}

bool BlockSet::has(BlockPtr block) const {
    return blocks_.find(block) != blocks_.end();
}

void BlockSet::clear() {
    BOOST_FOREACH (BlockPtr block, *this) {
        delete block;
    }
    blocks_.clear();
}

BlockPtr BlockSet::front() const {
    return empty() ? BlockPtr() : *(begin());
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

static struct FragmentCompare {
    bool operator()(const FragmentPtr& f1, const FragmentPtr& f2) const {
        return *f1 < *f2;
    }
} fragment_compare;

void BlockSet::connect_fragments() {
    typedef std::vector<FragmentPtr> Fs;
    typedef std::map<Sequence*, Fs> Seq2Fs;
    Seq2Fs seq2fs;
    BOOST_FOREACH (BlockPtr block, *this) {
        BOOST_FOREACH (FragmentPtr fragment, *block) {
            seq2fs[fragment->seq()].push_back(fragment);
        }
    }
    BOOST_FOREACH (Seq2Fs::value_type& seq_and_fs, seq2fs) {
        Fs& fs = seq_and_fs.second;
        std::sort(fs.begin(), fs.end(), fragment_compare);
        for (int i = 1; i < fs.size(); i++) {
            Fragment::connect(fs[i - 1], fs[i]);
        }
    }
}

void BlockSet::filter(int min_fragment_length, int min_block_size) {
    std::vector<BlockPtr> block_set_copy(begin(), end());
    BOOST_FOREACH (BlockPtr block, block_set_copy) {
        block->filter(min_fragment_length);
        if (block->size() < min_block_size) {
            erase(block);
        }
    }
}

static struct BlockGreater {
    bool operator()(const BlockPtr& b1, const BlockPtr& b2) const {
        return b1->size() > b2->size();
    }
} block_greater;

static BlockPtr neighbor_block(BlockPtr b, int ori) {
    BlockPtr result = 0;
    FragmentPtr f = b->front();
    if (f) {
        FragmentPtr neighbor_f = ori == 1 ? f->next() : f->prev();
        if (neighbor_f) {
            result = neighbor_f->block();
        }
    }
    return result;
}

void BlockSet::join(JoinApprover* j) {
    std::vector<BlockPtr> bs(begin(), end());
    std::sort(bs.begin(), bs.end(), block_greater);
    BOOST_FOREACH (BlockPtr block, bs) {
        if (has(block)) {
            for (int ori = -1; ori <= 1; ori += 2) {
                while (BlockPtr other_block = neighbor_block(block, ori)) {
                    BlockPtr new_block = Block::try_join(block, other_block, j);
                    if (new_block) {
                        erase(block);
                        erase(other_block);
                        insert(new_block);
                        block = new_block;
                    } else {
                        break;
                    }
                }
            }
        }
    }
}

void BlockSet::expand_blocks(PairAligner* aligner, int batch,
                             int ori, int max_overlap) {
    aligner = aligner ? : PairAligner::default_aligner();
    std::vector<BlockPtr> bs(begin(), end());
    std::sort(bs.begin(), bs.end(), block_greater);
    BOOST_FOREACH (BlockPtr block, bs) {
        block->expand(aligner, batch, ori, max_overlap);
    }
}

bool BlockSet::overlaps() const {
    BOOST_FOREACH (BlockPtr block, *this) {
        BOOST_FOREACH (FragmentPtr fragment, *block) {
            for (int ori = -1; ori <= 1; ori += 2) {
                FragmentPtr neighbor = fragment->neighbor(ori);
                if (neighbor && fragment->common_positions(*neighbor)) {
                    return true;
                }
            }
        }
    }
    return false;
}

struct BlockLess {
    BlockLess(BlockSet* block_set):
        block_set_(block_set)
    { }

    bool operator()(const BlockPtr& b1, const BlockPtr& b2) const {
        return (block_set_->has(b1) && block_set_->has(b2)) ?
               b1->size() < b2->size() : false;
    }

    BlockSet* block_set_;
};

typedef std::priority_queue<BlockPtr, std::vector<BlockPtr>, BlockLess> BQ;

static void treat_fragments(BlockSet* block_set, BQ& bs,
                            FragmentPtr x, FragmentPtr y) {
    BlockPtr x_block = x->block();
    BlockPtr y_block = y->block();
    if (x_block == y_block) {
        x_block->erase(x);
        return;
    }
    Fragment common = x->common_fragment(*y);
    BOOST_ASSERT(common.valid());
    if (*x == common && x->length() == y->length()) {
        BOOST_ASSERT(y_block);
        x->block()->merge(y_block);
        BOOST_ASSERT(y_block->empty());
        block_set->erase(y_block);
    } else {
        if (common == *x) {
            treat_fragments(block_set, bs, y, x);
        } else {
            size_t new_length;
            if (common.begin_pos() == x->begin_pos()) {
                new_length = common.length();
            } else {
                new_length = std::min(abs(x->begin_pos() - common.min_pos()),
                                      abs(x->begin_pos() - common.max_pos()));
            }
            BlockPtr new_block = x->block()->split(new_length);
            BOOST_ASSERT(new_block && !new_block->empty());
            bs.push(new_block);
            block_set->insert(new_block);
        }
    }
}

static bool treat_block(BlockSet* block_set, BQ& bs, BlockPtr block) {
    BOOST_FOREACH (FragmentPtr f, *block) {
        for (int ori = -1; ori <= 1; ori += 2) {
            FragmentPtr o_f = f->neighbor(ori);
            if (o_f && f->common_positions(*o_f)) {
                treat_fragments(block_set, bs, f, o_f);
                return true;
            }
        }
    }
    return false;
}

void BlockSet::resolve_overlaps() {
    BQ bs(begin(), end(), BlockLess(this));
    while (!bs.empty()) {
        BlockPtr block = bs.top();
        bs.pop();
        while (has(block) && treat_block(this, bs, block))
        { }
    }
#ifndef NDEBUG
    BOOST_ASSERT(!overlaps());
    connect_fragments();
    BOOST_ASSERT(!overlaps());
#endif
}

bool BlockSet::expand_blocks_by_fragments(PairAligner* aligner, int batch) {
    aligner = aligner ? : PairAligner::default_aligner();
    bool result = false;
    BOOST_FOREACH (BlockPtr block, *this) {
        BOOST_ASSERT(block);
        result |= block->expand_by_fragments(aligner, batch);
    }
    return result;
}

static void try_new_block(BlockSet& set, const Fragment& f, int ori,
                          FragmentPtr* prev) {
    FragmentPtr n = f.neighbor(ori);
    FragmentPtr new_f = Fragment::create_new(f.seq());
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
        BlockPtr block = Block::create_new();
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
    BOOST_FOREACH (BlockPtr block, *this) {
        BOOST_FOREACH (FragmentPtr f, *block) {
            Sequence* seq = f->seq();
            if (used.find(seq) == used.end()) {
                used.insert(seq);
                FragmentPtr prev = 0;
                while (FragmentPtr fr = f->neighbor(-1)) {
                    f = fr;
                }
                try_new_block(*result, *f, -1, &prev);
                while (FragmentPtr fr = f->neighbor(1)) {
                    f = fr;
                    try_new_block(*result, *f, -1, &prev);
                }
                try_new_block(*result, *f, 1, &prev);
            }
        }
    }
    return result;
}

void BlockSet::add_pangenome_options(po::options_description& desc) {
    // TODO
}

void BlockSet::make_pangenome(const po::variables_map& vm) {
    connect_fragments();
    resolve_overlaps();
    join(0);
    filter(10);
    expand_blocks_by_fragments();
    expand_blocks();
    filter(100);
    JoinApprover dist_1000(1000);
    join(&dist_1000);
}

void BlockSet::add_output_options(po::options_description& desc) {
    desc.add_options()
    ("out-file,o", po::value<std::string>(), "output file with all blocks")
    ("out-mask", po::value<std::string>(),
     "mask of output files (${block} is replaced with block name)");
}

void BlockSet::make_output(const po::variables_map& vm) {
    if (!vm["out-mask"].empty()) {
        std::string mask = vm["out-mask"].as<std::string>();
        BOOST_ASSERT(mask.find("${block}") != std::string::npos);
        BOOST_FOREACH (BlockPtr b, *this) {
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

class BlockSetFastaReader : public FastaReader {
public:
    BlockSetFastaReader(BlockSet& block_set, std::istream& input):
        FastaReader(input), block_set_(block_set) {
        BOOST_FOREACH (SequencePtr seq, block_set_.seqs_) {
            name2seq_[seq->name()] = seq.get();
        }
    }

    void new_sequence(const std::string& name, const std::string& description) {
        BOOST_ASSERT(!name.empty());
        size_t u1 = name.find('_');
        BOOST_ASSERT(u1 != std::string::npos);
        std::string seq_name = name.substr(0, u1);
        Sequence* seq = name2seq_[seq_name];
        BOOST_ASSERT(seq);
        BOOST_ASSERT(!seq_name.empty());
        size_t u2 = name.find('_', u1 + 1);
        BOOST_ASSERT(u2 != std::string::npos);
        std::string begin_pos_str = name.substr(u1 + 1, u2 - u1 - 1);
        size_t begin_pos = boost::lexical_cast<size_t>(begin_pos_str);
        std::string last_pos_str = name.substr(u2 + 1);
        size_t last_pos = boost::lexical_cast<size_t>(last_pos_str);
        FragmentPtr f = Fragment::create_new(seq);
        f->set_ori(begin_pos < last_pos ? 1 : -1);
        f->set_begin_pos(begin_pos);
        f->set_last_pos(last_pos);
        // block name
        size_t block_pos = description.find("block=");
        BOOST_ASSERT(block_pos != std::string::npos);
        size_t block_name_start = block_pos + std::string("block=").size();
        size_t space_pos = description.find(' ', block_name_start); // or npos
        std::string block_name = description.substr(block_name_start,
                                 space_pos - block_name_start);
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
    std::map<std::string, Sequence*> name2seq_;
    std::map<std::string, Block*> name2block_;
};

std::istream& operator>>(std::istream& input, BlockSet& block_set) {
    BlockSetFastaReader reader(block_set, input);
    reader.read_all_sequences();
    return input;
}

std::ostream& operator<<(std::ostream& o, const BlockSet& block_set) {
    BOOST_FOREACH (BlockPtr block, block_set) {
        o << *block;
        o << std::endl; // empty line
    }
    return o;
}

}

