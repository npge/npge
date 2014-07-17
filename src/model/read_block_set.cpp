/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "read_block_set.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "complement.hpp"
#include "FastaReader.hpp"
#include "key_value.hpp"
#include "to_s.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"
#include "thread_pool.hpp"
#include "global.hpp"

namespace npge {

typedef std::vector<BlockSet*> BlockSets;
typedef std::map<std::string, BlockSet*> Name2BlockSet;
typedef std::map<std::string, Block*> Name2Block;
typedef std::map<BlockSet*, Name2Block> Bs2Name2Block;
typedef std::map<std::string, SequencePtr> Name2Seq;

struct BSFRImpl {
    Name2BlockSet name2block_set_;
    Name2Seq name2seq_;
    std::istream* input_ptr_;
    RowType row_type_;
    SequenceType seq_type_;
    int workers_;
    bool unknown_bs_allowed_;

    BSFRImpl():
        workers_(1),
        unknown_bs_allowed_(true) {
    }
};

struct BlockSetFastaReader::Impl : public BSFRImpl {
};

BlockSetFastaReader::BlockSetFastaReader(BlockSet& block_set,
        std::istream& input, RowType row_type,
        SequenceType seq_type):
    impl_(new Impl) {
    impl_->input_ptr_ = &input;
    impl_->row_type_ = row_type;
    impl_->seq_type_ = seq_type;
    set_block_set("target", &block_set);
}

BlockSetFastaReader::~BlockSetFastaReader() {
    delete impl_;
}

void BlockSetFastaReader::set_block_set(const std::string& name,
                                        BlockSet* block_set) {
    impl_->name2block_set_[name] = block_set;
    BOOST_FOREACH (SequencePtr seq, block_set->seqs()) {
        impl_->name2seq_[seq->name()] = seq;
    }
}

static BlockSet* get_blockset(const std::string& name,
                              BSFRImpl* impl) {
    const Name2BlockSet& n2bs = impl->name2block_set_;
    Name2BlockSet::const_iterator it = n2bs.find(name);
    if (it != n2bs.end()) {
        return it->second;
    } else {
        return 0;
    }
}

BlockSet* BlockSetFastaReader::get_block_set(
    const std::string& name) const {
    return get_blockset(name, impl_);
}

bool BlockSetFastaReader::unknown_bs_allowed() const {
    return impl_->unknown_bs_allowed_;
}

void BlockSetFastaReader::set_unknown_bs_allowed(
    bool unknown_bs_allowed) {
    impl_->unknown_bs_allowed_ = unknown_bs_allowed;
}

int BlockSetFastaReader::workers() const {
    return impl_->workers_;
}

void BlockSetFastaReader::set_workers(int workers) {
    impl_->workers_ = workers;
}

// fasta reader

struct FastaValue {
    std::string data_;
    std::string descr_;
    std::string block_name_;
    SequencePtr s_;
    Fragment* f_;
    BlockSets bss_;

    FastaValue():
        f_(0) {
    }
};

typedef std::pair<std::string, FastaValue> FastaItem;
typedef std::vector<FastaItem> FastaMap;

class SimpleReader : public FastaReader {
public:
    FastaMap sequences_;
    FastaMap fragments_;

    SimpleReader(std::istream& input):
        FastaReader(input),
        v_(0) {
    }

    void new_sequence(const std::string& name,
                      const std::string& description) {
        bool is_fr = (name.find('_') != std::string::npos);
        FastaMap& map = is_fr ? fragments_ : sequences_;
        map.push_back(FastaItem(name, FastaValue()));
        v_ = &(map.back().second);
        v_->descr_ = description;
    }

    void grow_sequence(const std::string& data) {
        ASSERT_TRUE(v_);
        v_->data_ += data;
    }

private:
    FastaValue* v_;
};

// find blocksets list by fasta description

typedef std::map<std::string, BlockSets> Name2BlockSets;

static void checked_add(BlockSets& bss, BSFRImpl* impl,
                        const std::string& name) {
    BlockSet* bs = get_blockset(name, impl);
    if (bs) {
        bss.push_back(bs);
    } else if (!impl->unknown_bs_allowed_) {
        throw Exception("Unknown block set '" + name + "'");
    }
}

static const BlockSets& name2bss(const std::string& name,
                                 BSFRImpl* impl,
                                 Name2BlockSets& cache) {
    Name2BlockSets::iterator it = cache.find(name);
    if (it != cache.end()) {
        return it->second;
    } else {
        BlockSets& bss = cache[name];
        if (name.empty()) {
            checked_add(bss, impl, "target");
        } else if (name == "all") {
            typedef Name2BlockSet::value_type N2BSI;
            const Name2BlockSet& n2bs = impl->name2block_set_;
            BOOST_FOREACH (const N2BSI& n2bsi, n2bs) {
                BlockSet* bs = n2bsi.second;
                bss.push_back(bs);
            }
        } else {
            Strings names;
            using namespace boost::algorithm;
            split(names, name, is_any_of(","));
            BOOST_FOREACH (const std::string& n, names) {
                checked_add(bss, impl, n);
            }
        }
        return bss;
    }
}

class BSWorker : public ThreadWorker {
public:
    Name2BlockSets cache_;

    BSWorker(ThreadGroup* g):
        ThreadWorker(g) {
    }
};

class BSTG : public ReusingThreadGroup {
public:
    FastaMap& map_;
    FastaMap::iterator it_;
    BSFRImpl* impl_;

    BSTG(FastaMap& map, BSFRImpl* impl):
        map_(map), it_(map.begin()), impl_(impl) {
        set_workers(impl->workers_);
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl() {
        return new BSWorker(this);
    }
};

class BSTask : public ThreadTask {
public:
    BSTask(FastaItem* item, ThreadWorker* worker):
        ThreadTask(worker), item_(item) {
    }

    void run_impl() {
        FastaValue& v = item_->second;
        BSTG* g = D_CAST<BSTG*>(thread_group());
        BSWorker* w = D_CAST<BSWorker*>(worker());
        std::string sets = extract_value(v.descr_, "set");
        v.bss_ = name2bss(sets, g->impl_, w->cache_);
    }

private:
    FastaItem* item_;
};

ThreadTask* BSTG::create_task_impl(ThreadWorker* worker) {
    if (it_ != map_.end()) {
        FastaItem* item = &(*it_);
        it_++;
        return new BSTask(item, worker);
    } else {
        return 0;
    }
}

// full sequences

class STG : public ReusingThreadGroup {
public:
    FastaMap& map_;
    FastaMap::iterator it_;
    BSFRImpl* impl_;
    SequenceType type_;

    STG(FastaMap& map, BSFRImpl* impl):
        map_(map), it_(map.begin()),
        impl_(impl), type_(impl->seq_type_) {
        set_workers(impl->workers_);
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);
};

class SequenceCreator : public ThreadTask {
public:
    SequenceCreator(FastaItem* item, ThreadWorker* worker):
        ThreadTask(worker), item_(item) {
    }

    void run_impl() {
        STG* g = D_CAST<STG*>(thread_group());
        SequenceType type = g->type_;
        SequencePtr s = Sequence::new_sequence(type);
        const std::string& name = item_->first;
        FastaValue& v = item_->second;
        v.s_ = s;
        s->read_from_string(v.data_);
        s->set_name(name);
        s->set_description(v.descr_);
    }

private:
    FastaItem* item_;
};

ThreadTask* STG::create_task_impl(ThreadWorker* worker) {
    if (it_ != map_.end()) {
        FastaItem* item = &(*it_);
        it_++;
        return new SequenceCreator(item, worker);
    } else {
        return 0;
    }
}

static void add_sequences(const FastaMap& sequences,
                          BSFRImpl* impl) {
    Name2Seq& name2seq = impl->name2seq_;
    BOOST_FOREACH (const FastaItem& item, sequences) {
        const std::string& name = item.first;
        const FastaValue& v = item.second;
        const SequencePtr& s = v.s_;
        BOOST_FOREACH (BlockSet* bs, v.bss_) {
            bs->add_sequence(s);
        }
        name2seq[name] = s;
    }
}

static void add_sequences_from_fragments(
    const FastaMap& fragments,
    BSFRImpl* impl) {
    Name2Seq& name2seq = impl->name2seq_;
    SequenceType type = impl->seq_type_;
    BOOST_FOREACH (const FastaItem& item, fragments) {
        const std::string& name = item.first;
        Strings parts;
        using namespace boost::algorithm;
        split(parts, name, is_any_of("_"));
        const std::string& seq_name = parts[0];
        if (name2seq.find(seq_name) == name2seq.end()) {
            const FastaValue& v = item.second;
            SequencePtr s = Sequence::new_sequence(type);
            s->set_name(seq_name);
            name2seq[seq_name] = s;
            BOOST_FOREACH (BlockSet* bs, v.bss_) {
                bs->add_sequence(s);
            }
        }
    }
}

// fragments

class FTG : public ReusingThreadGroup {
public:
    FastaMap& map_;
    FastaMap::iterator it_;
    BSFRImpl* impl_;
    RowType type_;

    FTG(FastaMap& map, BSFRImpl* impl):
        map_(map), it_(map.begin()),
        impl_(impl), type_(impl->row_type_) {
        set_workers(impl->workers_);
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);
};

static bool has_norow(const std::string& description) {
    std::string spaced_d = " " + description + " ";
    return spaced_d.find("norow") != std::string::npos;
}

class FragmentCreator : public ThreadTask {
public:
    FragmentCreator(FastaItem* item, ThreadWorker* worker):
        ThreadTask(worker), item_(item) {
    }

    void run_impl() {
        FTG* g = D_CAST<FTG*>(thread_group());
        const std::string& name = item_->first;
        FastaValue& v = item_->second;
        Strings parts;
        using namespace boost::algorithm;
        split(parts, name, is_any_of("_"));
        const std::string& seq_name = parts[0];
        int begin_pos = L_CAST<int>(parts[1]);
        int last_pos = L_CAST<int>(parts[2]);
        BSFRImpl* impl = g->impl_;
        const Name2Seq& name2seq = impl->name2seq_;
        Name2Seq::const_iterator it = name2seq.find(seq_name);
        ASSERT_MSG(it != name2seq.end(), seq_name.c_str());
        const SequencePtr& s = it->second;
        Fragment* f = new Fragment(s);
        v.f_ = f;
        f->set_ori(begin_pos <= last_pos ? 1 : -1);
        f->set_begin_pos(begin_pos);
        f->set_last_pos(last_pos);
        v.block_name_ = extract_value(v.descr_, "block");
        // row
        bool norow = has_norow(v.descr_);
        if (!norow) {
            RowType type = g->type_;
            AlignmentRow* row = AlignmentRow::new_row(type);
            if (s->size() > 0) {
                // check contents of row equals
                // corresponding letters from sequence
                f->set_row(row);
                row->grow(v.data_);
            } else {
                // do not check
                // sequence is built from fragments
                row->grow(v.data_);
                f->set_row(row);
            }
        }
    }

private:
    FastaItem* item_;
};

ThreadTask* FTG::create_task_impl(ThreadWorker* worker) {
    if (it_ != map_.end()) {
        FastaItem* item = &(*it_);
        it_++;
        return new FragmentCreator(item, worker);
    } else {
        return 0;
    }
}

static Block* get_block(const std::string& name,
                        BlockSet* bs, Bs2Name2Block& bs2n2b) {
    Block*& block = bs2n2b[bs][name];
    if (!block) {
        block = new Block;
        block->set_name(name);
        bs->insert(block);
    }
    return block;
}

static void add_blocks(const FastaMap& fragments,
                       BSFRImpl* impl) {
    Bs2Name2Block b2;
    BOOST_FOREACH (const FastaItem& item, fragments) {
        const FastaValue& v = item.second;
        Fragment* f = v.f_;
        const std::string& block_name = v.block_name_;
        const BlockSets& bss = v.bss_;
        if (bss.empty()) {
            // strange: useless fragment
            delete f;
        } else {
            Block* block = get_block(block_name, bss[0], b2);
            block->insert(f);
            // now copy to other blocksets, if any
            for (int i = 1; i < bss.size(); i++) {
                block = get_block(block_name, bss[i], b2);
                block->insert(f->clone());
            }
        }
    }
}

// sequences from fragments

typedef std::vector<FastaValue*> FastaValues;
typedef std::map<Sequence*, FastaValues> S2FV;
typedef S2FV::value_type SItem;

static void make_s2fv(S2FV& result,
                      FastaMap& fragments) {
    BOOST_FOREACH (FastaItem& item, fragments) {
        FastaValue& v = item.second;
        Fragment* f = v.f_;
        Sequence* s = f->seq();
        if (s->size() == 0) {
            result[s].push_back(&v);
        }
    }
}

class S2TG : public ReusingThreadGroup {
public:
    S2FV& map_;
    S2FV::iterator it_;

    S2TG(S2FV& map, BSFRImpl* impl):
        map_(map), it_(map.begin()) {
        set_workers(impl->workers_);
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);
};

// FIXME DRY
struct SimpleFragmentCmp {
    bool operator()(FastaValue* x, FastaValue* y) {
        Fragment* a = x->f_;
        Fragment* b = y->f_;
        typedef boost::tuple<int, int> Tie;
        return Tie(a->min_pos(), a->max_pos()) <
               Tie(b->min_pos(), b->max_pos());
    }
};

class SequenceTrace : public ThreadTask {
public:
    SequenceTrace(SItem* item, ThreadWorker* worker):
        ThreadTask(worker), item_(item) {
    }

    void run_impl() {
        Sequence* s = item_->first;
        FastaValues& values = item_->second;
        // sort
        std::sort(values.begin(), values.end(),
                  SimpleFragmentCmp());
        // grow string
        std::string data;
        BOOST_FOREACH (FastaValue* v, values) {
            Fragment* f = v->f_;
            ASSERT_EQ(f->seq(), s);
            if (f->max_pos() < data.size()) {
                continue;
            }
            std::string& d = v->data_;
            Sequence::to_atgcn(d);
            if (f->ori() == -1) {
                complement(d);
            }
            if (f->min_pos() > data.size()) {
                // gap
                int gap_length = f->min_pos() - data.size();
                std::string gap(gap_length, 'N');
                data += gap;
            }
            const char* part = d.c_str();
            if (f->min_pos() < data.size()) {
                int overlap = data.size() - f->min_pos();
                ASSERT_LT(overlap, d.size());
                part += overlap;
            }
            data += part;
            ASSERT_EQ(data.length(), f->max_pos() + 1);
        }
        // set sequence
        s->read_from_string(data);
    }

private:
    SItem* item_;
};

ThreadTask* S2TG::create_task_impl(ThreadWorker* worker) {
    if (it_ != map_.end()) {
        SItem* item = &(*it_);
        it_++;
        return new SequenceTrace(item, worker);
    } else {
        return 0;
    }
}

// main function

void BlockSetFastaReader::run() {
    SimpleReader reader(*impl_->input_ptr_);
    reader.read_all_sequences();
    {
        BSTG bstg(reader.sequences_, impl_);
        bstg.perform();
    }
    {
        BSTG bstg(reader.fragments_, impl_);
        bstg.perform();
    }
    {
        STG stg(reader.sequences_, impl_);
        stg.perform();
    }
    add_sequences(reader.sequences_, impl_);
    reader.sequences_.clear();
    add_sequences_from_fragments(reader.fragments_, impl_);
    {
        FTG ftg(reader.fragments_, impl_);
        ftg.perform();
    }
    add_blocks(reader.fragments_, impl_);
    S2FV s2fv;
    make_s2fv(s2fv, reader.fragments_);
    {
        S2TG s2tg(s2fv, impl_);
        s2tg.perform();
    }
}

}

