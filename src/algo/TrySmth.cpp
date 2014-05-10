/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <map>
#include <set>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "TrySmth.hpp"
#include "Union.hpp"
#include "Move.hpp"
#include "Rest.hpp"
#include "DeConSeq.hpp"
#include "MetaProcessor.hpp"
#include "Clear.hpp"
#include "FragmentCollection.hpp"
#include "Align.hpp"
#include "Filter.hpp"
#include "UniqueNames.hpp"
#include "RemoveNames.hpp"
#include "SizeLimits.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Sequence.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"
#include "convert_position.hpp"
#include "global.hpp"

namespace bloomrepeats {

typedef std::set<Fragment*, FragmentCompare> FragmentsSet;
typedef FragmentCollection<Fragment*, FragmentsSet> S2F;

struct BlockLengthLess {
    bool operator()(Block* a, Block* b) const {
        typedef boost::tuple<int, int, const std::string&> Tie;
        return Tie(b->size(), b->alignment_length(), b->name())
               <
               Tie(a->size(), a->alignment_length(), a->name());
    }
};

class SmthUnion : public Processor {
private:
    BlockLengthLess bll;
    mutable S2F s2f;
    mutable BlockSetPtr cons_bs, rest_bs, subblocks;
    mutable Rest rest;
    mutable DeConSeq deconseq;
    mutable SequencePtr cons;

protected:
    void before_run() const {
        BlockSet& t = *block_set();
        BlockSet& o = *other();
        s2f.clear();
        s2f.add_bs(t);
        cons.reset(new DummySequence);
        rest.set_bs("other", new_bs());
        rest.set_bs("target", new_bs());
        cons_bs = rest.other();
        rest_bs = rest.block_set();
        deconseq.set_bs("other", rest_bs);
        deconseq.set_bs("target", new_bs());
        subblocks = deconseq.block_set();
        subblocks->add_sequences(o.seqs());
        cons_bs->add_sequence(cons);
        rest_bs->add_sequence(cons);
    }

    void run_impl() const {
        before_run();
        BlockSet& o = *other();
        Blocks blocks(o.begin(), o.end());
        std::sort(blocks.begin(), blocks.end(), bll);
        BOOST_FOREACH (Block* block, blocks) {
            process_block(block);
            ASSERT_FALSE(o.has(block));
        }
        ASSERT_TRUE(o.empty());
        o.swap(*subblocks);
        after_run();
    }

    void after_run() const {
        s2f.clear();
        cons_bs->clear();
        rest_bs->clear();
        subblocks->clear();
        cons.reset();
    }

    void process_block(Block* block) const {
        if (!s2f.block_has_overlap(block)) {
            move_block(block);
        } else {
            split_block(block);
        }
    }

    void move_block(Block* block) const {
        BlockSet& t = *block_set();
        BlockSet& o = *other();
        s2f.add_block(block);
        o.detach(block);
        t.insert(block);
    }

    void split_block(Block* block) const {
        BlockSet& o = *other();
        int length = block->alignment_length();
        cons->set_size(length);
        cons->set_block(block, /* set_consensus */ false);
        Block* cons_b = new Block;
        cons_bs->insert(cons_b);
        Sequence* cons_ptr = cons.get();
        BOOST_FOREACH (Fragment* fragment, *block) {
            std::vector<Fragment> overlaps;
            s2f.find_overlaps(overlaps, fragment);
            BOOST_FOREACH (const Fragment& ol, overlaps) {
                ASSERT_NE(ol, Fragment::INVALID);
                ASSERT_GT(ol.length(), 0);
                mark_bad(ol, fragment, cons_b, cons_ptr);
            }
        }
        rest.run();
        deconseq.run();
        o.erase(block);
        cons_bs->clear_blocks();
        rest_bs->clear_blocks();
    }

    void mark_bad(const Fragment& ol,
                  Fragment* fragment,
                  Block* cons_b,
                  Sequence* cons) const {
        int seq_a = ol.min_pos();
        int seq_b = ol.max_pos();
        int fragment_a = seq_to_frag(fragment, seq_a);
        int fragment_b = seq_to_frag(fragment, seq_b);
        int length = cons->size();
        int block_a = block_pos(fragment, fragment_a, length);
        int block_b = block_pos(fragment, fragment_b, length);
        int min_pos = std::min(block_a, block_b);
        int max_pos = std::max(block_a, block_b);
        ASSERT_LTE(0, min_pos);
        ASSERT_LTE(min_pos, max_pos);
        ASSERT_LT(max_pos, length);
        cons_b->insert(new Fragment(cons, min_pos, max_pos));
    }
};

class AddingLoop : public Pipe {
public:
    AddingLoop() {
        add(new Align, "target=other");
        add(new SmthUnion);
    }
};

class AddingLoopBySize : public Processor {
public:
    AddingLoopBySize() {
        al_ = new AddingLoop;
        al_->set_parent(this);
        al_->set_options("target=target other=other", this);
    }

protected:
    void run_impl() const {
        while (!other()->empty()) {
            al_->run();
        }
    }

private:
    AddingLoop* al_;
};

TrySmth::TrySmth() {
    add(new Union, "target=smth-copy other=target");
    add(new UniqueNames, "target=smth-copy");
    add(new MetaProcessor, "prefix|smth-");
    add(new RemoveNames, "target=target --remove-seqs-names:=0 "
        " --remove-blocks-names:=1");
    add(new Filter, "target=target");
    add(new Move, "target=smth-copy other=target");
    add(new AddingLoopBySize, "target=target other=smth-copy");
    add(new UniqueNames, "target=target");
    add(new Clear, "target=smth-copy --clear-seqs:=1 no_options");
    declare_bs("target", "Target blockset");
}

const char* TrySmth::name_impl() const {
    return "Try to do something, align (and filter), "
           "restore original if worse";
}

}

