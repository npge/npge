/*
 * NPG-explorer, Nucleotide PanGenome explorer
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
#include "MetaProcessor.hpp"
#include "Clear.hpp"
#include "FragmentCollection.hpp"
#include "Align.hpp"
#include "UniqueNames.hpp"
#include "RemoveNames.hpp"
#include "SizeLimits.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"
#include "convert_position.hpp"
#include "global.hpp"

namespace npge {

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
    mutable BlockSet* subblocks;

protected:
    void run_impl() const {
        BlockSet& t = *block_set();
        BlockSet& o = *other();
        s2f.clear();
        s2f.add_bs(t);
        Blocks blocks(o.begin(), o.end());
        std::sort(blocks.begin(), blocks.end(), bll);
        BlockSetPtr subblocks_ptr = new_bs();
        subblocks = subblocks_ptr.get();
        subblocks->add_sequences(o.seqs());
        BOOST_FOREACH (Block* block, blocks) {
            process_block(block);
            ASSERT_FALSE(o.has(block));
        }
        ASSERT_TRUE(o.empty());
        o.swap(*subblocks);
        s2f.clear();
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

    typedef std::vector<bool> GoodPos;
    mutable GoodPos good_pos;

    void split_block(Block* block) const {
        BlockSet& o = *other();
        int length = block->alignment_length();
        good_pos.clear();
        good_pos.resize(length, true);
        BOOST_FOREACH (Fragment* fragment, *block) {
            std::vector<Fragment> overlaps;
            s2f.find_overlaps(overlaps, fragment);
            BOOST_FOREACH (const Fragment& ol, overlaps) {
                ASSERT_NE(ol, Fragment::INVALID);
                ASSERT_GT(ol.length(), 0);
                mark_bad(ol, fragment);
            }
        }
        add_subblocks(block);
        o.erase(block);
    }

    void mark_bad(const Fragment& ol, Fragment* fragment) const {
        int seq_a = ol.min_pos();
        int seq_b = ol.max_pos();
        int fragment_a = seq_to_frag(fragment, seq_a);
        int fragment_b = seq_to_frag(fragment, seq_b);
        int length = good_pos.size();
        int block_a = block_pos(fragment, fragment_a, length);
        int block_b = block_pos(fragment, fragment_b, length);
        int min_pos = std::min(block_a, block_b);
        int max_pos = std::max(block_a, block_b);
        ASSERT_LTE(0, min_pos);
        ASSERT_LTE(min_pos, max_pos);
        ASSERT_LT(max_pos, length);
        for (int col = min_pos; col <= max_pos; col++) {
            good_pos[col] = false;
        }
    }

    void add_subblocks(Block* block) const {
        int length = good_pos.size();
        int first_good = -1;
        for (int col = 0; col <= length; col++) {
            if (good_pos[col] && first_good == -1) {
                first_good == col;
            } else if (!good_pos[col] && first_good != -1) {
                int last_good = col - 1;
                add_subblock(block, first_good, last_good);
                first_good = -1;
            }
        }
        if (first_good != -1) {
            add_subblock(block, first_good, length - 1);
        }
    }

    void add_subblock(Block* block, int min_pos, int max_pos) const {
        subblocks->insert(block->slice(min_pos, max_pos));
    }
};

class AddingLoop : public Pipe {
public:
    AddingLoop() {
        add(new Align, "target=other");
        add(new SmthUnion);
    }
};

AddingLoopBySize::AddingLoopBySize() {
    al_ = new AddingLoop;
    al_->set_parent(this);
    al_->set_options("target=target other=other", this);
    declare_bs("other", "source blockset");
    declare_bs("target", "destination blockset");
}

void AddingLoopBySize::run_impl() const {
    while (!other()->empty()) {
        al_->run();
    }
}

const char* AddingLoopBySize::name_impl() const {
    return "Align and move overlapless from other to target";
}

TrySmth::TrySmth() {
    add(new Union, "target=smth-copy other=target");
    add(new UniqueNames, "target=smth-copy");
    add(new MetaProcessor, "prefix|smth-");
    add(new RemoveNames, "target=target --remove-seqs-names:=0 "
        " --remove-blocks-names:=1");
    add(new Align, "target=target");
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

