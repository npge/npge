/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

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
#include "Filter.hpp"
#include "UniqueNames.hpp"
#include "RemoveNames.hpp"
#include "SizeLimits.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "block_hash.hpp"
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
    mutable S2F s2f;
    mutable std::set<Block*> o_b;
    BlockLengthLess bll;

protected:
    void run_impl() const {
        BlockSet& t = *block_set();
        BlockSet& o = *other();
        s2f.clear();
        s2f.add_bs(t);
        Blocks blocks(o.begin(), o.end());
        std::sort(blocks.begin(), blocks.end(), bll);
        BOOST_FOREACH (Block* block, blocks) {
            process_block(block);
        }
    }

    void process_block(Block* block) const {
        bool swap, remove;
        test_block(block, swap, remove);
        if (swap) {
            swap_block(block);
        } else if (remove) {
            remove_block(block);
        }
    }

    void test_block(Block* block, bool& swap, bool& remove) const {
        std::vector<Fragment*> o_f;
        BOOST_FOREACH (Fragment* fragment, *block) {
            s2f.find_overlap_fragments(o_f, fragment);
        }
        o_b.clear();
        BOOST_FOREACH (Fragment* f, o_f) {
            o_b.insert(f->block());
        }
        swap = true;
        remove = false;
        BOOST_FOREACH (Block* b, o_b) {
            if (!has_alignment(b)) {
                swap = false;
                break;
            }
            if (!bll(block, b)) {
                remove = true;
                swap = false;
                break;
            }
        }
    }

    void swap_block(Block* block) const {
        BlockSet& t = *block_set();
        BlockSet& o = *other();
        BOOST_FOREACH (Block* b, o_b) {
            s2f.remove_block(b);
            t.detach(b);
            o.insert(b);
        }
        s2f.add_block(block);
        o.detach(block);
        t.insert(block);
    }

    void remove_block(Block* block) const {
        BlockSet& t = *block_set();
        BlockSet& o = *other();
        o.erase(block);
    }
};

class AddingLoop : public Pipe {
public:
    AddingLoop() {
        add(new SmthUnion);
        add(new Align);
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

