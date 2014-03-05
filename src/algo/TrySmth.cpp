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

namespace bloomrepeats {

typedef std::set<Fragment*, FragmentCompare> FragmentsSet;
typedef FragmentCollection<Fragment*, FragmentsSet> S2F;
typedef std::vector<Block*> Blocks;

struct BlockLengthLess {
    bool operator()(Block* a, Block* b) const {
        typedef boost::tuple<int, int, const std::string&> Tie;
        return Tie(b->size(), b->alignment_length(), b->name())
               <
               Tie(a->size(), a->alignment_length(), a->name());
    }
};

static bool has_alignment(const Block* block) {
    BOOST_FOREACH (const Fragment* fragment, *block) {
        if (!fragment->row()) {
            return false;
        }
    }
    return true;
}

class SmthUnion : public Processor {
protected:
    bool run_impl() const {
        BlockSet& t = *block_set();
        BlockSet& o = *other();
        S2F s2f;
        s2f.add_bs(t);
        Blocks blocks(o.begin(), o.end());
        BlockLengthLess bll;
        std::sort(blocks.begin(), blocks.end(), bll);
        BOOST_FOREACH (Block* block, blocks) {
            std::vector<Fragment*> o_f;
            BOOST_FOREACH (Fragment* fragment, *block) {
                s2f.find_overlap_fragments(o_f, fragment);
            }
            std::set<Block*> o_b;
            BOOST_FOREACH (Fragment* f, o_f) {
                o_b.insert(f->block());
            }
            bool swap = true;
            bool remove = false;
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
            if (swap) {
                BOOST_FOREACH (Block* b, o_b) {
                    s2f.remove_block(b);
                    t.detach(b);
                    o.insert(b);
                }
                s2f.add_block(block);
                o.detach(block);
                t.insert(block);
            } else if (remove) {
                o.erase(block);
            }
        }
        return true;
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
    bool run_impl() const {
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
    add(new Clear, "target=smth-copy --clear-seqs:=1");
}

}

