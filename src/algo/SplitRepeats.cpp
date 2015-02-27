/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "SplitRepeats.hpp"
#include "block_stat.hpp"
#include "char_to_size.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"

namespace npge {

SplitRepeats::SplitRepeats():
    BlocksJobs("other") {
    declare_bs("target", "Destination for weak blocks");
    declare_bs("other", "Input blocks (const)");
}

class SplitRepeatsData : public ThreadData {
public:
    Blocks blocks_to_insert;
};

ThreadData* SplitRepeats::before_thread_impl() const {
    return new SplitRepeatsData;
}


static bool has_repeats(const Fragments& fragments) {
    typedef std::set<std::string> StringSet;
    StringSet genomes;
    BOOST_FOREACH (Fragment* f, fragments) {
        ASSERT_TRUE(f->seq());
        std::string genome = f->seq()->genome();
        if (genomes.find(genome) != genomes.end()) {
            return true;
        }
        genomes.insert(genome);
    }
    return false;
}

typedef std::vector<int> Ints;

// returns string of 0 and 1 starting with 0
void fragmentsToClade(std::string& clade,
                      int pos, const Fragments& all) {
    ASSERT_GTE(all.size(), 2);
    clade.resize(all.size(), '0');
    char first_letter = all[0]->alignment_at(pos);
    for (int i = 0; i < all.size(); i++) {
        Fragment* f = all[i];
        char c = f->alignment_at(pos);
        clade[i] = (c == first_letter) ? '0' : '1';
    }
}

// 1 - ident, 2 - 2 variants, or 0
void buildStatus(Ints& status, const Fragments& all) {
    int length = status.size();
    for (int pos = 0; pos < length; pos++) {
        char first_letter = all[0]->alignment_at(pos);
        char second_letter = 0;
        bool gap = false;
        bool more_than_3 = false;
        BOOST_FOREACH (Fragment* f, all) {
            char c = f->alignment_at(pos);
            if (c == '-') {
                gap = true;
                break;
            } else if (c != first_letter) {
                if (!second_letter) {
                    second_letter = c;
                } else if (second_letter != c) {
                    more_than_3 = true;
                    break;
                }
            }
        }
        if (!gap && !more_than_3) {
            status[pos] = second_letter ? 2 : 1;
        }
    }
}

void findDiag(Ints& diag_pos, const Ints& status) {
    int length = status.size();
    for (int pos = 1; pos < length - 1; pos++) {
        int prev = status[pos - 1];
        int curr = status[pos];
        int next = status[pos + 1];
        if (prev == 1 && curr == 2 && next == 1) {
            diag_pos.push_back(pos);
        }
    }
}

void findDiagnostic(Ints& result, const Fragments& all) {
    int length = all[0]->alignment_length();
    Ints status((length)); // 1 - ident, 2 - 2 variants, or 0
    buildStatus(status, all);
    findDiag(result, status);
}

typedef std::map<std::string, int> S2I;

struct CladeCmp {
    const S2I& clade_score_;

    CladeCmp(const S2I& clade_score):
        clade_score_(clade_score) {
    }

    bool operator()(const std::string& a,
                    const std::string& b) const {
        S2I::const_iterator a_it = clade_score_.find(a);
        ASSERT_TRUE(a_it != clade_score_.end());
        S2I::const_iterator b_it = clade_score_.find(b);
        ASSERT_TRUE(b_it != clade_score_.end());
        return a_it->second > b_it->second;
    }
};

void findClades(Strings& clades, const Fragments& all) {
    Ints diag_pos;
    findDiagnostic(diag_pos, all);
    S2I clade_score;
    BOOST_FOREACH (int pos, diag_pos) {
        std::string clade;
        fragmentsToClade(clade, pos, all);
        clade_score[clade] += 1;
    }
    BOOST_FOREACH (const S2I::value_type& kv, clade_score) {
        clades.push_back(kv.first);
    }
    std::sort(clades.begin(), clades.end(),
              CladeCmp(clade_score));
}

void cladeToFragments(Fragments& group,
                      const std::string& clade,
                      const Fragments& all, char status) {
    int n = clade.length();
    ASSERT_EQ(all.size(), n);
    for (int i = 0; i < n; i++) {
        if (clade[i] == status) {
            group.push_back(all[i]);
        }
    }
}

bool groupIsGood(const Fragments& group) {
    return group.size() >= 2 && !has_repeats(group);
}

Block* makeBlock(const Fragments& group) {
    Block* new_block = new Block;
    new_block->set_weak(true);
    BOOST_FOREACH (Fragment* f, group) {
        new_block->insert(f);
    }
    return new_block;
}

Block* findGoodClade(const Strings& clades,
                     const Fragments& all) {
    BOOST_FOREACH (const std::string& clade, clades) {
        for (int i = 0; i <= 1; i++) {
            // i = 0, 1 for two sets produced by a partition
            Fragments group;
            cladeToFragments(group, clade, all, '0' + i);
            if (groupIsGood(group)) {
                return makeBlock(group);
            }
        }
    }
    return 0;
}

void excludeClade(Fragments& all, const Block* new_block) {
    Fragments all2;
    BOOST_FOREACH (Fragment* f, all) {
        if (!new_block->has(f)) {
            all2.push_back(f);
        }
    }
    all.swap(all2);
}

void SplitRepeats::process_block_impl(Block* block,
                                      ThreadData* data) const {
    SplitRepeatsData* d;
    d = boost::polymorphic_downcast<SplitRepeatsData*>(data);
    Blocks& new_blocks = d->blocks_to_insert;
    Fragments all(block->begin(), block->end());
    if (!has_repeats(all)) {
        return;
    }
    int n = 0;
    while (all.size() > 3) {
        Strings clades;
        findClades(clades, all);
        Block* new_block = findGoodClade(clades, all);
        if (new_block) {
            n += 1;
            new_block->set_name(block->name() + "g" + TO_S(n));
            new_blocks.push_back(new_block);
            excludeClade(all, new_block);
            if (groupIsGood(all)) {
                Block* new_block2 = makeBlock(all);
                n += 1;
                new_block2->set_name(block->name() +
                        "g" + TO_S(n));
                new_blocks.push_back(new_block2);
                excludeClade(all, new_block2);
            }
        } else {
            break;
        }
    }
}

void SplitRepeats::after_thread_impl(ThreadData* data) const {
    BlockSet& target = *block_set();
    SplitRepeatsData* d;
    d = boost::polymorphic_downcast<SplitRepeatsData*>(data);
    Blocks& new_blocks = d->blocks_to_insert;
    BOOST_FOREACH (Block* new_block, new_blocks) {
        target.insert(new_block);
    }
}

const char* SplitRepeats::name_impl() const {
    return "Find splittable blocks-repeats and split them";
}

}
