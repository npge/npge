/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Connector.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

namespace npge {

Connector::Connector() {
    declare_bs("target", "Target blockset");
    add_opt("connect-circular", "Connect last and first fragments "
            "of circular chromosomes", false);
}

static struct FragmentCompare {
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        return *f1 < *f2;
    }
} fragment_compare;

void Connector::run_impl() const {
    bool circular = opt_value("connect-circular").as<bool>();
    typedef std::vector<Fragment*> Fs;
    typedef std::map<Sequence*, Fs> Seq2Fs;
    Seq2Fs seq2fs;
    BOOST_FOREACH (Block* block, *block_set()) {
        BOOST_FOREACH (Fragment* fragment, *block) {
            seq2fs[fragment->seq()].push_back(fragment);
            fragment->disconnect();
        }
    }
    BOOST_FOREACH (Seq2Fs::value_type& seq_and_fs, seq2fs) {
        Sequence* seq = seq_and_fs.first;
        Fs& fs = seq_and_fs.second;
        std::sort(fs.begin(), fs.end(), fragment_compare);
        for (int i = 1; i < fs.size(); i++) {
            Fragment::connect(fs[i - 1], fs[i]);
        }
        if (circular && seq->circular()) {
            Fragment::connect(fs.back(), fs.front());
        }
    }
}

const char* Connector::name_impl() const {
    return "Connect all the fragments (prev-next)";
}

}

