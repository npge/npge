/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "SliceNless.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "cast.hpp"

namespace npge {

SliceNless::SliceNless() {
    add_gopt("max-ns", "Maximum number of subsequent N's "
             "in consensus", "MAX_NS");
}

struct SNData : public ThreadData {
    Blocks to_delete_;
    Blocks to_insert_;
    std::string pattern_;
    std::string repl_;
    int max_ns_;

    SNData(const Processor* p) {
        max_ns_ = p->opt_value("max-ns").as<int>();
        pattern_ = std::string(max_ns_, 'N');
        repl_ = std::string(max_ns_, '*');
    }
};

ThreadData* SliceNless::before_thread_impl() const {
    return new SNData(this);
}

void SliceNless::process_block_impl(Block* block,
                                    ThreadData* d) const {
    SNData* data = D_CAST<SNData*>(d);
    std::string ccc = block->consensus_string();
    int l = ccc.size();
    if (ccc.find(data->pattern_) == std::string::npos) {
        return;
    }
    data->to_delete_.push_back(block);
    // replace NNN with ***
    using namespace boost::algorithm;
    replace_all(ccc, data->pattern_, data->repl_);
    // replace **NN with ****
    for (int i = 1; i < l; i++) {
        if (ccc[i - 1] == '*' && ccc[i] == 'N') {
            ccc[i] = '*';
        }
    }
    // find uncensored parts
    int region_start = -1;
    for (int i = 0; i < l; i++) {
        if (ccc[i] == '*') {
            if (region_start != -1) {
                Block* bn = block->slice(region_start, i - 1);
                data->to_insert_.push_back(bn);
            }
            region_start = -1;
        } else if (region_start == -1) {
            region_start = i;
        }
    }
    if (region_start != -1) {
        Block* bn = block->slice(region_start, l - 1);
        data->to_insert_.push_back(bn);
    }
}

void SliceNless::after_thread_impl(ThreadData* d) const {
    SNData* data = D_CAST<SNData*>(d);
    BlockSet& bs = *block_set();
    BOOST_FOREACH (Block* b, data->to_delete_) {
        bs.erase(b);
    }
    BOOST_FOREACH (Block* b, data->to_insert_) {
        bs.insert(b);
    }
}

const char* SliceNless::name_impl() const {
    return "Find blocks with long N's and slice them to "
           "good parts";
}

}

