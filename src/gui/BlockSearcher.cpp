/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Meta.hpp"
#include "BlockSearcher.hpp"

using namespace npge;

void BlockSearcher::run() {
    try {
        MetaThreadKeeper mtk(meta_);
        hits_constructor_();
        filtered_blocks_->clear();
        BOOST_FOREACH (const Block* block, *blocks_) {
            if (block_checker_(block)) {
                filtered_blocks_->push_back(block);
            }
        }
        filtered_blocks_->sort();
        emit searchingFinished("");
    } catch (std::exception& e) {
        emit searchingFinished(e.what());
    } catch (...) {
        emit searchingFinished("Unknown error");
    }
}

