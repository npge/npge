/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Meta.hpp"
#include "BlockSearcher.hpp"
#include "BlockSetModel.hpp"

using namespace npge;

void BlockSearcher::run() {
    try {
        MetaThreadKeeper mtk(meta_);
        model_->construct_hits();
        filtered_blocks_->clear();
        BOOST_FOREACH (const Block* block, *blocks_) {
            if (model_->check_block(block)) {
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

