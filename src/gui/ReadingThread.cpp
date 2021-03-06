/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "ReadingThread.hpp"
#include "Meta.hpp"
#include "BlockSet.hpp"
#include "bsa_algo.hpp"
#include "name_to_stream.hpp"
#include "Read.hpp"
#include "cast.hpp"

using namespace npge;

typedef boost::shared_ptr<std::istream> IPtr;

static void read_bs(BlockSetPtr bs, std::string name,
                    bool required) {
    if (!required) {
        try {
            // test file
            name_to_istream(name);
        } catch (...) {
            return;
        }
    }
    Read p_in;
    p_in.set_block_set(bs);
    p_in.set_opt_value("in-blocks", name);
    p_in.run();
}

ReadingThread::ReadingThread(Meta* meta, GuiBSs* bss,
                             std::string fname,
                             MainWindow* parent):
    QThread(parent), meta_(meta), bss_(bss),
    fname_(fname) {
}

void ReadingThread::run() {
    try {
        run_impl();
        emit readingFinished("");
    } catch (std::exception& e) {
        emit readingFinished(e.what());
    } catch (...) {
        emit readingFinished("Unknown error");
    }
}

static void read_ba(BlockSetPtr bs, std::string name) {
    IPtr test_bsaln;
    try {
        test_bsaln = name_to_istream(name);
    } catch (...) {
    }
    if (test_bsaln) {
        bsa_input(*bs, *test_bsaln);
    }
}

void ReadingThread::run_impl() {
    MetaThreadKeeper mtk(meta_);
    BlockSetPtr& pangenome_bs = bss_->pangenome_bs_;
    BlockSetPtr& genes_bs = bss_->genes_bs_;
    BlockSetPtr& split_parts = bss_->split_parts_;
    BlockSetPtr& low_similarity = bss_->low_similarity_;
    BlockSetPtr& global_blocks = bss_->global_blocks_;
    pangenome_bs = new_bs();
    if (!fname_.empty()) {
        read_bs(pangenome_bs, fname_, true);
    } else {
        read_bs(pangenome_bs, "pangenome/pangenome.bs", true);
        //
        genes_bs = new_bs();
        genes_bs->add_sequences(pangenome_bs->seqs());
        read_bs(genes_bs, "genes/features.bs", false);
        //
        split_parts = new_bs();
        split_parts->add_sequences(pangenome_bs->seqs());
        read_bs(split_parts, "extra-blocks/split.bs", false);
        //
        low_similarity = new_bs();
        low_similarity->add_sequences(pangenome_bs->seqs());
        read_bs(low_similarity, "extra-blocks/low.bs", false);
        //
        read_ba(pangenome_bs, "pangenome/pangenome.ba");
        //
        global_blocks = new_bs();
        global_blocks->add_sequences(pangenome_bs->seqs());
        read_bs(global_blocks, "global-blocks/blocks.bs",
                false);
        //
        read_ba(global_blocks, "global-blocks/blocks.ba");
    }
}

