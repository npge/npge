/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>

#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSetWidget.hpp"
#include "BlockSet.hpp"
#include "bsa_algo.hpp"
#include "name_to_stream.hpp"
#include "In.hpp"

using namespace npge;

BlockSetPtr pangenome_bs;
BlockSetPtr genes_bs;
BlockSetPtr split_parts;
BlockSetPtr low_similarity;

typedef boost::shared_ptr<std::istream> IPtr;

static void read_bs(BlockSetPtr bs, std::string name) {
    In p_in;
    p_in.set_block_set(bs);
    p_in.set_opt_value("in-blocks", name);
    p_in.run();
}

MainWindow::MainWindow(int argc, char** argv,
                       QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    showMaximized();
    //
    pangenome_bs = new_bs();
    if (argc >= 2) {
        read_bs(pangenome_bs, argv[1]);
    } else {
        read_bs(pangenome_bs, "pangenome.bs");
        //
        genes_bs = new_bs();
        genes_bs->add_sequences(pangenome_bs->seqs());
        read_bs(genes_bs, "features.bs");
        //
        split_parts = new_bs();
        split_parts->add_sequences(pangenome_bs->seqs());
        read_bs(split_parts, "split.bs");
        //
        low_similarity = new_bs();
        low_similarity->add_sequences(pangenome_bs->seqs());
        read_bs(low_similarity, "low.bs");
        //
        IPtr test_bsaln = name_to_istream("pangenome.bsa");
        bsa_input(*pangenome_bs, *test_bsaln);
    }
    BlockSetWidget* bsw = new BlockSetWidget(pangenome_bs);
    bsw->set_genes(genes_bs);
    bsw->set_split_parts(split_parts);
    bsw->set_low_similarity(low_similarity);
    ui->verticalLayout_2->addWidget(bsw);
}

MainWindow::~MainWindow() {
    delete ui;
}

