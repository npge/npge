#include <fstream>

#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSetWidget.hpp"
#include "BlockSet.hpp"
#include "block_set_alignment.hpp"
#include "Block.hpp"
#include "config.hpp"
#include "AddGenes.hpp"

using namespace bloomrepeats;

BlockSetPtr test_bs;
BlockSetPtr genes_bs;
BlockSetPtr split_parts;
BlockSetPtr low_similarity;

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    showMaximized();
    //
    test_bs = new_bs();
    std::ifstream test_file("test.fasta");
    test_file >> *test_bs;
    //
    AddGenes ag;
    genes_bs = ag.block_set();
    genes_bs->add_sequences(test_bs->seqs());
    ag.set_opt_value("in-genes", std::string("test.genes"));
    ag.run();
    //
    split_parts = new_bs();
    split_parts->add_sequences(test_bs->seqs());
    std::ifstream split_file("split.fasta");
    split_file >> *split_parts;
    //
    low_similarity = new_bs();
    low_similarity->add_sequences(test_bs->seqs());
    std::ifstream low_file("low.fasta");
    low_file >> *low_similarity;
    //
    std::ifstream test_bsaln("test.bsaln");
    bsa_input(*test_bs, test_bsaln);
    //
    BlockSetWidget* bsw = new BlockSetWidget(test_bs);
    bsw->set_genes(genes_bs);
    bsw->set_split_parts(split_parts);
    bsw->set_low_similarity(low_similarity);
    ui->verticalLayout_2->addWidget(bsw);
}

MainWindow::~MainWindow() {
    delete ui;
}

