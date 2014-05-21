#include <fstream>

#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSetWidget.hpp"
#include "BlockSet.hpp"
#include "bsa_algo.hpp"

using namespace bloomrepeats;

BlockSetPtr pangenome_bs;
BlockSetPtr genes_bs;
BlockSetPtr split_parts;
BlockSetPtr low_similarity;

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    showMaximized();
    //
    pangenome_bs = new_bs();
    std::ifstream pangenome_file("pangenome-merged.fasta");
    pangenome_file >> *pangenome_bs;
    //
    genes_bs = new_bs();
    genes_bs->add_sequences(pangenome_bs->seqs());
    std::ifstream genes_file("genes.fasta");
    genes_file >> *genes_bs;
    //
    split_parts = new_bs();
    split_parts->add_sequences(pangenome_bs->seqs());
    std::ifstream split_file("pangenome-merged-split.fasta");
    split_file >> *split_parts;
    //
    low_similarity = new_bs();
    low_similarity->add_sequences(pangenome_bs->seqs());
    std::ifstream low_file("pangenome-merged-low.fasta");
    low_file >> *low_similarity;
    //
    std::ifstream test_bsaln("pangenome-merged.bsa");
    bsa_input(*pangenome_bs, test_bsaln);
    //
    BlockSetWidget* bsw = new BlockSetWidget(pangenome_bs);
    bsw->set_genes(genes_bs);
    bsw->set_split_parts(split_parts);
    bsw->set_low_similarity(low_similarity);
    ui->verticalLayout_2->addWidget(bsw);
}

MainWindow::~MainWindow() {
    delete ui;
}

