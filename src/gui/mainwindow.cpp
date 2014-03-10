#include <fstream>

#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSetWidget.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "config.hpp"
#include "AddGenes.hpp"

using namespace bloomrepeats;

BlockSetPtr test_bs;
BlockSetPtr genes_bs;

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    test_bs = new_bs();
    std::ifstream test_file("test.fasta");
    test_file >> *test_bs;
    AddGenes ag;
    genes_bs = ag.block_set();
    genes_bs->add_sequences(test_bs->seqs());
    std::vector<std::string> test_genes;
    test_genes.push_back("test.genes");
    ag.set_opt_value("in-genes", test_genes);
    ag.run();
    BlockSetWidget* bsw = new BlockSetWidget(test_bs);
    bsw->set_genes(genes_bs);
    ui->verticalLayout_2->addWidget(bsw);
}

MainWindow::~MainWindow() {
    delete ui;
}

