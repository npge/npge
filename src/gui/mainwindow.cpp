#include <fstream>

#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSetWidget.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "config.hpp"

using namespace bloomrepeats;

BlockSetPtr test_bs;

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    test_bs = new_bs();
    std::ifstream test_file("test.fasta");
    test_file >> *test_bs;
    ui->verticalLayout_2->addWidget(new BlockSetWidget(test_bs));
}

MainWindow::~MainWindow() {
    delete ui;
}

