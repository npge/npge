#include <fstream>

#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"

using namespace bloomrepeats;

BlockSet test_bs;

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    std::ifstream test_file("test.fasta");
    test_file >> test_bs;
    if (!test_bs.empty()) {
        Block* block = test_bs.front();
        AlignmentModel* model = new AlignmentModel(block, this);
        AlignmentView* view = new AlignmentView;
        view->setModel(model);
        ui->verticalLayout_2->addWidget(view);
    }
}

MainWindow::~MainWindow() {
    delete ui;
}

