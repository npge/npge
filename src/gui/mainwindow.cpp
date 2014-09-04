/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSetWidget.hpp"
#include "ReadingThread.hpp"
#include "Exception.hpp"

using namespace npge;

MainWindow::MainWindow(int argc, char** argv,
                       QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    showMaximized();
    //
    std::string fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    bsw_ = new BlockSetWidget(new_bs());
    ui->verticalLayout_2->addWidget(bsw_);
    ReadingThread* thread = new ReadingThread(&meta_, &bss_,
            fname, this);
    connect(thread, SIGNAL(readingFinished(QString)),
            this, SLOT(onReadingFinished(QString)),
            Qt::QueuedConnection);
    thread->start();
}

void MainWindow::onReadingFinished(QString message) {
    if (message.isEmpty()) {
        bsw_->set_block_set(bss_.pangenome_bs_);
        bsw_->set_genes(bss_.genes_bs_);
        bsw_->set_split_parts(bss_.split_parts_);
        bsw_->set_low_similarity(bss_.low_similarity_);
    } else {
        throw Exception(message.toStdString());
    }
}

MainWindow::~MainWindow() {
    delete ui;
}

