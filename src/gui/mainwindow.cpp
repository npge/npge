/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
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
    setWindowIcon(QIcon(":/npge.ico"));
    connect(this, SIGNAL(exceptionThrown(QString)),
            this, SLOT(onExceptionThrown(QString)),
            Qt::QueuedConnection);
    ui->setupUi(this);
    showMaximized();
    loading_ = new QLabel("Reading files...", this);
    loading_->setAlignment(
        Qt::AlignHCenter | Qt::AlignVCenter);
    ui->verticalLayout_2->addWidget(loading_);
    //
    std::string fname;
    if (argc >= 2) {
        fname = argv[1];
    }
    ReadingThread* thread = new ReadingThread(&meta_, &bss_,
            fname, this);
    connect(thread, SIGNAL(readingFinished(QString)),
            this, SLOT(onReadingFinished(QString)),
            Qt::QueuedConnection);
    thread->start();
}

void MainWindow::onReadingFinished(QString message) {
    ui->verticalLayout_2->removeWidget(loading_);
    if (message.isEmpty()) {
        BlockSetWidget* bsw, *gb; // gb is global blocks
        bsw = new BlockSetWidget(bss_.pangenome_bs_);
        if (bss_.global_blocks_) {
            // FIXME memory leak (gb)
            gb = new BlockSetWidget(bss_.global_blocks_);
            BlockSetWidget::moveBsaWidget(bsw, gb);
        }
        ui->verticalLayout_2->addWidget(bsw);
        bsw->set_block_set(bss_.pangenome_bs_);
        if (bss_.genes_bs_) {
            bsw->set_genes(bss_.genes_bs_);
        }
        if (bss_.split_parts_) {
            bsw->set_split_parts(bss_.split_parts_);
        }
        if (bss_.low_similarity_) {
            bsw->set_low_similarity(bss_.low_similarity_);
        }
    } else {
        emit exceptionThrown(message);
    }
}

void MainWindow::onExceptionThrown(QString message) {
    throw Exception(message.toStdString());
}

MainWindow::~MainWindow() {
    delete ui;
}

