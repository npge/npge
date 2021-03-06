/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef READING_THREAD_HPP
#define READING_THREAD_HPP

#include <QtCore>

#ifndef Q_MOC_RUN
#include "global.hpp"
#include "GuiBSs.hpp"
#include "mainwindow.hpp"
#endif

class ReadingThread : public QThread {
    Q_OBJECT

public:
    ReadingThread(Meta* meta, GuiBSs* bss,
                  std::string fname, MainWindow* parent = 0);

    void run();

signals:
    void readingFinished(QString message);

private:
    Meta* meta_;
    GuiBSs* bss_;
    std::string fname_;

    void run_impl();
};

#endif

