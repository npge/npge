/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui>

#include "Meta.hpp"
#include "GuiBSs.hpp"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(int argc, char** argv,
                        QWidget* parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow* ui;
    QLabel* loading_;
    npge::Meta meta_;
    GuiBSs bss_;

private slots:
    void onReadingFinished(QString message);
};

#endif // MAINWINDOW_H

