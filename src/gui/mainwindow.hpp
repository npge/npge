/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui>

#ifndef Q_MOC_RUN
#include "Meta.hpp"
#include "GuiBSs.hpp"
#endif

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

signals:
    void exceptionThrown(QString message);

private slots:
    void onReadingFinished(QString message);

    void onExceptionThrown(QString message);
};

#endif // MAINWINDOW_H

