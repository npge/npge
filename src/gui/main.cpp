/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <exception>
#include <QtGui>

#include "mainwindow.hpp"
#include "htmlencode.hpp"

class MyApplication : public QApplication {
public:
    MyApplication(int argc, char* argv[]):
        QApplication(argc, argv) {
    }

    bool notify(QObject* receiver, QEvent* e) {
        try {
            QApplication::notify(receiver, e);
        } catch (const std::exception& e) {
            using namespace npge;
            QString what = QString::fromStdString(htmlencode(e.what()));
            QString error = "<b>The error occured</b>.<br><br>"
                            "Description for developers:<br><br>" + what;
            QErrorMessage::qtHandler()->resize(400, 300);
            QErrorMessage::qtHandler()->showMessage(error);
            return false;
        }
    }
};

int main(int argc, char* argv[]) {
    MyApplication a(argc, argv);
    MainWindow w(argc, argv);
    w.show();
    return a.exec();
}

