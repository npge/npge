/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <exception>
#include <QtGui>

#include "mainwindow.hpp"
#include "htmlencode.hpp"

class MyApplication : public QApplication {
public:
    // Pass argc by reference
    // Otherwise we give reference to stack variable
    // Which results in memory error (segfault)
    // When breaking Qt's  main loop
    // (e.g. exec() of dialogs or message boxes)
    MyApplication(int& argc, char* argv[]):
        QApplication(argc, argv) {
    }

    bool notify(QObject* receiver, QEvent* event) {
        try {
            return QApplication::notify(receiver, event);
        } catch (const std::exception& error) {
            using namespace npge;
            QString what = QString::fromStdString(htmlencode(error.what()));
            QString message = "<b>The error occurred</b>."
                              "<br/><br/>"
                              "Description for developers:"
                              "<br/><br/>" + what;
            QErrorMessage::qtHandler()->resize(400, 300);
            QErrorMessage::qtHandler()->showMessage(message);
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

