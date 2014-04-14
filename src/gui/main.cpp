#include <exception>
#include <QtGui>
#include "mainwindow.hpp"

class MyApplication : public QApplication {
public:
    MyApplication(int argc, char* argv[]):
        QApplication(argc, argv)
    { }

    bool notify(QObject* receiver, QEvent* e) {
        try {
            QApplication::notify(receiver, e);
        } catch (const std::exception& e) {
            QString what = QString::fromStdString(e.what());
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
    MainWindow w;
    w.show();
    return a.exec();
}

