#ifndef ALIGNMENTVIEW_HPP
#define ALIGNMENTVIEW_HPP

#include <QTableView>

class AlignmentView : public QTableView {
    Q_OBJECT
public:
    explicit AlignmentView(QWidget* parent = 0);

    void keyPressEvent(QKeyEvent* event);

signals:

public slots:

};

#endif // ALIGNMENTVIEW_HPP

