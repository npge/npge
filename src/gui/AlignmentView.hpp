#ifndef ALIGNMENTVIEW_HPP
#define ALIGNMENTVIEW_HPP

#include <QTableView>

#include "global.hpp"

using namespace bloomrepeats;

class AlignmentView : public QTableView {
    Q_OBJECT
public:
    explicit AlignmentView(QWidget* parent = 0);

    void keyPressEvent(QKeyEvent* event);

signals:
    void jump_to(Fragment* to, int column);

public slots:

};

#endif // ALIGNMENTVIEW_HPP

