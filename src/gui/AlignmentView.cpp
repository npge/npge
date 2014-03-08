#include <QtGui>

#include "AlignmentView.hpp"

class HorizontalHeader : public QHeaderView {
public:
    HorizontalHeader(QWidget* parent):
        QHeaderView(Qt::Horizontal, parent) {
        setDefaultSectionSize(17);
        setStyleSheet("::section { border : 0px; }");
        setResizeMode(QHeaderView::Fixed);
    }

    void paintSection(QPainter* painter, const QRect& rect,
                      int logicalIndex) const {
        QColor bg = model()->headerData(logicalIndex, Qt::Horizontal,
                                        Qt::BackgroundRole).value<QColor>();
        QColor fg = model()->headerData(logicalIndex, Qt::Horizontal,
                                        Qt::ForegroundRole).value<QColor>();
        QString text = model()->headerData(logicalIndex, Qt::Horizontal,
                                           Qt::DisplayRole).value<QString>();
        QRect rect_bottom(rect);
        rect_bottom.setTop((rect.top() + rect.bottom()) / 2);
        painter->fillRect(rect_bottom, bg);
        QString t1 = text.left(1) + "\n", t2 = "\n" + text.right(1);
        painter->drawText(rect, Qt::AlignCenter, t1);
        painter->setPen(fg);
        painter->drawText(rect, Qt::AlignCenter, t2);
    }
};

AlignmentView::AlignmentView(QWidget* parent) :
    QTableView(parent) {
    setShowGrid(false);
    setHorizontalHeader(new HorizontalHeader(this));
    verticalHeader()->setResizeMode(QHeaderView::Fixed);
    verticalHeader()->setDefaultSectionSize(20);
}

