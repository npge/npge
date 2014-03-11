#include <set>
#include <QtGui>

#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "throw_assert.hpp"

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

void AlignmentView::keyPressEvent(QKeyEvent* e) {
    bool ctrl = e->modifiers().testFlag(Qt::ControlModifier);
    bool up_down = e->key() == Qt::Key_Up || e->key() == Qt::Key_Down;
    QItemSelectionModel* sm = selectionModel();
    int r = sm->currentIndex().row();
    bool right = e->key() == Qt::Key_Right;
    bool left = e->key() == Qt::Key_Left;
    bool r_jump = right && sm->currentIndex().column() ==
                  model()->columnCount() - 1;
    bool l_jump = left && sm->currentIndex().column() == 0;
    if (ctrl && up_down) {
        std::set<int> rows_set;
        foreach (QModelIndex index, selectedIndexes()) {
            rows_set.insert(index.row());
        }
        std::vector<int> rows(rows_set.begin(), rows_set.end());
        AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
        BOOST_ASSERT(m);
        m->move_rows(rows, e->key() == Qt::Key_Up);
        QItemSelectionModel* sm = selectionModel();
        sm->clear();
        foreach (int row, rows) {
            sm->select(m->index(row, 0), QItemSelectionModel::Select
                       | QItemSelectionModel::Rows);
        }
    } else if (ctrl && (left || right)) {
        QModelIndex index = currentIndex();
        int row = index.row();
        int col = index.column();
        AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
        BOOST_ASSERT(m);
        bool _;
        const Fragment* current_gene = m->test_genes(index, _, _, _);
        while (true) {
            if (left) {
                col -= 1;
            } else if (right) {
                col += 1;
            }
            if (col <= -1 || col >= m->columnCount()) {
                return;
            }
            index = m->index(row, col);
            bool gap = m->test_gap(index);
            bool gene = m->test_genes(index, _, _, _) != current_gene;
            if (!gap && gene) {
                // gene changed
                break;
            }
        }
        selectionModel()->clearSelection();
        setCurrentIndex(index);
        scrollTo(index);
    } else if (r_jump || l_jump) {
        AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
        const Fragment* f = m->fragment_at(r);
        int ori = r_jump ? 1 : -1;
        Fragment* neighbour = f->logical_neighbor(ori);
        if (neighbour) {
            int col = (f->ori() * neighbour->ori() * ori == 1) ? 0 :
                      (neighbour->block()->alignment_length() - 1);
            qDebug() << col;
            emit jump_to(neighbour, col);
        }
    } else {
        QTableView::keyPressEvent(e);
    }
}

