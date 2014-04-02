#include <set>
#include <boost/bind.hpp>
#include <QtGui>

#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "Block.hpp"
#include "move_rows.hpp"
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
    bool c = e->key() == Qt::Key_C;
    if (ctrl && up_down) {
        AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
        BOOST_ASSERT(m);
        move_view_rows(this, e->key() == Qt::Key_Up,
                       boost::bind(&AlignmentModel::move_rows,
                                   m, _1, _2));
    } else if (ctrl && (left || right)) {
        QModelIndex index = currentIndex();
        int row = index.row();
        int col = index.column();
        AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
        BOOST_ASSERT(m);
        GeneInfo go;
        const Fragment* current_gene = m->test_genes(index, &go);
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
            bool gene = m->test_genes(index, &go) != current_gene;
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
        Sequence* seq = f->seq();
        bool circular = seq->circular();
        if (!neighbour && circular && f == seq2first_[seq]) {
            neighbour = seq2last_[seq];
        }
        if (!neighbour && circular && f == seq2last_[seq]) {
            neighbour = seq2first_[seq];
        }
        if (neighbour) {
            int col = (f->ori() * neighbour->ori() * ori == 1) ? 0 :
                      (neighbour->block()->alignment_length() - 1);
            emit jump_to(neighbour, col);
        }
    } else if (ctrl && c) {
        QModelIndexList selected  = sm->selectedIndexes();
        qSort(selected);
        QModelIndex prev;
        QString text;
        foreach (QModelIndex index, selected) {
            if (prev.isValid() && index.row() != prev.row()) {
                text += "\n";
            } else if (prev.isValid() &&
                       index.column() != prev.column() + 1) {
                text += " ";
            }
            text += model()->data(index).toString();
            prev = index;
        }
        QApplication::clipboard()->setText(text);
    } else {
        QTableView::keyPressEvent(e);
    }
}

void AlignmentView::set_model(AlignmentModel* new_model) {
    QAbstractItemModel* m = model();
    QItemSelectionModel* sm = selectionModel();
    setModel(new_model);
    if (m) {
        m->deleteLater();
    }
    if (sm) {
        sm->deleteLater();
    }
    connect(selectionModel(),
            SIGNAL(currentChanged(QModelIndex, QModelIndex)),
            this, SLOT(clicked_f(QModelIndex)));
}

void AlignmentView::set_first_last(const Seq2Fragment& first,
                                   const Seq2Fragment& last) {
    seq2first_ = first;
    seq2last_ = last;
}

void AlignmentView::select_fragment(Fragment* fragment) {
    AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
    int row = m->fragment_index(fragment);
    selectionModel()->clearSelection();
    selectionModel()->select(m->index(row, 0),
                             QItemSelectionModel::Select |
                             QItemSelectionModel::Rows);
    int col = columnAt(0);
    scrollTo(m->index(row, col));
}

void AlignmentView::clicked_f(const QModelIndex& index) {
    AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
    const Fragment* fragment = m->fragment_at(index.row());
    Fragment* f = const_cast<Fragment*>(fragment); // FIXME
    emit fragment_selected(f, index.column());
}

