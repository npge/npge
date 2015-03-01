/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <sstream>
#include <boost/scoped_ptr.hpp>
#include <boost/bind.hpp>
#include <QtGui>

#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "Block.hpp"
#include "move_rows.hpp"
#include "block_hash.hpp"
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
        bool ls = model()->headerData(logicalIndex, Qt::Horizontal,
                                      Qt::UserRole).value<bool>();
        QString text = model()->headerData(logicalIndex, Qt::Horizontal,
                                           Qt::DisplayRole).value<QString>();
        if (ls) {
            QRect rect_top(rect);
            rect_top.setBottom((rect.top() + rect.bottom()) / 2);
            painter->fillRect(rect_top, Qt::red);
        }
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
    if (!e) {
        return;
    }
    bool ctrl = e->modifiers().testFlag(Qt::ControlModifier);
    bool shift = e->modifiers().testFlag(Qt::ShiftModifier);
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
        ASSERT_TRUE(m);
        move_view_rows(this, e->key() == Qt::Key_Up,
                       boost::bind(&AlignmentModel::move_rows,
                                   m, _1, _2));
    } else if (ctrl && (left || right)) {
        QModelIndex index = currentIndex();
        int row = index.row();
        int col = index.column();
        AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
        ASSERT_TRUE(m);
        GeneInfo go;
        Fragment* current_gene = m->test_genes(index, &go);
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
    } else if (shift && (left || right)) {
        QModelIndex index = currentIndex();
        int row = index.row();
        int col = index.column();
        AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
        ASSERT_TRUE(m);
        bool low_col = m->is_low_col(col);
        while (true) {
            if (left) {
                col -= 1;
            } else if (right) {
                col += 1;
            }
            if (col <= -1 || col >= m->columnCount()) {
                return;
            }
            bool low_col1 = m->is_low_col(col);
            if (low_col1 != low_col) {
                // region start/end
                break;
            }
        }
        index = m->index(row, col);
        selectionModel()->clearSelection();
        setCurrentIndex(index);
        scrollTo(index);
    } else if (r_jump || l_jump) {
        AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
        Fragment* f = m->fragment_at(r);
        int ori = r_jump ? 1 : -1;
        Fragment* neighbour = m->logical_neighbor(f, ori);
        Sequence* seq = f->seq();
        if (neighbour) {
            int col = (f->ori() * neighbour->ori() * ori == 1) ? 0 :
                      (neighbour->block()->alignment_length() - 1);
            emit jump_to(neighbour, col);
        }
    } else if (ctrl && c) {
        Blocks selected_blocks = make_selected_blocks();
        QString text;
        if (selected_blocks.size() == 1 &&
                selected_blocks.front()->size() == 1) {
            Fragment* f = selected_blocks.front()->front();
            text = QString::fromStdString(f->str());
        } else {
            std::stringstream ss;
            foreach (Block* block, selected_blocks) {
                ss << *block;
                delete block;
            }
            text = QString::fromStdString(ss.str());
        }
        QApplication::clipboard()->setText(text);
    } else {
        QTableView::keyPressEvent(e);
    }
}

Blocks AlignmentView::make_selected_blocks() const {
    AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
    int genomes = genomes_number(*(m->block_set()));
    Blocks result;
    QItemSelectionModel* sm = selectionModel();
    foreach (const QItemSelectionRange& range,
            sm->selection()) {
        Block* slice = make_selected_block(range);
        set_canonical_name(slice, genomes);
        result.push_back(slice);
    }
    return result;
}

Block* AlignmentView::make_selected_block(
    const QItemSelectionRange& range) const {
    AlignmentModel* m = dynamic_cast<AlignmentModel*>(model());
    boost::scoped_ptr<Block> rows_copy(new Block);
    rows_copy->set_weak(true);
    int top = range.top();
    int bottom = range.bottom();
    for (int i = top; i <= bottom; i++) {
        rows_copy->insert(m->fragment_at(i));
    }
    int left = range.left();
    int right = range.right();
    bool alignment = true;
    return rows_copy->slice(left, right, alignment);
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
    Fragment* fragment = m->fragment_at(index.row());
    Fragment* f = const_cast<Fragment*>(fragment); // FIXME
    emit fragment_selected(f, index.column());
}

