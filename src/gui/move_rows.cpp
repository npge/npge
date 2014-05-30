/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>

#include "move_rows.hpp"

struct Range {
    int row, left, right;

    Range(int rw, int lft, int rght):
        row(rw), left(lft), right(rght) {
    }
};

void move_view_rows(QTableView* view, bool up,
                    ModelAction model_action) {
    QAbstractItemModel* m = view->model();
    QItemSelectionModel* sm = view->selectionModel();
    std::vector<Range> old_selection;
    foreach (const QItemSelectionRange& range, sm->selection()) {
        if (range.height() == 1) {
            Range r1(range.top(), range.left(), range.right());
            old_selection.push_back(r1);
        } else {
            for (int r = range.top(); r <= range.bottom(); r++) {
                Range r1(r, range.left(), range.right());
                old_selection.push_back(r1);
            }
        }
    }
    std::set<int> rows_set;
    foreach (const Range& range, old_selection) {
        rows_set.insert(range.row);
    }
    std::vector<int> rows(rows_set.begin(), rows_set.end());
    std::vector<int> rows_orig = rows;
    QModelIndexList selected  = sm->selectedIndexes();
    QModelIndex cur = view->currentIndex();
    model_action(rows, up);
    std::map<int, int> old2new;
    for (int i = 0; i < rows.size(); i++) {
        old2new[rows_orig[i]] = rows[i];
    }
    sm->clear();
    QModelIndex x = m->index(old2new[cur.row()], cur.column());
    view->setCurrentIndex(x);
    QItemSelection new_selection;
    foreach (const Range& range, old_selection) {
        int old_row = range.row;
        int new_row = old2new[old_row];
        QModelIndex tl = m->index(new_row, range.left);
        QModelIndex br = m->index(new_row, range.right);
        new_selection << QItemSelectionRange(tl, br);
    }
    sm->select(new_selection, QItemSelectionModel::Select);
}

