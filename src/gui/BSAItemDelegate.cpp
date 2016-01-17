/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "BSAItemDelegate.hpp"

// See http://stackoverflow.com/a/24805735

BSAItemDelegate::BSAItemDelegate(QTableView* view):
    QStyledItemDelegate(view), view_(view) {
}

void BSAItemDelegate::initStyleOption(
    QStyleOptionViewItem* option,
    const QModelIndex& index) const {
    QStyledItemDelegate::initStyleOption(option, index);
    QStyleOptionViewItemV4 *v4 =
        qstyleoption_cast<QStyleOptionViewItemV4 *>(option);
    QItemSelectionModel* m = view_->selectionModel();
    if (m->isSelected(index)) {
        //v4->backgroundBrush = QBrush(Qt::yellow);
        v4->font.setBold(true);
    }
}
