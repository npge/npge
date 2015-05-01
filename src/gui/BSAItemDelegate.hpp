/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BSAITEMDELEGATE_HPP
#define BSAITEMDELEGATE_HPP

#include <QtGui>

// See http://stackoverflow.com/a/24805735

class BSAItemDelegate : public QStyledItemDelegate {
    Q_OBJECT
public:
    explicit BSAItemDelegate(QTableView* view);

    void initStyleOption(QStyleOptionViewItem* option,
                         const QModelIndex& index) const;

private:
    QTableView* view_;
};

#endif // BSAITEMDELEGATE_HPP
