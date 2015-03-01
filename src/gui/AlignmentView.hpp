/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef ALIGNMENTVIEW_HPP
#define ALIGNMENTVIEW_HPP

#include <QTableView>

#ifndef Q_MOC_RUN
#include "global.hpp"
#include "gui-global.hpp"
#endif

using namespace npge;

class AlignmentView : public QTableView {
    Q_OBJECT
public:
    explicit AlignmentView(QWidget* parent = 0);

    void set_model(AlignmentModel* model);

    typedef std::map<Sequence*, Fragment*> Seq2Fragment;

    void keyPressEvent(QKeyEvent* event);

    Blocks make_selected_blocks() const;

    Block* make_selected_block(
        const QItemSelectionRange& range) const;

public slots:
    void select_fragment(Fragment* fragment);

signals:
    void fragment_selected(Fragment* fragment, int column);
    void jump_to(Fragment* to, int column);

private slots:
    void clicked_f(const QModelIndex& index);
};

#endif // ALIGNMENTVIEW_HPP

