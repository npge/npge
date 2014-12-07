/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef MOVE_ROWS_HPP
#define MOVE_ROWS_HPP

#ifndef Q_MOC_RUN
#include <boost/function.hpp>
#endif
#include <QtGui>

typedef std::vector<int> Rows;
typedef boost::function<void(Rows&, bool)> ModelAction;

/** Move rows of view up or down, preserving selection */
void move_view_rows(QTableView* view, bool up,
                    ModelAction model_action);

#endif

