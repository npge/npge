#ifndef MOVE_ROWS_HPP
#define MOVE_ROWS_HPP

#include <boost/function.hpp>
#include <QtGui>

typedef std::vector<int> Rows;
typedef boost::function<void(Rows&, bool)> ModelAction;

/** Move rows of view up or down, preserving selection */
void move_view_rows(QTableView* view, bool up,
                    ModelAction model_action);

#endif

