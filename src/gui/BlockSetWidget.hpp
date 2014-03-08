#ifndef BLOCKSETWIDGET_HPP
#define BLOCKSETWIDGET_HPP

#include <QtGui>

#include "gui-global.hpp"
#include "global.hpp"

using namespace bloomrepeats;

namespace Ui {
class BlockSetWidget;
}

class BlockSetModel;

class BlockSetWidget : public QWidget {
    Q_OBJECT

public:
    explicit BlockSetWidget(BlockSetPtr block_set = BlockSetPtr(),
                            QWidget* parent = 0);
    ~BlockSetWidget();

    void set_block_set(BlockSetPtr block_set);

private:
    Ui::BlockSetWidget* ui;
    AlignmentView* alignment_view_;
    AlignmentModel* alignment_model_;
    BlockSetModel* block_set_model_;
    QSortFilterProxyModel* proxy_model_;

private slots:
    void clicked_f(const QModelIndex& index);
};

#endif // BLOCKSETWIDGET_HPP

