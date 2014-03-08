#include <QtGui>

#include "BlockSetWidget.hpp"
#include "ui_BlockSetWidget.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "block_stat.hpp"

enum {
    FRAGMENTS_C, COLUMNS_C,
    IDENT_NOGAP_C, IDENT_GAP_C,
    NOIDENT_NOGAP_C, NOIDENT_GAP_C,
    PURE_GAP_C,
    IDENTITY_C, GC_C
};

class BlockSetModel : public QAbstractTableModel {
public:
    explicit BlockSetModel(QObject* parent = 0):
        QAbstractTableModel(parent) {
        columns_ << tr("fragments") << tr("columns");
        columns_ << tr("ident-nogap") << tr("ident-gap");
        columns_ << tr("noident-nogap") << tr("noident-gap");
        columns_ << tr("pure gap");
        columns_ << tr("identity") << tr("GC");
    }

    QVariant data(const QModelIndex& index,
                  int role = Qt::DisplayRole) const {
        if (role == Qt::DisplayRole) {
            const Block* block = blocks_[index.row()];
            if (index.column() == FRAGMENTS_C) {
                return int(block->size());
            } else if (index.column() == COLUMNS_C) {
                return int(block->alignment_length());
            } else {
                AlignmentStat* stat = stats_[index.row()];
                if (stat == 0) {
                    stat = new AlignmentStat;
                    stats_[index.row()] = stat;
                    make_stat(*stat, block);
                }
                if (index.column() == IDENT_NOGAP_C) {
                    return stat->ident_nogap();
                } else if (index.column() == IDENT_GAP_C) {
                    return stat->ident_gap();
                } else if (index.column() == NOIDENT_NOGAP_C) {
                    return stat->noident_nogap();
                } else if (index.column() == NOIDENT_GAP_C) {
                    return stat->noident_gap();
                } else if (index.column() == PURE_GAP_C) {
                    return stat->pure_gap();
                } else if (index.column() == IDENTITY_C) {
                    return block_identity(*stat);
                } else if (index.column() == GC_C) {
                    return stat->gc();
                }
            }
        }
        return QVariant();
    }

    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const {
        if (role == Qt::DisplayRole) {
            if (orientation == Qt::Vertical) {
                return QString::fromStdString(blocks_[section]->name());
            } else if (orientation == Qt::Horizontal) {
                return columns_[section];
            }
        }
        return QAbstractTableModel::headerData(section, orientation, role);
    }

    int rowCount(const QModelIndex& parent = QModelIndex()) const {
        return blocks_.size();
    }

    int columnCount(const QModelIndex& parent = QModelIndex()) const {
        return columns_.size();
    }

    const Block* block_at(int row) const {
        return blocks_[row];
    }

public slots:
    void set_block_set(BlockSetPtr block_set) {
        beginResetModel();
        block_set_ = block_set;
        if (block_set_) {
            std::vector<const Block*> blocks(block_set_->begin(),
                                             block_set_->end());
            blocks_.swap(blocks);
        } else {
            blocks_.clear();
        }
        stats_.clear();
        stats_.resize(blocks_.size(), 0);
        endResetModel();
    }

private:
    BlockSetPtr block_set_;
    std::vector<const Block*> blocks_;
    mutable std::vector<AlignmentStat*> stats_;
    QStringList columns_;
};

BlockSetWidget::BlockSetWidget(BlockSetPtr block_set, QWidget* parent) :
    QWidget(parent),
    ui(new Ui::BlockSetWidget) {
    ui->setupUi(this);
    alignment_view_ = new AlignmentView(this);
    alignment_model_ = new AlignmentModel(0, this);
    alignment_view_->setModel(alignment_model_);
    ui->BlockSetView_layout->addWidget(alignment_view_);
    QSplitter* splitter = new QSplitter(Qt::Vertical, this);
    ui->BlockSetView_layout->addWidget(splitter);
    splitter->addWidget(ui->blocksetview);
    splitter->addWidget(alignment_view_);
    block_set_model_ = new BlockSetModel(this);
    ui->blocksetview->setModel(block_set_model_);
    ui->blocksetview->setSelectionBehavior(QAbstractItemView::SelectRows);
    set_block_set(block_set);
    connect(ui->blocksetview, SIGNAL(clicked(QModelIndex)),
            this, SLOT(clicked_f(QModelIndex)));
    connect(ui->blocksetview->verticalHeader(), SIGNAL(sectionClicked(int)),
            this, SLOT(clicked_h(int)));
}

BlockSetWidget::~BlockSetWidget() {
    delete ui;
}

void BlockSetWidget::set_block_set(BlockSetPtr block_set) {
    block_set_model_->set_block_set(block_set);
}

void BlockSetWidget::clicked_f(const QModelIndex& index) {
    clicked_h(index.row());
}

void BlockSetWidget::clicked_h(int section) {
    const Block* block = block_set_model_->block_at(section);
    alignment_model_->set_block(block);
}

