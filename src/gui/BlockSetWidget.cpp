#include <boost/foreach.hpp>
#include <QtGui>

#include "BlockSetWidget.hpp"
#include "ui_BlockSetWidget.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "block_stat.hpp"
#include "FragmentCollection.hpp"
#include "Connector.hpp"
#include "throw_assert.hpp"

enum {
    FRAGMENTS_C, COLUMNS_C,
    IDENT_NOGAP_C, IDENT_GAP_C,
    NOIDENT_NOGAP_C, NOIDENT_GAP_C,
    PURE_GAP_C,
    IDENTITY_C, GC_C
};

typedef std::vector<Fragment*> Fragments;
typedef FragmentCollection<Fragment*, Fragments> S2F;

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
        if (role == Qt::UserRole && index.column() == FRAGMENTS_C) {
            // for filter
            int section = index.row();
            //qDebug() << section << ' ' << blocks_.size();
            return QString::fromStdString(blocks_[section]->name());
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

    int block_index(const Block* block) const {
        for (int i = 0; i < blocks_.size(); i++) {
            if (blocks_[i] == block) {
                return i;
            }
        }
        return -1;
    }

    const QPoint& xy_of(int row) const {
        return alignment_xy_[row];
    }

    void set_xy_of(int row, const QPoint& xy) {
        alignment_xy_[row] = xy;
    }

public slots:
    void set_block_set(BlockSetPtr block_set) {
        beginResetModel();
        block_set_ = block_set;
        Connector c;
        c.apply(block_set_);
        if (block_set_) {
            std::vector<const Block*> blocks(block_set_->begin(),
                                             block_set_->end());
            blocks_.swap(blocks);
        } else {
            blocks_.clear();
        }
        stats_.clear();
        stats_.resize(blocks_.size(), 0);
        alignment_xy_.clear();
        alignment_xy_.resize(blocks_.size());
        endResetModel();
    }

    void set_genes(BlockSetPtr genes) {
        genes_ = genes;
        genes_s2f_.clear();
        if (genes_) {
            genes_s2f_.add_bs(*genes_);
        }
        stats_.clear();
        stats_.resize(blocks_.size(), 0);
        alignment_xy_.clear();
        alignment_xy_.resize(blocks_.size());
    }

    void find_genes(std::vector<Fragment*>& overlap_genes,
                    Fragment* f) const {
        genes_s2f_.find_overlap_fragments(overlap_genes, f);
    }

private:
    BlockSetPtr block_set_;
    std::vector<const Block*> blocks_;
    mutable std::vector<AlignmentStat*> stats_;
    mutable std::vector<QPoint> alignment_xy_;
    QStringList columns_;
    BlockSetPtr genes_;
    S2F genes_s2f_;
};

BlockSetWidget::BlockSetWidget(BlockSetPtr block_set, QWidget* parent) :
    QWidget(parent),
    ui(new Ui::BlockSetWidget) {
    ui->setupUi(this);
    alignment_view_ = new AlignmentView(this);
    alignment_model_ = new AlignmentModel(0, this);
    alignment_view_->setModel(alignment_model_);
    ui->AlignmentView_layout->addWidget(alignment_view_);
    block_set_model_ = new BlockSetModel(this);
    proxy_model_ = new QSortFilterProxyModel(this);
    proxy_model_->setSourceModel(block_set_model_);
    ui->blocksetview->setModel(proxy_model_);
    ui->blocksetview->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->blocksetview->setSortingEnabled(true);
    set_block_set(block_set);
    connect(ui->blocksetview->selectionModel(),
            SIGNAL(currentChanged(QModelIndex, QModelIndex)),
            this, SLOT(clicked_f(QModelIndex)));
    connect(alignment_model_, SIGNAL(modelReset()),
            this, SLOT(update_gene_layout()));
    connect(ui->showGenesCheckBox, SIGNAL(clicked(bool)),
            alignment_model_, SLOT(set_show_genes(bool)));
    alignment_model_->set_show_genes(ui->showGenesCheckBox->isChecked());
    update_gene_layout();
    connect(alignment_view_, SIGNAL(clicked(QModelIndex)),
            this, SLOT(alignment_clicked(QModelIndex)));
    connect(alignment_view_->selectionModel(),
            SIGNAL(currentChanged(QModelIndex, QModelIndex)),
            this, SLOT(alignment_clicked(QModelIndex)));
    connect(alignment_view_, SIGNAL(jump_to(Fragment*, int)),
            this, SLOT(jump_to_f(Fragment*, int)));
    ui->blocksetview->addAction(ui->actionCopy_block_name);
    ui->blocksetview->setContextMenuPolicy(Qt::ActionsContextMenu);
    alignment_view_->addAction(ui->actionCopy_fragment_id);
    alignment_view_->setContextMenuPolicy(Qt::ActionsContextMenu);
}

BlockSetWidget::~BlockSetWidget() {
    delete ui;
}

void BlockSetWidget::set_block_set(BlockSetPtr block_set) {
    block_set_model_->set_block_set(block_set);
    prev_row_ = -1;
}

void BlockSetWidget::set_genes(BlockSetPtr genes) {
    block_set_model_->set_genes(genes);
}

void BlockSetWidget::set_block(const Block* block) {
    if (prev_row_ != -1) {
        int col = alignment_view_->columnAt(0);
        int row = alignment_view_->rowAt(0);
        block_set_model_->set_xy_of(prev_row_, QPoint(col, row));
        const Block* prev_block = block_set_model_->block_at(prev_row_);
        fragments_[prev_block] = alignment_model_->fragments();
    }
    int section = block_set_model_->block_index(block);
    QModelIndex index = block_set_model_->index(section, 0);
    QModelIndex index_in_proxy = proxy_model_->mapFromSource(index);
    ui->blocksetview->selectionModel()->clearSelection();
    ui->blocksetview->selectionModel()->select(index_in_proxy,
            QItemSelectionModel::Select |
            QItemSelectionModel::Rows);
    ui->blocksetview->scrollTo(index_in_proxy);
    alignment_model_->set_block(block);
    if (fragments_.find(block) != fragments_.end()) {
        alignment_model_->set_fragments(fragments_[block]);
    }
    QPoint xy = block_set_model_->xy_of(section);
    QModelIndex rb, target;
    rb = alignment_model_->index(alignment_model_->rowCount() - 1,
                                 alignment_model_->columnCount() - 1);
    target = alignment_model_->index(xy.y(), xy.x());
    alignment_view_->scrollTo(rb);
    alignment_view_->scrollTo(target);
    prev_row_ = section;
    // genes
    BOOST_FOREACH (Fragment* f, *block) {
        std::vector<Fragment*> overlap_genes;
        block_set_model_->find_genes(overlap_genes, f);
        alignment_model_->add_genes(f, overlap_genes);
    }
    ui->geneNameLineEdit->setText("");
}

void BlockSetWidget::clicked_f(const QModelIndex& index) {
    int section = proxy_model_->mapToSource(index).row();
    const Block* block = block_set_model_->block_at(section);
    set_block(block);
}

void BlockSetWidget::jump_to_f(Fragment* fragment, int col) {
    BOOST_ASSERT(fragment->block());
    set_block(fragment->block());
    int row = alignment_model_->fragment_index(fragment);
    QModelIndex index = alignment_model_->index(row, col);
    alignment_view_->selectionModel()->clearSelection();
    alignment_view_->setCurrentIndex(index);
    alignment_view_->scrollTo(index);
}

void BlockSetWidget::on_nonunique_stateChanged(int state) {
    if (state == Qt::Checked) {
        proxy_model_->setFilterRegExp(QRegExp("[^1]|.{2,}"));
        proxy_model_->setFilterKeyColumn(FRAGMENTS_C);
        proxy_model_->setFilterRole(Qt::DisplayRole);
    } else {
        proxy_model_->setFilterRegExp(QRegExp(""));
    }
}

void BlockSetWidget::update_gene_layout() {
    ui->genesWidget->setHidden(!alignment_model_->has_genes());
    ui->geneNameLineEdit->setText("");
}

void BlockSetWidget::alignment_clicked(const QModelIndex& index) {
    QVariant tip = alignment_model_->data(index, Qt::ToolTipRole);
    ui->geneNameLineEdit->setText(tip.toString());
}

void BlockSetWidget::on_blockNameLineEdit_editingFinished() {
    QString pattern = ui->blockNameLineEdit->text();
    if (pattern.isEmpty()) {
        // re-enable filter by number of fragments
        on_nonunique_stateChanged(ui->nonunique->checkState());
    } else {
        proxy_model_->setFilterWildcard(pattern);
        proxy_model_->setFilterKeyColumn(FRAGMENTS_C);
        proxy_model_->setFilterRole(Qt::UserRole);
    }
}

void BlockSetWidget::on_actionCopy_block_name_triggered() {
    QModelIndex index = ui->blocksetview->currentIndex();
    if (!index.isValid()) {
        return;
    }
    int section = index.row();
    QString name = proxy_model_->headerData(section, Qt::Vertical).toString();
    QApplication::clipboard()->setText(name);
}

void BlockSetWidget::on_actionCopy_fragment_id_triggered() {
    QModelIndex index = alignment_view_->currentIndex();
    if (!index.isValid()) {
        return;
    }
    int section = index.row();
    QString name = alignment_model_->headerData(section,
                   Qt::Vertical).toString();
    QApplication::clipboard()->setText(name);
}

