#include <boost/foreach.hpp>
#include <QtGui>

#include "BlockSetWidget.hpp"
#include "ui_BlockSetWidget.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "block_stat.hpp"
#include "FragmentCollection.hpp"
#include "block_set_alignment.hpp"
#include "block_hash.hpp"
#include "Connector.hpp"
#include "throw_assert.hpp"

enum {
    FRAGMENTS_C, COLUMNS_C,
    IDENTITY_C, GC_C
};

typedef std::vector<Fragment*> Fragments;
typedef FragmentCollection<Fragment*, Fragments> S2F;

struct SeqCmp {
    bool operator()(const Sequence* s1, const Sequence* s2) const {
        return s1->name() < s2->name();
    }
};

class BlockSetModel : public QAbstractTableModel {
public:
    explicit BlockSetModel(QObject* parent = 0):
        QAbstractTableModel(parent) {
        columns_ << tr("fragments") << tr("columns");
        columns_ << tr("% identity") << tr("% GC");
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
                if (index.column() == IDENTITY_C) {
                    double id = block_identity(*stat) * 100;
                    return QString::number(id, 'f', 1);
                } else if (index.column() == GC_C) {
                    double gc = stat->gc() * 100;
                    return QString::number(gc, 'f', 1);
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

class BSAModel : public QAbstractTableModel {
public:
    explicit BSAModel(QObject* parent = 0):
        QAbstractTableModel(parent) {
    }

    QVariant data(const QModelIndex& index,
                  int role = Qt::DisplayRole) const {
        if (!block_set_) {
            return QVariant();
        }
        if (role == Qt::DisplayRole) {
            const BSRow& bsrow = index2bsrow(index);
            Fragment* fragment = index2fragment(index);
            if (fragment) {
                Block* block = fragment->block();
                QString str = QString::number(block->size()) + "x";
                str += QString::number(block->alignment_length());
                int ori = fragment->ori() * bsrow.ori;
                str += " ";
                str += (ori == 1) ? ">" : "<";
                return str;
            } else {
                return "-";
            }
        } else if (role == Qt::BackgroundRole) {
            Fragment* fragment = index2fragment(index);
            if (fragment) {
                Block* block = fragment->block();
                uint32_t hash = block_hash(block);
                hash |= 0xFF808080; // alpha = FF, first bit = 1
                return QColor(QRgb(hash));
            }
        }
        return QVariant();
    }

    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const {
        if (!block_set_) {
            return QVariant();
        }
        if (role == Qt::DisplayRole && orientation == Qt::Vertical) {
            Sequence* seq = bsa2seqs_[bsa_name_][section];
            const BSRow& bsrow = block_set_->bsa(bsa_name_)[seq];
            std::string ori = (bsrow.ori == 1) ? "+" : "-";
            return QString::fromStdString(ori + seq->name());
        }
        return QAbstractTableModel::headerData(section, orientation,
                                               role);
    }

    int rowCount(const QModelIndex& parent = QModelIndex()) const {
        if (!block_set_) {
            return 0;
        }
        return bsa2seqs_[bsa_name_].size();
    }

    int columnCount(const QModelIndex& parent = QModelIndex()) const {
        if (!block_set_) {
            return 0;
        }
        return bsa_length(block_set_->bsa(bsa_name_));
    }

    std::string seq2bsa(Sequence* seq) const {
        if (!block_set_) {
            return "";
        }
        return seq2bsa_[seq];
    }

    const BSRow& index2bsrow(const QModelIndex& index) const {
        Sequence* seq = bsa2seqs_[bsa_name_][index.row()];
        return block_set_->bsa(bsa_name_)[seq];
    }

    Fragment* index2fragment(const QModelIndex& index) const {
        Sequence* seq = bsa2seqs_[bsa_name_][index.row()];
        const BSRow& bsrow = block_set_->bsa(bsa_name_)[seq];
        Fragment* fragment = bsrow.fragments[index.column()];
        return fragment;
    }

    QModelIndex fragment2index(Fragment* fragment) const {
        Sequence* seq = fragment->seq();
        int row = seq2int_[seq];
        int column = fragment2int_[fragment];
        return index(row, column);
    }

public slots:
    void set_block_set(BlockSetPtr block_set) {
        beginResetModel();
        block_set_ = block_set;
        bsa2seqs_.clear();
        seq2int_.clear();
        fragment2int_.clear();
        seq2bsa_.clear();
        bsa_name_ = "";
        BOOST_FOREACH (std::string bsa_name, block_set->bsas()) {
            bsa_name_ = bsa_name; // one of
            const BSA& bsa = block_set->bsa(bsa_name);
            Seqs& seqs = bsa2seqs_[bsa_name];
            BOOST_FOREACH (const BSA::value_type& seq_and_row, bsa) {
                Sequence* seq = seq_and_row.first;
                seqs.push_back(seq);
                seq2bsa_[seq] = bsa_name;
                const BSRow& bsrow = seq_and_row.second;
                for (int i = 0; i < bsrow.fragments.size(); i++) {
                    Fragment* fragment = bsrow.fragments[i];
                    fragment2int_[fragment] = i;
                }
            }
            std::sort(seqs.begin(), seqs.end(), SeqCmp());
            for (int seq_i = 0; seq_i < seqs.size(); seq_i++) {
                Sequence* seq = seqs[seq_i];
                seq2int_[seq] = seq_i;
            }
        }
        endResetModel();
    }

    void set_bsa(const std::string& bsa_name) {
        if (bsa_name != bsa_name_) {
            beginResetModel();
            bsa_name_ = bsa_name;
            endResetModel();
        }
    }

private:
    BlockSetPtr block_set_;
    std::string bsa_name_;
    typedef std::vector<Sequence*> Seqs;
    typedef std::map<std::string, Seqs> Bsa2Seqs;
    mutable Bsa2Seqs bsa2seqs_;
    typedef std::map<Sequence*, int> Seq2Int;
    mutable Seq2Int seq2int_;
    typedef std::map<Fragment*, int> Fragment2Int;
    mutable Fragment2Int fragment2int_;
    typedef std::map<Sequence*, std::string> Seq2Bsa;
    mutable Seq2Bsa seq2bsa_;
};

BlockSetWidget::BlockSetWidget(BlockSetPtr block_set, QWidget* parent) :
    QWidget(parent),
    ui(new Ui::BlockSetWidget) {
    ui->setupUi(this);
    alignment_view_ = new AlignmentView(this);
    alignment_model_ = new AlignmentModel(0, this);
    alignment_view_->set_model(alignment_model_);
    ui->AlignmentView_layout->addWidget(alignment_view_);
    block_set_model_ = new BlockSetModel(this);
    proxy_model_ = new QSortFilterProxyModel(this);
    proxy_model_->setSourceModel(block_set_model_);
    ui->blocksetview->setModel(proxy_model_);
    ui->blocksetview->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->blocksetview->setSortingEnabled(true);
    ui->blocksetview->horizontalHeader()
    ->setResizeMode(QHeaderView::Stretch);
    ui->blocksetview->horizontalHeader()->setMinimumSectionSize(40);
    bsa_model_ = new BSAModel(this);
    ui->bsaView->setModel(bsa_model_);
    connect(bsa_model_, SIGNAL(modelReset()),
            ui->bsaView, SLOT(resizeColumnsToContents()));
    QHeaderView* vh = ui->bsaView->verticalHeader();
    vh->setResizeMode(QHeaderView::Fixed);
    vh->setDefaultSectionSize(vh->fontInfo().pixelSize() + 5);
    ui->bsaView->setModel(bsa_model_);
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
    connect(alignment_view_, SIGNAL(fragment_selected(Fragment*, int)),
            this, SLOT(fragment_selected_f(Fragment*, int)));
    connect(ui->bsaView->selectionModel(),
            SIGNAL(currentChanged(QModelIndex, QModelIndex)),
            this, SLOT(bsa_clicked(QModelIndex)));
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
    bsa_model_->set_block_set(block_set);
    ui->bsaComboBox->clear();
    BOOST_FOREACH (std::string bsa_name, block_set->bsas()) {
        ui->bsaComboBox->addItem(QString::fromStdString(bsa_name));
    }
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

void BlockSetWidget::bsa_clicked(const QModelIndex& index) {
    Fragment* fragment = bsa_model_->index2fragment(index);
    if (fragment) {
        set_block(fragment->block());
        alignment_view_->select_fragment(fragment);
    }
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

void BlockSetWidget::fragment_selected_f(Fragment* fragment, int col) {
    Sequence* seq = fragment->seq();
    std::string bsa_name = bsa_model_->seq2bsa(seq);
    bsa_model_->set_bsa(bsa_name);
    QComboBox* cb = ui->bsaComboBox;
    int row = cb->findText(QString::fromStdString(bsa_name));
    cb->setCurrentIndex(row);
    QModelIndex index = bsa_model_->fragment2index(fragment);
    ui->bsaView->selectionModel()->clearSelection();
    ui->bsaView->selectionModel()->select(index,
                                          QItemSelectionModel::Select);
    ui->bsaView->scrollTo(index);
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

