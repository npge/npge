#include <set>
#include <boost/bind.hpp>
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
#include "move_rows.hpp"
#include "throw_assert.hpp"

enum {
    FRAGMENTS_C, COLUMNS_C,
    IDENTITY_C, GC_C,
    GENES_C, SPLIT_C, LOW_C
};

typedef std::vector<Fragment*> S2F_Fragments;
typedef FragmentCollection<Fragment*, S2F_Fragments> S2F;

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
        columns_ << tr("genes") << tr("split parts");
        columns_ << tr("low similarity regions");
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
                    int id = block_identity(*stat) * 1000;
                    return double(id) / 10.0;
                } else if (index.column() == GC_C) {
                    int gc = stat->gc() * 1000;
                    return double(gc) / 10.0;
                } else if (index.column() == GENES_C) {
                    Fragments genes;
                    find_genes(genes, block);
                    return int(genes.size());
                } else if (index.column() == SPLIT_C) {
                    Blocks bb;
                    find_split_parts(bb, block);
                    return int(bb.size());
                } else if (index.column() == LOW_C) {
                    Blocks bb;
                    find_low_similarity(bb, block);
                    return int(bb.size());
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
        find_first_last();
        endResetModel();
    }

    void set_genes(BlockSetPtr genes) {
        genes_ = genes;
        Connector c;
        c.set_opt_value("connect-circular", true);
        c.apply(genes_);
        genes_s2f_.clear();
        if (genes_) {
            genes_s2f_.add_bs(*genes_);
        }
        stats_.clear();
        stats_.resize(blocks_.size(), 0);
        alignment_xy_.clear();
        alignment_xy_.resize(blocks_.size());
    }

    void find_genes(Fragments& overlap_genes, Fragment* f) const {
        genes_s2f_.find_overlap_fragments(overlap_genes, f);
    }

    void find_genes(Fragments& overlap_genes,
                    const Block* block) const {
        BOOST_FOREACH (Fragment* f, *block) {
            find_genes(overlap_genes, f);
        }
    }

    void set_split_parts(BlockSetPtr split_parts) {
        split_parts_ = split_parts;
        split_s2f_.clear();
        if (split_parts_) {
            split_s2f_.add_bs(*split_parts_);
        }
    }

    void find_split_parts(Fragments& ff, Fragment* f) const {
        split_s2f_.find_overlap_fragments(ff, f);
    }

    void find_split_parts(Blocks& bb, const Block* block) const {
        std::set<Block*> split_parts_set;
        BOOST_FOREACH (Fragment* f, *block) {
            Fragments ff;
            find_split_parts(ff, f);
            BOOST_FOREACH (Fragment* f1, ff) {
                split_parts_set.insert(f1->block());
            }
        }
        BOOST_FOREACH (Block* b, split_parts_set) {
            bb.push_back(b);
        }
    }

    void set_low_similarity(BlockSetPtr low_similarity) {
        low_similarity_ = low_similarity;
        low_s2f_.clear();
        if (low_similarity_) {
            low_s2f_.add_bs(*low_similarity_);
        }
    }

    void find_low_similarity(Fragments& ff, Fragment* f) const {
        low_s2f_.find_overlap_fragments(ff, f);
    }

    void find_low_similarity(Blocks& bb, const Block* block) const {
        std::set<Block*> low_similarity_set;
        BOOST_FOREACH (Fragment* f, *block) {
            Fragments ff;
            find_low_similarity(ff, f);
            BOOST_FOREACH (Fragment* f1, ff) {
                low_similarity_set.insert(f1->block());
            }
        }
        BOOST_FOREACH (Block* b, low_similarity_set) {
            bb.push_back(b);
        }
    }

    void find_first_last() {
        seq2first_.clear();
        seq2last_.clear();
        BOOST_FOREACH (Block* block, *block_set_) {
            BOOST_FOREACH (Fragment* f, *block) {
                Sequence* seq = f->seq();
                BOOST_ASSERT(seq);
                if (seq2first_.find(seq) == seq2first_.end()) {
                    seq2first_[seq] = f;
                    seq2last_[seq] = f;
                } else {
                    if (*f < * (seq2first_[seq])) {
                        seq2first_[seq] = f;
                    }
                    if (*(seq2last_[seq]) < *f) {
                        seq2last_[seq] = f;
                    }
                }
            }
        }
    }

    typedef std::map<Sequence*, Fragment*> Seq2Fragment;

    const Seq2Fragment& seq2first() const {
        return seq2first_;
    }

    const Seq2Fragment& seq2last() const {
        return seq2last_;
    }

private:
    BlockSetPtr block_set_;
    std::vector<const Block*> blocks_;
    mutable std::vector<AlignmentStat*> stats_;
    mutable std::vector<QPoint> alignment_xy_;
    QStringList columns_;
    BlockSetPtr genes_;
    BlockSetPtr split_parts_;
    BlockSetPtr low_similarity_;
    S2F genes_s2f_;
    S2F split_s2f_;
    S2F low_s2f_;
    mutable Seq2Fragment seq2first_;
    mutable Seq2Fragment seq2last_;
};

class BSAModel : public QAbstractTableModel {
public:
    explicit BSAModel(QObject* parent = 0):
        QAbstractTableModel(parent) {
    }

    QVariant data(const QModelIndex& index,
                  int role = Qt::DisplayRole) const {
        if (!block_set_ || !block_set_->has_bsa(bsa_name_)) {
            return QVariant();
        }
        if (role == Qt::TextAlignmentRole) {
            return Qt::AlignCenter;
        }
        if (role == Qt::DisplayRole) {
            const BSRow& bsrow = index2bsrow(index);
            Fragment* fragment = index2fragment(index);
            if (fragment) {
                Block* block = fragment->block();
                QString str = QString::number(block->size()) + "x";
                str += QString::number(fragment->length());
                int ori = fragment->ori() * bsrow.ori;
                str += " ";
                str += (ori == 1) ? ">" : "<";
                int end_ori = (!fragment->prev()) ? -1 :
                              (!fragment->next()) ? 1 : 0;
                end_ori *= bsrow.ori;
                if (end_ori == -1) {
                    str = "[ " + str;
                } else if (end_ori == 1) {
                    str += " ]";
                }
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
        if (!block_set_ || !block_set_->has_bsa(bsa_name_)) {
            return QVariant();
        }
        if (role == Qt::DisplayRole && orientation == Qt::Vertical) {
            Sequence* seq = bsa2seqs_[bsa_name_][section];
            if (seq) {
                const BSRow& bsrow = block_set_->bsa(bsa_name_)[seq];
                std::string ori = (bsrow.ori == 1) ? "+" : "-";
                return QString::fromStdString(ori + seq->name());
            }
        }
        return QAbstractTableModel::headerData(section, orientation,
                                               role);
    }

    int rowCount(const QModelIndex& parent = QModelIndex()) const {
        if (!block_set_ || !block_set_->has_bsa(bsa_name_)) {
            return 0;
        }
        return bsa2seqs_[bsa_name_].size();
    }

    int columnCount(const QModelIndex& parent = QModelIndex()) const {
        if (!block_set_ || !block_set_->has_bsa(bsa_name_)) {
            return 0;
        }
        return bsa_length(block_set_->bsa(bsa_name_));
    }

    std::string seq2bsa(Sequence* seq) const {
        if (!block_set_ || !block_set_->has_bsa(bsa_name_)) {
            return "";
        }
        return seq2bsa_[seq];
    }

    const BSRow& index2bsrow(const QModelIndex& index) const {
        Sequence* seq = bsa2seqs_[bsa_name_][index.row()];
        return block_set_->bsa(bsa_name_)[seq];
    }

    Fragment* index2fragment(const QModelIndex& index) const {
        if (!block_set_ || !block_set_->has_bsa(bsa_name_)) {
            return 0;
        }
        if (!index.isValid()) {
            return 0;
        }
        Sequence* seq = bsa2seqs_[bsa_name_][index.row()];
        if (!seq) {
            return 0;
        }
        const BSRow& bsrow = block_set_->bsa(bsa_name_)[seq];
        if (index.column() >= bsrow.fragments.size()) {
            return 0;
        }
        Fragment* fragment = bsrow.fragments[index.column()];
        return fragment;
    }

    QModelIndex fragment2index(Fragment* fragment) const {
        Sequence* seq = fragment->seq();
        Seq2Int::const_iterator it = seq2int_.find(seq);
        if (it == seq2int_.end()) {
            return QModelIndex();
        }
        int row = it->second;
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
        BOOST_FOREACH (std::string bsa_name, block_set_->bsas()) {
            const BSA& bsa = block_set_->bsa(bsa_name);
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
            update_seq2int(bsa_name);
        }
        endResetModel();
        if (!block_set_->bsas().empty()) {
            set_bsa(block_set_->bsas()[0]);
        }
    }

    void set_bsa(const std::string& bsa_name) {
        if (bsa_name != bsa_name_ &&
                bsa2seqs_.find(bsa_name) != bsa2seqs_.end()) {
            beginResetModel();
            bsa_name_ = bsa_name;
            endResetModel();
        }
    }

    void update_seq2int(const std::string& bsa_name) {
        const BSA& bsa = block_set_->bsa(bsa_name);
        const Seqs& seqs = bsa2seqs_[bsa_name];
        for (int seq_i = 0; seq_i < seqs.size(); seq_i++) {
            Sequence* seq = seqs[seq_i];
            seq2int_[seq] = seq_i;
        }
    }

    void move_seqs(std::vector<int>& rows, bool up) {
        beginResetModel();
        Seqs& seqs = bsa2seqs_[bsa_name_];
        if (up) {
            int prev_row = -100;
            BOOST_FOREACH (int& row, rows) {
                if (row > 0 && row - 1 != prev_row) {
                    std::swap(seqs[row], seqs[row - 1]);
                    row -= 1;
                }
                prev_row = row;
            }
        } else {
            int prev_row = -100;
            BOOST_REVERSE_FOREACH (int& row, rows) {
                if (row < seqs.size() - 1 && row + 1 != prev_row) {
                    std::swap(seqs[row], seqs[row + 1]);
                    row += 1;
                }
                prev_row = row;
            }
        }
        update_seq2int(bsa_name_);
        endResetModel();
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

class BSAView : public QTableView {
public:
    BSAView(QWidget* parent):
        QTableView(parent) {
        QHeaderView* vh = verticalHeader();
        vh->setResizeMode(QHeaderView::Fixed);
        vh->setDefaultSectionSize(vh->fontInfo().pixelSize() + 5);
        QHeaderView* hh = horizontalHeader();
        hh->setResizeMode(QHeaderView::Fixed);
        hh->setDefaultSectionSize(100);
    }

    void keyPressEvent(QKeyEvent* e) {
        bool ctrl = e->modifiers().testFlag(Qt::ControlModifier);
        bool up_down = e->key() == Qt::Key_Up ||
                       e->key() == Qt::Key_Down;
        if (ctrl && up_down) {
            BSAModel* m = dynamic_cast<BSAModel*>(model());
            BOOST_ASSERT(m);
            move_view_rows(this, e->key() == Qt::Key_Up,
                           boost::bind(&BSAModel::move_seqs,
                                       m, _1, _2));
        } else {
            QTableView::keyPressEvent(e);
        }
    }
};

BlockSetWidget::BlockSetWidget(BlockSetPtr block_set, QWidget* parent) :
    QWidget(parent),
    ui(new Ui::BlockSetWidget) {
    ui->setupUi(this);
    //
    ui->splitter_2->setStretchFactor(0, 1);
    ui->splitter_2->setStretchFactor(1, 3);
    //
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
    bsa_view_ = new BSAView(this);
    ui->BSAView_layout->addWidget(bsa_view_);
    bsa_view_->setModel(bsa_model_);
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
    connect(bsa_view_->selectionModel(),
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
    alignment_view_->set_first_last(block_set_model_->seq2first(),
                                    block_set_model_->seq2last());
}

void BlockSetWidget::set_genes(BlockSetPtr genes) {
    block_set_model_->set_genes(genes);
}

void BlockSetWidget::set_split_parts(BlockSetPtr split_parts) {
    block_set_model_->set_split_parts(split_parts);
}

void BlockSetWidget::set_low_similarity(BlockSetPtr low_similarity) {
    block_set_model_->set_low_similarity(low_similarity);
}

void BlockSetWidget::set_block(const Block* block) {
    if (alignment_model_->block() == block) {
        return;
    }
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
    if (index_in_proxy.isValid()) {
        ui->blocksetview->selectionModel()->clearSelection();
        ui->blocksetview->selectionModel()->select(index_in_proxy,
                QItemSelectionModel::Select |
                QItemSelectionModel::Rows);
        ui->blocksetview->scrollTo(index_in_proxy);
    }
    alignment_model_->set_block(block);
    QPoint xy = block_set_model_->xy_of(section);
    QModelIndex rb, target;
    rb = alignment_model_->index(alignment_model_->rowCount() - 1,
                                 alignment_model_->columnCount() - 1);
    target = alignment_model_->index(xy.y(), xy.x());
    alignment_view_->scrollTo(rb);
    alignment_view_->scrollTo(target);
    prev_row_ = section;
    // split_parts
    Blocks split_parts;
    block_set_model_->find_split_parts(split_parts, block);
    alignment_model_->set_split_parts(split_parts);
    // low_similarity
    Blocks low_similarity;
    block_set_model_->find_low_similarity(low_similarity, block);
    alignment_model_->set_low_similarity(low_similarity);
    //
    if (fragments_.find(block) != fragments_.end()) {
        alignment_model_->set_fragments(fragments_[block]);
    }
    // genes
    BOOST_FOREACH (Fragment* f, *block) {
        Fragments overlap_genes;
        block_set_model_->find_genes(overlap_genes, f);
        alignment_model_->add_genes(f, overlap_genes);
    }
    //
    ui->geneNameLineEdit->setText("");
}

void BlockSetWidget::clicked_f(const QModelIndex& index) {
    QModelIndex index_in_proxy = proxy_model_->mapToSource(index);
    if (index_in_proxy.isValid()) {
        int section = proxy_model_->mapToSource(index).row();
        const Block* block = block_set_model_->block_at(section);
        set_block(block);
    }
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
    if (index.isValid()) {
        QItemSelectionModel* sm = bsa_view_->selectionModel();
        sm->clearSelection();
        sm->select(index, QItemSelectionModel::Select);
        bsa_view_->scrollTo(index);
    }
}

void BlockSetWidget::on_bsaComboBox_activated(QString bsa_name) {
    bsa_model_->set_bsa(bsa_name.toStdString());
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

