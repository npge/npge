/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <boost/bind.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <QtGui>

#include "BlockSetWidget.hpp"
#include "BlockSetModel.hpp"
#include "ui_BlockSetWidget.h"
#include "AlignmentView.hpp"
#include "AlignmentModel.hpp"
#include "Meta.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "block_stat.hpp"
#include "block_set_alignment.hpp"
#include "block_hash.hpp"
#include "move_rows.hpp"
#include "Processor.hpp"
#include "BlockSearcher.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "global.hpp"

struct SeqCmp {
    bool operator()(const Sequence* s1,
                    const Sequence* s2) const {
        return s1->name() < s2->name();
    }
};

BlockSetModel::BlockSetModel(QObject* parent):
    QAbstractTableModel(parent), more_than_1_(false) {
    columns_ << tr("fragments") << tr("columns");
    columns_ << tr("% identity") << tr("% GC");
    columns_ << tr("genes") << tr("split parts");
    columns_ << tr("low similarity regions");
    connect(this, SIGNAL(exceptionThrown(QString)),
            this, SLOT(onExceptionThrown(QString)),
            Qt::QueuedConnection);
}

BlockSetPtr BlockSetModel::block_set() const {
    return block_set_;
}

QVariant BlockSetModel::data(const QModelIndex& index,
                             int role) const {
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
                Decimal id = block_identity(*stat);
                return (id * 100).to_d();
            } else if (index.column() == GC_C) {
                Decimal gc = stat->gc();
                return (gc * 100).to_d();
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
        const Block* block = blocks_[index.row()];
        if (filtered_blocks_.has_elem(block)) {
            return QString("yes");
        } else {
            return QString();
        }
    }
    return QVariant();
}

QVariant BlockSetModel::headerData(
    int section, Qt::Orientation orientation,
    int role) const {
    if (role == Qt::DisplayRole) {
        if (orientation == Qt::Vertical) {
            return QString::fromStdString(
                       blocks_[section]->name());
        } else if (orientation == Qt::Horizontal) {
            return columns_[section];
        }
    }
    return QAbstractTableModel::headerData(section,
                                           orientation, role);
}

int BlockSetModel::rowCount(const QModelIndex& parent) const {
    return blocks_.size();
}

int BlockSetModel::columnCount(
    const QModelIndex& parent) const {
    return columns_.size();
}

const Block* BlockSetModel::block_at(int row) const {
    return blocks_[row];
}

int BlockSetModel::block_index(const Block* block) const {
    for (int i = 0; i < blocks_.size(); i++) {
        if (blocks_[i] == block) {
            return i;
        }
    }
    return -1;
}

const QPoint& BlockSetModel::xy_of(int row) const {
    return alignment_xy_[row];
}

void BlockSetModel::set_xy_of(int row, const QPoint& xy) {
    alignment_xy_[row] = xy;
}

void BlockSetModel::set_block_set(BlockSetPtr block_set) {
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
    alignment_xy_.clear();
    alignment_xy_.resize(blocks_.size());
    endResetModel();
}

void BlockSetModel::set_genes(BlockSetPtr genes) {
    genes_ = genes;
    genes_s2f_.clear();
    if (genes_) {
        genes_s2f_.add_bs(*genes_);
        genes_s2f_.prepare();
    }
    stats_.clear();
    stats_.resize(blocks_.size(), 0);
    alignment_xy_.clear();
    alignment_xy_.resize(blocks_.size());
}

void BlockSetModel::find_genes(Fragments& overlap_genes,
                               Fragment* f) const {
    genes_s2f_.find_overlap_fragments(overlap_genes, f);
}

void BlockSetModel::find_genes(Fragments& overlap_genes,
                               const Block* block) const {
    BOOST_FOREACH (Fragment* f, *block) {
        find_genes(overlap_genes, f);
    }
}

Fragments BlockSetModel::return_genes(
    const Block* block) const {
    Fragments genes;
    find_genes(genes, block);
    return genes;
}

void BlockSetModel::set_split_parts(BlockSetPtr split_parts) {
    split_parts_ = split_parts;
    split_s2f_.clear();
    if (split_parts_) {
        split_s2f_.add_bs(*split_parts_);
        split_s2f_.prepare();
    }
}

void BlockSetModel::find_split_parts(
    Fragments& ff, Fragment* f) const {
    split_s2f_.find_overlap_fragments(ff, f);
}

void BlockSetModel::find_split_parts(
    Blocks& bb, const Block* block) const {
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

void BlockSetModel::set_low_similarity(
    BlockSetPtr low_similarity) {
    low_similarity_ = low_similarity;
    low_s2f_.clear();
    if (low_similarity_) {
        low_s2f_.add_bs(*low_similarity_);
        low_s2f_.prepare();
    }
}

void BlockSetModel::find_low_similarity(
    Fragments& ff, Fragment* f) const {
    low_s2f_.find_overlap_fragments(ff, f);
}

void BlockSetModel::find_low_similarity(
    Blocks& bb, const Block* block) const {
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

void BlockSetModel::onExceptionThrown(QString message) {
    throw Exception(message.toStdString());
}

void BlockSetModel::set_more_than_1(bool more_than_1) {
    if (more_than_1 != more_than_1_) {
        more_than_1_ = more_than_1;
        update_filter();
    }
}

void BlockSetModel::set_pattern(const std::string& pattern) {
    if (pattern != pattern_) {
        pattern_ = pattern;
        update_filter();
    }
}

void BlockSetModel::construct_hits() {
    hits_bs_ = new_bs();
    BOOST_FOREACH (char c, pattern_) {
        c = toupper(c);
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C' &&
                c != 'N' && c != '-') {
            return;
        }
    }
    std::string pattern = pattern_;
    Sequence::to_atgcn(pattern);
    if (!pattern.empty()) {
        Meta* m = Meta::instance();
        ASSERT_TRUE(m);
        SharedProcessor finder = m->get("FragmentFinder");
        finder->set_opt_value("pattern", pattern);
        hits_bs_->add_sequences(block_set()->seqs());
        finder->set_block_set(hits_bs_);
        finder->run();
    }
    hits_s2f_.clear();
    hits_s2f_.add_bs(*hits_bs_);
    hits_s2f_.prepare();
}

static bool name_matches(const std::string& name,
                         const std::string& pattern,
                         bool starting_only,
                         bool ending_only) {
    using namespace boost::algorithm;
    if (starting_only && ending_only) {
        return name == pattern;
    } else if (starting_only) {
        return starts_with(name, pattern);
    } else if (ending_only) {
        return ends_with(name, pattern);
    } else {
        return name.find(pattern) != std::string::npos;
    }
}

bool BlockSetModel::check_block(const Block* block) const {
    const size_t npos = std::string::npos;
    // FIXME no const Blocks, Fragments etc here at all!
    Block* b = const_cast<Block*>(block);
    if (block->size() <= 1 && more_than_1_) {
        return false;
    }
    std::string pattern = pattern_;
    bool starting_only = false;
    if (!pattern.empty() && pattern[0] == '^') {
        pattern = pattern.substr(1);
        starting_only = true;
    }
    bool ending_only = false;
    if (!pattern.empty() &&
            pattern[pattern.size() - 1] == '$') {
        pattern = pattern.substr(0, pattern.size() - 1);
        ending_only = true;
    }
    if (pattern.empty()) {
        return true;
    } else if (name_matches(block->name(), pattern,
                            starting_only, ending_only)) {
        return true;
    } else {
        // genes
        BOOST_FOREACH (Fragment* gene, return_genes(block)) {
            ASSERT_TRUE(gene->block());
            std::string gene_name = gene->block()->name();
            if (name_matches(gene_name, pattern,
                             starting_only, ending_only)) {
                return true;
            }
        }
        // sequence hits
        if (hits_s2f_.block_has_overlap(b)) {
            return true;
        }
    }
    return false;
}

void BlockSetModel::update_filter() {
    BlockSearcher* searcher = new BlockSearcher;
    searcher->meta_ = Meta::instance();
    ASSERT_TRUE(searcher->meta_);
    searcher->blocks_ = &blocks_;
    searcher->model_ = this;
    searcher->filtered_blocks_ = &filtered_blocks_;
    connect(searcher, SIGNAL(searchingFinished(QString)),
            this, SLOT(onSearchingFinished(QString)),
            Qt::QueuedConnection);
    QThreadPool::globalInstance()->start(searcher);
    emit searchStarted();
}

void BlockSetModel::onSearchingFinished(QString message) {
    emit searchFinished();
    reset();
    if (!message.isEmpty()) {
        emit exceptionThrown(message);
    }
}

class BSAModel : public QAbstractTableModel {
public:
    explicit BSAModel(QObject* parent = 0):
        QAbstractTableModel(parent) {
    }

    BlockSetPtr block_set() const {
        return block_set_;
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
                QString str;
                if (block->name().size() >= 2 &&
                        isdigit(block->name()[1])) {
                    str += block->name()[0];
                }
                str += QString::number(block->size()) + "x";
                str += QString::number(fragment->length());
                int ori = fragment->ori() * bsrow.ori;
                str += " ";
                str += (ori == 1) ? ">" : "<";
                int end_ori = find_end_ori(fragment);
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
                hash_t hash = block_hash(block);
                hash &= 0x00FFFFFF;
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

    const BSRow& seq2bsrow(Sequence* seq) const {
        return block_set_->bsa(bsa_name_)[seq];
    }

    const BSRow& index2bsrow(const QModelIndex& index) const {
        Sequence* seq = bsa2seqs_[bsa_name_][index.row()];
        return seq2bsrow(seq);
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

    Fragment* first_fragment(Sequence* seq) const {
        const Fragments& ff = s2f_.fragments_of(seq);
        ASSERT_GT(ff.size(), 0);
        return ff[0];
    }

    Fragment* last_fragment(Sequence* seq) const {
        const Fragments& ff = s2f_.fragments_of(seq);
        ASSERT_GT(ff.size(), 0);
        return ff[ff.size() - 1];
    }

    int find_end_ori(Fragment* f) const {
        if (f == first_fragment(f->seq())) {
            return -1;
        }
        if (f == last_fragment(f->seq())) {
            return 1;
        }
        return 0;
    }

public slots:
    void set_block_set(BlockSetPtr block_set) {
        beginResetModel();
        block_set_ = block_set;
        //
        s2f_.clear();
        s2f_.add_bs(*block_set_);
        s2f_.prepare();
        //
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
    VectorFc s2f_;
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
        hh->setDefaultSectionSize(130);
    }

    Block* selected_block() const {
        Block* result = new Block;
        BSAModel* m = dynamic_cast<BSAModel*>(model());
        ASSERT_TRUE(m);
        QItemSelectionModel* sm = selectionModel();
        typedef std::pair<int, int> FL;
        typedef std::map<Sequence*, FL> Seq2FL;
        Seq2FL seq2fl;
        foreach (const QModelIndex& index,
                sm->selectedIndexes()) {
            Fragment* f = m->index2fragment(index);
            if (!f) {
                continue;
            }
            ASSERT_TRUE(f);
            Sequence* seq = f->seq();
            ASSERT_TRUE(seq);
            if (seq2fl.find(seq) == seq2fl.end()) {
                seq2fl[seq] = FL(f->min_pos(), f->max_pos());
            } else {
                int& min_pos = seq2fl[seq].first;
                min_pos = std::min(min_pos, int(f->min_pos()));
                int& max_pos = seq2fl[seq].second;
                max_pos = std::max(max_pos, int(f->max_pos()));
            }
        }
        BOOST_FOREACH (const Seq2FL::value_type& vt, seq2fl) {
            Sequence* seq = vt.first;
            const FL& first_last = vt.second;
            int min_pos = first_last.first;
            int max_pos = first_last.second;
            Fragment* f = new Fragment(seq, min_pos, max_pos);
            const BSRow& bsrow = m->seq2bsrow(seq);
            if (bsrow.ori == -1) {
                f->inverse();
            }
            result->insert(f);
        }
        int genomes = genomes_number(*(m->block_set()));
        set_canonical_name(result, genomes);
        return result;
    }

    void keyPressEvent(QKeyEvent* e) {
        bool ctrl = e->modifiers().testFlag(Qt::ControlModifier);
        bool up_down = e->key() == Qt::Key_Up ||
                       e->key() == Qt::Key_Down;
        bool c = e->key() == Qt::Key_C;
        if (ctrl && up_down) {
            BSAModel* m = dynamic_cast<BSAModel*>(model());
            ASSERT_TRUE(m);
            move_view_rows(this, e->key() == Qt::Key_Up,
                           boost::bind(&BSAModel::move_seqs,
                                       m, _1, _2));
        } else if (ctrl && c) {
            boost::scoped_ptr<Block> s(selected_block());
            std::stringstream ss;
            ss << *s;
            QString text = QString::fromStdString(ss.str());
            QApplication::clipboard()->setText(text);
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
    ui->splitter_2->setStretchFactor(0, 3);
    ui->splitter_2->setStretchFactor(1, 4);
    //
    alignment_view_ = new AlignmentView(this);
    alignment_model_ = new AlignmentModel(0, this);
    alignment_view_->set_model(alignment_model_);
    ui->AlignmentView_layout->addWidget(alignment_view_);
    block_set_model_ = new BlockSetModel(this);
    connect(block_set_model_, SIGNAL(searchStarted()),
            this, SLOT(onSearchStarted()));
    connect(block_set_model_, SIGNAL(searchFinished()),
            this, SLOT(onSearchFinished()));
    proxy_model_ = new QSortFilterProxyModel(this);
    proxy_model_->setSourceModel(block_set_model_);
    proxy_model_->setFilterFixedString("yes");
    proxy_model_->setFilterKeyColumn(FRAGMENTS_C);
    proxy_model_->setFilterRole(Qt::UserRole);
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
    ui->blocksetview->sortByColumn(FRAGMENTS_C);
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
    alignment_model_->set_block_set(block_set);
    ui->bsaComboBox->clear();
    BOOST_FOREACH (std::string bsa_name, block_set->bsas()) {
        ui->bsaComboBox->addItem(QString::fromStdString(bsa_name));
    }
    if (block_set->bsas().empty()) {
        ui->bsaWidget->hide();
    } else {
        ui->bsaWidget->show();
    }
    prev_row_ = -1;
    block_set_model_->update_filter();
}

BlockSetPtr BlockSetWidget::block_set() const {
    return block_set_model_->block_set();
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

void BlockSetWidget::moveBsaWidget(
    BlockSetWidget* dst,
    BlockSetWidget* src) {
    QWidget* bsaWidget = src->ui->bsaWidget;
    bsaWidget->setParent(dst);
    dst->ui->global_bsa_layout->addWidget(bsaWidget);
    connect(src, SIGNAL(blockClicked(QString)),
            dst, SLOT(onblockClicked(QString)));
}

void BlockSetWidget::onblockClicked(QString name) {
    ui->tabWidget->setCurrentWidget(ui->common_bsa);
    set_bsa(name.toStdString());
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
    //
    bsa_view_->selectionModel()->setCurrentIndex(QModelIndex(),
            QItemSelectionModel::Clear);
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
    if (!index.isValid()) {
        return;
    }
    Fragment* fragment = bsa_model_->index2fragment(index);
    if (fragment) {
        set_block(fragment->block());
        alignment_view_->select_fragment(fragment);
        bsa_view_->setCurrentIndex(index);
        std::string name = fragment->block()->name();
        emit blockClicked(QString::fromStdString(name));
    }
}

void BlockSetWidget::jump_to_f(Fragment* fragment, int col) {
    ASSERT_TRUE(fragment->block());
    set_block(fragment->block());
    int row = alignment_model_->fragment_index(fragment);
    QModelIndex index = alignment_model_->index(row, col);
    alignment_view_->selectionModel()->clearSelection();
    alignment_view_->setCurrentIndex(index);
    alignment_view_->scrollTo(index);
}

void BlockSetWidget::set_bsa(std::string bsa_name) {
    bsa_model_->set_bsa(bsa_name);
    QComboBox* cb = ui->bsaComboBox;
    int row = cb->findText(QString::fromStdString(bsa_name));
    cb->setCurrentIndex(row);
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
    block_set_model_->set_more_than_1(state == Qt::Checked);
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
    block_set_model_->set_pattern(pattern.toStdString());
}

void BlockSetWidget::on_clearBlockNameButton_clicked() {
    ui->blockNameLineEdit->setText("");
    on_blockNameLineEdit_editingFinished();
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

void BlockSetWidget::onSearchStarted() {
    ui->nonunique->setEnabled(false);
    ui->blockNameLineEdit->setEnabled(false);
    ui->clearBlockNameButton->setEnabled(false);
}

void BlockSetWidget::onSearchFinished() {
    ui->nonunique->setEnabled(true);
    ui->blockNameLineEdit->setEnabled(true);
    ui->clearBlockNameButton->setEnabled(true);
}

