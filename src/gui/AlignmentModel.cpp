/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <QtGui>

#include "AlignmentModel.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "AlignmentRow.hpp"
#include "char_to_size.hpp"
#include "block_stat.hpp"
#include "throw_assert.hpp"
#include "convert_position.hpp"

AlignmentModel::AlignmentModel(const Block* block, QObject* parent) :
    QAbstractTableModel(parent),
    genes_s2f_(0),
    has_genes_(false), show_genes_(true) {
    set_block(block);
}

QColor colors_[5] = {
    QRgb(0xFF64F73F), // green
    QRgb(0xFF3C88EE), // blue
    QRgb(0xFFEB413C), // red
    QRgb(0xFFFFB340), // orange
    Qt::red
};

QVariant AlignmentModel::data(const QModelIndex& index, int role) const {
    if (!index.isValid()) {
        return QVariant();
    }
    if (role == Qt::TextAlignmentRole) {
        return Qt::AlignCenter;
    }
    if (role == Qt::FontRole) {
        GeneInfo go;
        test_genes(index, &go);
        if (go.is_reverse) {
            QFont font;
            font.setUnderline(true);
            return font;
        }
    }
    if (role == Qt::DisplayRole) {
        Fragment* f = fragments_[index.row()];
        return QChar(f->alignment_at(index.column()) ? : '-');
    } else if (role == Qt::BackgroundRole) {
        GeneInfo go;
        test_genes(index, &go);
        if (go.is_start) {
            return Qt::black;
        } else if (go.is_stop) {
            return Qt::gray;
        } else if (go.gene_overlap) {
            return Qt::magenta;
        }
        Fragment* f = fragments_[index.row()];
        char c = f->alignment_at(index.column());
        if (c == 0) {
            // gap
            return Qt::white;
        }
        size_t s = char_to_size(c);
        if (s < 5) {
            return colors_[s];
        }
    } else if (role == Qt::ForegroundRole) {
        GeneInfo go;
        test_genes(index, &go);
        if (go.is_gene) {
            return Qt::white;
        }
    } else if (role == Qt::ToolTipRole) {
        GeneInfo go;
        test_genes(index, &go);
        QStringList gene_texts;
        BOOST_FOREACH (Fragment* gene, go.genes) {
            if (gene && gene->block()) {
                QString gene_text = QString("%1, %2 bp %3")
                       .arg(QString::fromStdString(gene->block()->name()))
                       .arg(gene->length())
                       .arg(go.is_reverse ? "<" : ">");
                gene_texts << gene_text;
            }
        }
        if (!go.genes.empty()) {
            return gene_texts.join(" %< ");
        }
    }
    return QVariant();
}

static int column_digit(int section, int length) {
    int col = section + 1;
    int next_10 = ((col + 9) / 10) * 10;
    if (next_10 >= length) {
        return -1;
    }
    int digits_shift = next_10 - col;
    if (digits_shift == 0) {
        return 0;
    }
    for (int i = 0; i < digits_shift; i++) {
        next_10 /= 10;
    }
    if (next_10 == 0) {
        return -1;
    }
    return next_10 % 10;
}

QVariant AlignmentModel::headerData(int section, Qt::Orientation orientation,
                                    int role) const {
    if (orientation == Qt::Vertical) {
        if (role == Qt::DisplayRole) {
            QString result;
            result = QString::fromStdString(fragments_[section]->id());
            Fragment* f = fragments_[section];
            int part = split_parts_[f];
            if (part) {
                result = QString::number(part) + " " + result;
            }
            return result;
        }
    } else if (orientation == Qt::Horizontal) {
        if (role == Qt::DisplayRole) {
            int digit = column_digit(section, length_);
            QString h = (digit == -1) ? "" : QString::number(digit);
            h += "\n";
            h += consensus_char(section);
            return h;
        } else if (role == Qt::BackgroundRole) {
            bool ident, gap;
            test_col(section, ident, gap);
            return (!ident) ? Qt::white : gap ?
                   Qt::gray : Qt::black;
        } else if (role == Qt::ForegroundRole) {
            bool ident, gap;
            test_col(section, ident, gap);
            return ident && !gap ? Qt::white : Qt::black;
        } else if (role == Qt::UserRole) {
            if (is_low_col(section)) {
                return true;
            }
        }
    }
    return QAbstractTableModel::headerData(section, orientation, role);
}

int AlignmentModel::rowCount(const QModelIndex&) const {
    return fragments_.size();
}

int AlignmentModel::columnCount(const QModelIndex&) const {
    return length_;
}

Fragment* AlignmentModel::fragment_at(int row) const {
    return fragments_[row];
}

int AlignmentModel::fragment_index(Fragment* f) const {
    for (int i = 0; i < fragments_.size(); i++) {
        if (fragments_[i] == f) {
            return i;
        }
    }
    return -1;
}

bool AlignmentModel::test_gap(const QModelIndex& index) const {
    Fragment* f = fragments_[index.row()];
    return f->alignment_at(index.column()) == 0;
}

bool AlignmentModel::is_low_col(int col) const {
    return low_similarity_.find(col) != low_similarity_.end();
}

BlockSetPtr AlignmentModel::block_set() const {
    return block_set_;
}

typedef std::map<Fragment*, int> Fragment2Int;

struct SeqComp {
    SeqComp(Fragment2Int& split_parts):
        split_parts_(split_parts) {
    }

    bool operator()(Fragment* f1, Fragment* f2) const {
        typedef boost::tuple < int, const std::string&,
                Fragment& > Tie;
        return Tie(split_parts_[f1], f1->seq()->name(), *f1) <
               Tie(split_parts_[f2], f2->seq()->name(), *f2);
    }

private:
    mutable Fragment2Int split_parts_;
};

struct ByPosInBlockCmp {
    ByPosInBlockCmp(const Block* block):
        block_(block) {
    }

    bool operator()(Fragment* a, Fragment* b) const {
        int block_length = block_->alignment_length();
        int a1 = block_pos(a, 0, block_length);
        int a2 = block_pos(a, a->length() - 1, block_length);
        int b1 = block_pos(b, 0, block_length);
        int b2 = block_pos(b, b->length() - 1, block_length);
        int a_min = std::min(a1, a2);
        int a_max = std::max(a1, a2);
        int b_min = std::min(b1, b2);
        int b_max = std::max(b1, b2);
        typedef boost::tuple<int, int> Tie;
        return Tie(a_min, a_max) < Tie(b_min, b_max);
    }

private:
    const Block* block_;
};

void AlignmentModel::set_block_set(BlockSetPtr block_set) {
    block_set_ = block_set;
    s2f_.clear();
    s2f_.add_bs(*block_set_);
    s2f_.prepare();
}

void AlignmentModel::set_block(const Block* block) {
    beginResetModel();
    block_ = block;
    if (block_) {
        length_ = block_->alignment_length();
        std::vector<Fragment*> fragments(block_->begin(),
                                         block_->end());
        std::sort(fragments.begin(), fragments.end(),
                  SeqComp(split_parts_));
        fragments_.swap(fragments);
    } else {
        length_ = 0;
        fragments_.clear();
    }
    genes_.clear();
    genes_.resize(fragments_.size());
    has_genes_ = false;
    endResetModel();
}

void AlignmentModel::set_fragments(const Fragments& ff) {
    beginResetModel();
    fragments_ = ff;
    genes_.clear();
    genes_.resize(fragments_.size());
    has_genes_ = false;
    endResetModel();
}

void AlignmentModel::move_rows(std::vector<int>& rows, bool up) {
    beginResetModel();
    if (up) {
        int prev_row = -100;
        BOOST_FOREACH (int& row, rows) {
            if (row > 0 && row - 1 != prev_row) {
                std::swap(fragments_[row], fragments_[row - 1]);
                genes_[row].swap(genes_[row - 1]);
                row -= 1;
            }
            prev_row = row;
        }
    } else {
        int prev_row = -100;
        BOOST_REVERSE_FOREACH (int& row, rows) {
            if (row < fragments_.size() - 1 && row + 1 != prev_row) {
                std::swap(fragments_[row], fragments_[row + 1]);
                genes_[row].swap(genes_[row + 1]);
                row += 1;
            }
            prev_row = row;
        }
    }
    endResetModel();
}

static struct FragmentCompareG {
    bool operator()(Fragment* f1, Fragment* f2) const {
        return *f1 < *f2;
    }
} fragment_compare_g;

void AlignmentModel::add_genes(Fragment* fragment,
                               const Fragments& genes) {
    beginResetModel();
    int fragment_id = -1;
    for (int i = 0; i < fragments_.size(); i++) {
        if (fragments_[i] == fragment) {
            fragment_id = i;
            break;
        }
    }
    ASSERT_NE(fragment_id, -1);
    std::vector<Fragment*>& g = genes_[fragment_id];
    BOOST_FOREACH (Fragment* gene, genes) {
        g.push_back(gene);
        has_genes_ = true;
    }
    std::sort(g.begin(), g.end(), fragment_compare_g);
    if (fragment->ori() == -1) {
        std::reverse(g.begin(), g.end());
    }
    endResetModel();
}

void AlignmentModel::set_split_parts(const Blocks& blocks) {
    split_parts_.clear();
    int part = 0;
    BOOST_FOREACH (Block* block, blocks) {
        part += 1;
        BOOST_FOREACH (Fragment* f, *block) {
            Fragment* orig_f = 0;
            BOOST_FOREACH (Fragment* f1, *block_) {
                if (*f1 == *f) {
                    orig_f = f1;
                }
            }
            ASSERT_TRUE(orig_f);
            split_parts_[orig_f] = part;
        }
    }
    std::sort(fragments_.begin(), fragments_.end(),
              SeqComp(split_parts_));
}

void AlignmentModel::set_low_similarity(const Blocks& blocks) {
    low_similarity_.clear();
    BOOST_FOREACH (Block* block, blocks) {
        int min_col, max_col;
        find_slice(min_col, max_col, block_, block);
        for (int col = min_col; col <= max_col; col++) {
            low_similarity_.insert(col);
        }
    }
}

void AlignmentModel::set_show_genes(bool show_genes) {
    beginResetModel();
    show_genes_ = show_genes;
    endResetModel();
}

Fragment* AlignmentModel::logical_neighbor(Fragment* f,
        int ori) const {
    return s2f_.logical_neighbor(f, ori);
}

void AlignmentModel::set_genes_s2f(const VectorFc* genes_s2f) {
    genes_s2f_ = genes_s2f;
}

// ori = -1 for start codon, 1 for stop codon
bool AlignmentModel::is_gene_start_stop(
    Fragment* gene, int ori) const {
    Block* gene_block = gene->block();
    using namespace boost::algorithm;
    if (!starts_with(gene_block->name(), "CDS")) {
        return false;
    }
    if (gene_block->size() == 1) {
        return true;
    }
    ASSERT_TRUE(genes_s2f_);
    Fragment* neighbour =
        genes_s2f_->logical_neighbor(gene, ori);
    if (neighbour && neighbour->block() == gene_block) {
        return false;
    }
    return true;
}

void AlignmentModel::test_genes(const QModelIndex& index,
                                GeneInfo* gene_info) const {
    gene_info->is_gene = false;
    gene_info->is_reverse = false;
    gene_info->is_start = false;
    gene_info->is_stop = false;
    gene_info->gene_overlap = false;
    if (!has_genes_ || !show_genes_) {
        return;
    }
    Fragment* f = fragments_[index.row()];
    const AlignmentRow* row = f->row();
    int f_pos;
    if (row) {
        f_pos = row->map_to_fragment(index.column());
    } else if (index.column() < f->length()) {
        f_pos = index.column();
    } else {
        return;
    }
    if (f_pos == -1) {
        return;
    }
    int s_pos = frag_to_seq(f, f_pos);
    BOOST_FOREACH (Fragment* gene, genes_[index.row()]) {
        if (gene->has(s_pos)) {
            gene_info->is_gene = true;
            gene_info->is_reverse = (gene->ori() != f->ori());
            int g_pos = seq_to_frag(gene, s_pos);
            if (g_pos < 3 && is_gene_start_stop(gene, -1)) {
                gene_info->is_start = true;
            }
            if (g_pos >= gene->length() - 3 &&
                    is_gene_start_stop(gene, 1)) {
                gene_info->is_stop = true;
            }
            gene_info->genes.push_back(gene);
        }
    }
    Fragments& genes = gene_info->genes;
    std::sort(genes.begin(), genes.end(),
              fragment_compare_g);
    if (f->ori() == -1) {
        std::reverse(genes.begin(), genes.end());
    }
    if (genes.size() >= 2) {
        gene_info->gene_overlap = true;
    }
}

void AlignmentModel::test_col(int col,
                              bool& ident, bool& gap) const {
    bool pure_gap;
    int atgc[LETTERS_NUMBER];
    test_column(block_, col, ident, gap, pure_gap, atgc);
}

char AlignmentModel::consensus_char(int col) const {
    return block_->consensus_char(col);
}

