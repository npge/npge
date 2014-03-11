#include <algorithm>
#include <boost/foreach.hpp>
#include <QtGui>

#include "AlignmentModel.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "char_to_size.hpp"
#include "block_stat.hpp"
#include "throw_assert.hpp"
#include "convert_position.hpp"

AlignmentModel::AlignmentModel(const Block* block, QObject* parent) :
    QAbstractTableModel(parent),
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
    bool is_gene, is_reverse, is_start;
    if (role == Qt::TextAlignmentRole) {
        return Qt::AlignCenter;
    }
    if (role == Qt::FontRole) {
        test_genes(index, is_gene, is_reverse, is_start);
        if (is_reverse) {
            QFont font;
            font.setUnderline(true);
            return font;
        }
    }
    if (role == Qt::DisplayRole) {
        const Fragment* f = fragments_[index.row()];
        return QChar(f->alignment_at(index.column()) ? : '-');
    } else if (role == Qt::BackgroundRole) {
        test_genes(index, is_gene, is_reverse, is_start);
        if (is_start) {
            return Qt::black;
        }
        const Fragment* f = fragments_[index.row()];
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
        test_genes(index, is_gene, is_reverse, is_start);
        if (is_gene) {
            return Qt::white;
        }
    } else if (role == Qt::ToolTipRole) {
        bool _;
        const Fragment* gene = test_genes(index, _, _, _);
        if (gene && gene->block()) {
            return QString("%1, %2 bp")
                   .arg(QString::fromStdString(gene->block()->name()))
                   .arg(gene->length());
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
            if (section < fragments_.size()) {
                return QString::fromStdString(fragments_[section]->id());
            } else if (section == fragments_.size()) {
                return tr("consensus");
            }
        }
    } else if (orientation == Qt::Horizontal) {
        if (role == Qt::DisplayRole) {
            int digit = column_digit(section, length_);
            QString h = (digit == -1) ? "" : QString::number(digit);
            h += "\n";
            h += consensus_[section];
            return h;
        } else if (role == Qt::BackgroundRole) {
            return (!ident_[section]) ? Qt::white : (gap_[section]) ?
                   Qt::gray : Qt::black;
        } else if (role == Qt::ForegroundRole) {
            return (ident_[section] && !gap_[section]) ? Qt::white : Qt::black;
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

const Fragment* AlignmentModel::fragment_at(int row) const {
    return fragments_[row];
}

int AlignmentModel::fragment_index(const Fragment* f) const {
    for (int i = 0; i < fragments_.size(); i++) {
        if (fragments_[i] == f) {
            return i;
        }
    }
    return -1;
}

bool AlignmentModel::test_gap(const QModelIndex& index) const {
    const Fragment* f = fragments_[index.row()];
    return f->alignment_at(index.column()) == 0;
}

void AlignmentModel::set_block(const Block* block) {
    beginResetModel();
    block_ = block;
    if (block_) {
        length_ = block_->alignment_length();
        std::vector<const Fragment*> fragments(block_->begin(), block_->end());
        fragments_.swap(fragments);
        ident_.resize(length_);
        gap_.resize(length_);
        for (int col = 0; col <= length_; col++) {
            bool ident, gap, pure_gap;
            int atgc[LETTERS_NUMBER];
            test_column(block_, col, ident, gap, pure_gap, atgc);
            ident_[col] = ident;
            gap_[col] = gap;
        }
        consensus_ = block_->consensus_string();
    } else {
        length_ = 0;
        fragments_.clear();
        ident_.clear();
        gap_.clear();
        consensus_.clear();
    }
    genes_.clear();
    genes_.resize(fragments_.size());
    has_genes_ = false;
    endResetModel();
}

void AlignmentModel::set_fragments(const std::vector<const Fragment*>&
        ff) {
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
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        return *f1 < *f2;
    }
} fragment_compare_g;

void AlignmentModel::add_genes(const Fragment* fragment,
                               const std::vector<Fragment*>& genes) {
    beginResetModel();
    int fragment_id = -1;
    for (int i = 0; i < fragments_.size(); i++) {
        if (fragments_[i] == fragment) {
            fragment_id = i;
            break;
        }
    }
    BOOST_ASSERT(fragment_id != -1);
    std::vector<Fragment*>& g = genes_[fragment_id];
    BOOST_FOREACH (Fragment* gene, genes) {
        g.push_back(gene);
        has_genes_ = true;
    }
    std::sort(g.begin(), g.end(), fragment_compare_g);
    endResetModel();
}

void AlignmentModel::set_show_genes(bool show_genes) {
    beginResetModel();
    show_genes_ = show_genes;
    endResetModel();
}

const Fragment* AlignmentModel::test_genes(
    const QModelIndex& index,
    bool& is_gene,
    bool& is_reverse,
    bool& is_start) const {
    is_gene = false;
    is_reverse = false;
    is_start = false;
    if (!has_genes_ || !show_genes_) {
        return 0;
    }
    const Fragment* f = fragments_[index.row()];
    const AlignmentRow* row = f->row();
    int f_pos;
    if (row) {
        f_pos = row->map_to_fragment(index.column());
    } else if (index.column() < f->length()) {
        f_pos = index.column();
    } else {
        return 0;
    }
    if (f_pos == -1) {
        return 0;
    }
    int s_pos = frag_to_seq(f, f_pos);
    BOOST_FOREACH (const Fragment* gene, genes_[index.row()]) {
        if (gene->has(s_pos)) {
            is_gene = true;
            if (gene->ori() != f->ori()) {
                is_reverse = true;
            }
            int g_pos = seq_to_frag(gene, s_pos);
            if (g_pos < 3) {
                is_start = true;
            }
            return gene;
        }
    }
    return 0;
}

