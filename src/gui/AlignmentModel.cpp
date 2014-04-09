#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
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
        const Fragment* f = fragments_[index.row()];
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
        GeneInfo go;
        test_genes(index, &go);
        if (go.is_gene) {
            return Qt::white;
        }
    } else if (role == Qt::ToolTipRole) {
        GeneInfo go;
        const Fragment* gene = test_genes(index, &go);
        if (gene && gene->block()) {
            return QString("%1, %2 bp %3")
                   .arg(QString::fromStdString(gene->block()->name()))
                   .arg(gene->length())
                   .arg(go.is_reverse ? "<" : ">");
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
            const Fragment* f = fragments_[section];
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
            if (low_similarity_.find(section) !=
                    low_similarity_.end()) {
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

typedef std::map<const Fragment*, int> Fragment2Int;

struct SeqComp {
    SeqComp(Fragment2Int& split_parts):
        split_parts_(split_parts)
    { }

    bool operator()(const Fragment* f1, const Fragment* f2) const {
        typedef boost::tuple<int, const std::string&,
                             const Fragment&> Tie;
        return Tie(split_parts_[f1], f1->seq()->name(), *f1) <
               Tie(split_parts_[f2], f2->seq()->name(), *f2);
    }

private:
    mutable Fragment2Int split_parts_;
};

void AlignmentModel::set_block(const Block* block) {
    beginResetModel();
    block_ = block;
    if (block_) {
        length_ = block_->alignment_length();
        std::vector<const Fragment*> fragments(block_->begin(),
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
                               const Fragments& genes) {
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
            BOOST_ASSERT(orig_f);
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

const Fragment* AlignmentModel::test_genes(const QModelIndex& index,
        GeneInfo* gene_info) const {
    gene_info->is_gene = false;
    gene_info->is_reverse = false;
    gene_info->is_start = false;
    gene_info->is_stop = false;
    gene_info->gene_overlap = false;
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
    const Fragment* result = 0;
    BOOST_FOREACH (const Fragment* gene, genes_[index.row()]) {
        if (gene->has(s_pos)) {
            if (result) {
                gene_info->gene_overlap = true;
            }
            gene_info->is_gene = true;
            gene_info->is_reverse = (gene->ori() != f->ori());
            int g_pos = seq_to_frag(gene, s_pos);
            if (g_pos < 3) {
                gene_info->is_start = true;
            }
            if (g_pos >= gene->length() - 3) {
                gene_info->is_stop = true;
            }
            result = gene;
        } else if (result) {
            break;
        }
    }
    return result;
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

