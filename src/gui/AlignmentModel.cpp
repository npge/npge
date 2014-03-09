#include <algorithm>
#include <boost/foreach.hpp>
#include <QtGui>

#include "AlignmentModel.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "char_to_size.hpp"
#include "block_stat.hpp"

AlignmentModel::AlignmentModel(const Block* block, QObject* parent) :
    QAbstractTableModel(parent) {
    set_block(block);
}

QColor colors_[4] = {
    QRgb(0xFF64F73F), // green
    QRgb(0xFF3C88EE), // blue
    QRgb(0xFFEB413C), // red
    QRgb(0xFFFFB340) // orange
};

QVariant AlignmentModel::data(const QModelIndex& index, int role) const {
    if (role == Qt::TextAlignmentRole) {
        return Qt::AlignCenter;
    }
    if (index.row() < fragments_.size()) {
        // fragments
        if (role == Qt::DisplayRole) {
            const Fragment* f = fragments_[index.row()];
            return QChar(f->alignment_at(index.column()) ? : '-');
        } else if (role == Qt::BackgroundRole) {
            const Fragment* f = fragments_[index.row()];
            char c = f->alignment_at(index.column());
            size_t s = char_to_size(c);
            if (s < 4) {
                return colors_[s];
            }
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
                return "consensus";
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
    endResetModel();
}

void AlignmentModel::move_rows(std::vector<int>& rows, bool up) {
    beginResetModel();
    if (up) {
        int prev_row = -100;
        BOOST_FOREACH (int& row, rows) {
            if (row > 0 && row - 1 != prev_row) {
                std::swap(fragments_[row], fragments_[row - 1]);
                row -= 1;
            }
            prev_row = row;
        }
    } else {
        int prev_row = -100;
        BOOST_REVERSE_FOREACH (int& row, rows) {
            if (row < fragments_.size() - 1 && row + 1 != prev_row) {
                std::swap(fragments_[row], fragments_[row + 1]);
                row += 1;
            }
            prev_row = row;
        }
    }
    endResetModel();
}

