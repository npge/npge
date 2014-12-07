/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef GUI_BLOCK_SET_MODEL_HPP_
#define GUI_BLOCK_SET_MODEL_HPP_

#include <QtCore>

#ifndef Q_MOC_RUN
#include "global.hpp"
#include "gui-global.hpp"
#endif

using namespace npge;

enum {
    FRAGMENTS_C, COLUMNS_C,
    IDENTITY_C, GC_C,
    GENES_C, SPLIT_C, LOW_C
};

class BlockSetModel : public QAbstractTableModel {
    Q_OBJECT

public:
    explicit BlockSetModel(QObject* parent = 0);

    BlockSetPtr block_set() const;

    QVariant data(const QModelIndex& index,
                  int role = Qt::DisplayRole) const;

    QVariant headerData(int section,
                        Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;

    int rowCount(const QModelIndex& parent =
                     QModelIndex()) const;

    int columnCount(
        const QModelIndex& parent = QModelIndex()) const;

    const Block* block_at(int row) const;

    int block_index(const Block* block) const;

    const QPoint& xy_of(int row) const;

    void set_xy_of(int row, const QPoint& xy);

public slots:
    void set_block_set(BlockSetPtr block_set);

    void set_genes(BlockSetPtr genes);

    void find_genes(Fragments& overlap_genes,
                    Fragment* f) const;

    void find_genes(Fragments& overlap_genes,
                    const Block* block) const;

    Fragments return_genes(const Block* block) const;

    void construct_hits();

    bool check_block(const Block* block) const;

    void set_split_parts(BlockSetPtr split_parts);

    void find_split_parts(Fragments& ff, Fragment* f) const;

    void find_split_parts(Blocks& bb,
                          const Block* block) const;

    void set_low_similarity(BlockSetPtr low_similarity);

    void find_low_similarity(Fragments& ff, Fragment* f) const;

    void find_low_similarity(Blocks& bb,
                             const Block* block) const;

    void set_more_than_1(bool more_than_1);

    void set_pattern(const std::string& pattern);

    void update_filter();

    void onSearchingFinished(QString message);

    void onExceptionThrown(QString message);

signals:
    void exceptionThrown(QString message);

    void searchStarted();
    void searchFinished();

private:
    BlockSetPtr block_set_;
    std::vector<const Block*> blocks_;
    mutable std::vector<AlignmentStat*> stats_;
    mutable std::vector<QPoint> alignment_xy_;
    QStringList columns_;
    BlockSetPtr genes_;
    BlockSetPtr split_parts_;
    BlockSetPtr low_similarity_;
    VectorFc genes_s2f_;
    VectorFc split_s2f_;
    VectorFc low_s2f_;
    // filter
    // violation of MVC :(
    bool more_than_1_;
    std::string pattern_;
    SortedBlocks filtered_blocks_;
    BlockSetPtr hits_bs_;
    VectorFc hits_s2f_;
};

#endif

