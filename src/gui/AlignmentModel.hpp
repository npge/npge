/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef ALIGNMENTMODEL_HPP
#define ALIGNMENTMODEL_HPP

#include <vector>
#include <set>
#include <QAbstractTableModel>

#include "global.hpp"
#include "FragmentCollection.hpp"

using namespace npge;

struct GeneInfo {
    bool is_gene : 1;
    bool is_reverse : 1;
    bool is_start : 1;
    bool is_stop : 1;
    bool gene_overlap : 1;
};

class AlignmentModel : public QAbstractTableModel {
    Q_OBJECT
public:
    explicit AlignmentModel(const Block* block = 0, QObject* parent = 0);

    QVariant data(const QModelIndex& index,
                  int role = Qt::DisplayRole) const;

    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;

    int rowCount(const QModelIndex& parent = QModelIndex()) const;

    int columnCount(const QModelIndex& parent = QModelIndex()) const;

    bool get_show_genes() const {
        return show_genes_;
    }

    bool has_genes() const {
        return has_genes_;
    }

    Fragment* fragment_at(int row) const;

    int fragment_index(Fragment* f) const;

    const std::vector<Fragment*>& fragments() const {
        return fragments_;
    }

    Fragment* test_genes(const QModelIndex& index,
                         GeneInfo* gene_info) const;

    void test_col(int col, bool& ident, bool& gap) const;

    char consensus_char(int col) const;

    bool test_gap(const QModelIndex& index) const;

    bool is_low_col(int col) const;

    BlockSetPtr block_set() const;

    const Block* block() const {
        return block_;
    }

signals:

public slots:
    void set_block_set(BlockSetPtr block_set);

    void set_block(const Block* block);

    /** Change order of fragments, lose genes */
    void set_fragments(const Fragments& ff);

    /** Move rows Up/Down.
    List of rows should be sorted. It is changed to
    new list of selected rows
    */
    void move_rows(std::vector<int>& rows, bool up);

    void add_genes(Fragment* fragment, const Fragments& genes);

    void set_split_parts(const Blocks& blocks);

    void set_low_similarity(const Blocks& blocks);

    void set_show_genes(bool show_genes);

    Fragment* logical_neighbor(Fragment* f, int ori) const;

    void set_genes_s2f(const VectorFc* genes_s2f);

private:
    std::vector<Fragment*> fragments_;
    std::vector<std::vector<Fragment*> > genes_;
    typedef std::map<Fragment*, int> Fragment2Int;
    mutable Fragment2Int split_parts_;
    typedef std::set<int> IntSet;
    IntSet low_similarity_;
    BlockSetPtr block_set_;
    VectorFc s2f_;
    const VectorFc* genes_s2f_;
    const Block* block_;
    int length_;
    bool has_genes_, show_genes_;

    bool is_gene_start_stop(Fragment* gene, int ori) const;
};

#endif // ALIGNMENTMODEL_HPP

