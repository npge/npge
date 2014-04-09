#ifndef ALIGNMENTMODEL_HPP
#define ALIGNMENTMODEL_HPP

#include <vector>
#include <set>
#include <QAbstractTableModel>

#include "global.hpp"

using namespace bloomrepeats;

struct GeneInfo {
    bool is_gene;
    bool is_reverse;
    bool is_start;
    bool is_stop;
    bool gene_overlap;
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

    const Fragment* fragment_at(int row) const;

    int fragment_index(const Fragment* f) const;

    const std::vector<const Fragment*>& fragments() const {
        return fragments_;
    }

    const Fragment* test_genes(const QModelIndex& index,
                               GeneInfo* gene_info) const;

    void test_col(int col, bool& ident, bool& gap) const;

    char consensus_char(int col) const;

    bool test_gap(const QModelIndex& index) const;

    const Block* block() const {
        return block_;
    }

signals:

public slots:
    void set_block(const Block* block);

    /** Change order of fragments, lose genes */
    void set_fragments(const std::vector<const Fragment*>& ff);

    /** Move rows Up/Down.
    List of rows should be sorted. It is changed to
    new list of selected rows
    */
    void move_rows(std::vector<int>& rows, bool up);

    void add_genes(const Fragment* fragment, const Fragments& genes);

    void set_split_parts(const Blocks& blocks);

    void set_low_similarity(const Blocks& blocks);

    void set_show_genes(bool show_genes);

private:
    std::vector<const Fragment*> fragments_;
    std::vector<std::vector<Fragment*> > genes_;
    typedef std::map<const Fragment*, int> Fragment2Int;
    mutable Fragment2Int split_parts_;
    typedef std::set<int> IntSet;
    IntSet low_similarity_;
    const Block* block_;
    int length_;
    bool has_genes_, show_genes_;
};

#endif // ALIGNMENTMODEL_HPP

