#ifndef ALIGNMENTMODEL_HPP
#define ALIGNMENTMODEL_HPP

#include <vector>
#include <QAbstractTableModel>

#include "global.hpp"

using namespace bloomrepeats;

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

signals:

public slots:
    void set_block(const Block* block);

    /** Move rows Up/Down.
    List of rows should be sorted. It is changed to
    new list of selected rows
    */
    void move_rows(std::vector<int>& rows, bool up);

private:
    const Block* block_;
    int length_;
    std::vector<bool> ident_;
    std::vector<bool> gap_;
    std::string consensus_;
    std::vector<const Fragment*> fragments_;
};

#endif // ALIGNMENTMODEL_HPP

