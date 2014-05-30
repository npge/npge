#ifndef ALIGNMENTVIEW_HPP
#define ALIGNMENTVIEW_HPP

#include <QTableView>

#include "global.hpp"
#include "gui-global.hpp"

using namespace npge;

class AlignmentView : public QTableView {
    Q_OBJECT
public:
    explicit AlignmentView(QWidget* parent = 0);

    void set_model(AlignmentModel* model);

    typedef std::map<Sequence*, Fragment*> Seq2Fragment;

    void set_first_last(const Seq2Fragment& first,
                        const Seq2Fragment& last);

    void keyPressEvent(QKeyEvent* event);

public slots:
    void select_fragment(Fragment* fragment);

signals:
    void fragment_selected(Fragment* fragment, int column);
    void jump_to(Fragment* to, int column);

private slots:
    void clicked_f(const QModelIndex& index);

private:
    Seq2Fragment seq2first_, seq2last_;
};

#endif // ALIGNMENTVIEW_HPP

