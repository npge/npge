/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef GUI_BLOCK_SEARCHER_HPP_
#define GUI_BLOCK_SEARCHER_HPP_

#include <QtCore>

#ifndef Q_MOC_RUN
#include "global.hpp"
#include "gui-global.hpp"
#endif

using namespace npge;

struct BlockSearcher : public QObject, public QRunnable {
    Q_OBJECT

public:
    Meta* meta_;
    BlockSetModel* model_;
    ConstBlocks* blocks_;
    SortedBlocks* filtered_blocks_;

    void run();

signals:
    void searchingFinished(QString message);
};

#endif

