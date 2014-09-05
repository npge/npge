/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef GUI_BLOCK_SEARCHER_HPP_
#define GUI_BLOCK_SEARCHER_HPP_

#include <boost/function.hpp>
#include <QtCore>

#include "global.hpp"
#include "gui-global.hpp"
#include "SortedVector.hpp"

using namespace npge;

struct BlockSearcher : public QObject, public QRunnable {
    Q_OBJECT

public:
    typedef boost::function<bool(const Block*)> BlockChecker;

    ConstBlocks* blocks_;
    SortedBlocks* filtered_blocks_;
    BlockChecker block_checker_;

    void run();

signals:
    void searchingFinished(QString message);
};

#endif

