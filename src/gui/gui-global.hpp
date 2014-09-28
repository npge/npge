/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GUI_GLOBAL_HPP_
#define NPGE_GUI_GLOBAL_HPP_

#include <map>

#include "global.hpp"
#ifndef Q_MOC_RUN
#include "SortedVector.hpp"
#include "FragmentCollection.hpp"
#endif

using namespace npge;

class AlignmentView;
class AlignmentModel;
class BSAModel;
class BlockSetWidget;
class BlockSetModel;
class BlockSearcher;
class ReadingThread;

typedef std::vector<const Block*> ConstBlocks;
typedef SortedVector<const Block*> SortedBlocks;

typedef std::map<Sequence*, Fragment*> Seq2Fragment;

#endif

