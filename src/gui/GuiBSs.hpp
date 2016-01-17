/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef GUI_BSS_HPP_
#define GUI_BSS_HPP_

#ifndef Q_MOC_RUN
#include "global.hpp"
#endif

using namespace npge;

struct GuiBSs {
    BlockSetPtr pangenome_bs_;
    BlockSetPtr genes_bs_;
    BlockSetPtr split_parts_;
    BlockSetPtr low_similarity_;
    BlockSetPtr global_blocks_;
};

#endif

