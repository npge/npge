/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "OldMakePangenome.hpp"
#include "RemoveNames.hpp"
#include "UniqueNames.hpp"
#include "Filter.hpp"
#include "Clear.hpp"
#include "Union.hpp"
#include "Joiner.hpp"
#include "Align.hpp"
#include "OverlaplessUnion.hpp"
#include "Rest.hpp"
#include "AddBlastBlocks.hpp"
#include "OneByOne.hpp"
#include "Info.hpp"
#include "Filter.hpp"

namespace npge {

OldMakePangenome::OldMakePangenome() {
    set_max_loops(-1);
    add(new Filter);
    add(new Clear, "target=backup");
    add(new Union, "target=backup other=target");
    add(new Joiner);
    add(new Filter);
    add(new Align);
    add(new Filter);
    add(new OverlaplessUnion, "other=backup");
    //
    add(new Clear, "target=hits");
    add(new Rest, "other=target");
    add(new AddBlastBlocks, "target=hits other=target");
    add(new OneByOne, "other=hits");
    add(new Filter);
    add(new Align);
    add(new Filter);
    add(new Rest, "other=target");
    //
    add(new RemoveNames, "--remove-seqs-names:=0");
    add(new UniqueNames);
    add(new Info);
    declare_bs("target", "Target blockset with pre-pangenome, "
               "which is evolved to pangenome");
}

const char* OldMakePangenome::name_impl() const {
    return "Run blast and Joiner until this blockset "
           "becomes pangenome (deprecated)";
}

}

