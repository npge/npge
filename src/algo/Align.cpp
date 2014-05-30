/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Align.hpp"
#include "MetaAligner.hpp"
#include "MoveGaps.hpp"
#include "CutGaps.hpp"
#include "Filter.hpp"
#include "SelfOverlapsResolver.hpp"

namespace npge {

class LiteAlignLoop : public Pipe {
public:
    LiteAlignLoop() {
        set_max_loops(-1);
        add(new MoveGaps);
        add(new CutGaps);
    }
};

LiteAlign::LiteAlign() {
    add(new MetaAligner);
    add(new LiteAlignLoop);
    declare_bs("target", "Aligned blockset");
}

const char* LiteAlign::name_impl() const {
    return "Align, move and cut gap";
}

class AlignLoop : public Pipe {
public:
    AlignLoop() {
        set_max_loops(-1);
        add(new MoveGaps);
        add(new CutGaps);
        add(new Filter);
    }
};

Align::Align() {
    add(new MetaAligner);
    add(new SelfOverlapsResolver);
    add(new MetaAligner);
    add(new AlignLoop);
    declare_bs("target", "Aligned blockset");
}

const char* Align::name_impl() const {
    return "LiteAlign + Filter + SelfOverlapsResolver";
}

}

