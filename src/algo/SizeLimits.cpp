/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SizeLimits.hpp"
#include "Processor.hpp"

namespace bloomrepeats {

void add_size_limits_options(Processor* p) {
    p->add_opt("min-fragment", "Minimal fragment length", 100);
    p->add_opt("max-fragment", "Maximal fragment length (-1 = all)", -1);
    p->add_opt("min-block", "Minimal block size", 2);
    p->add_opt("max-block", "Maximal block size (-1 = all)", -1);
    p->add_opt("min-spreading",
               "Minimal fragment length spreading ((max - min) / avg)", 0.0);
    p->add_opt("max-spreading", "Maximal fragment length spreading", 0.2);
    p->add_opt("min-identity",
               "Minimal block identity (only if alignment is known, "
               "columns without gaps as 1, columns with gaps as 0.5)", 0.9);
    p->add_opt("max-identity", "Maximal block identity", 0.9);
    p->add_opt("min-gaps",
               "Min gap columns percentage (only if alignment is known)", 0.0);
    p->add_opt("max-gaps", "Max gap columns percentage", 0.2);
    //
    p->add_opt_rule("min-fragment >= 0");
    p->add_opt_rule("max-fragment >= -1");
    p->add_opt_rule("min-block >= 0");
    p->add_opt_rule("max-block >= -1");
    p->add_opt_rule("min-spreading >= 0.0");
    p->add_opt_rule("max-spreading >= 0.0");
    p->add_opt_rule("min-identity >= 0.0");
    p->add_opt_rule("min-identity <= 1.0");
    p->add_opt_rule("max-identity >= 0.0");
    p->add_opt_rule("max-identity <= 1.0");
    p->add_opt_rule("min-gaps >= 0.0");
    p->add_opt_rule("min-gaps <= 1.0");
    p->add_opt_rule("max-gaps >= 0.0");
    p->add_opt_rule("max-gaps <= 1.0");
}

void allow_everything(Processor* p) {
    p->set_opt_value("min-fragment", 0);
    p->set_opt_value("max-fragment", -1);
    p->set_opt_value("min-block", 0);
    p->set_opt_value("max-block", -1);
    p->set_opt_value("min-spreading", 0.0);
    p->set_opt_value("max-spreading", 999999.9);
    p->set_opt_value("min-identity", 0.0);
    p->set_opt_value("max-identity", 1.0);
    p->set_opt_value("min-gaps", 0.0);
    p->set_opt_value("max-gaps", 1.0);
}

}

