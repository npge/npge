/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SizeLimits.hpp"
#include "Processor.hpp"
#include "Decimal.hpp"

namespace npge {

void add_lite_size_limits_options(Processor* p) {
    p->add_gopt("min-fragment", "Minimum fragment length",
                "MIN_LENGTH");
    p->add_opt("min-block", "Minimum block size", 2);
}

void add_size_limits_options(Processor* p) {
    add_lite_size_limits_options(p);
    p->add_opt("max-fragment", "Maximum fragment length (-1 = all)",
               -1);
    p->add_opt("max-block", "Maximum block size (-1 = all)", -1);
    p->add_opt("min-spreading",
               "Minimum fragment length spreading ((max - min) / avg)", D(0.0));
    p->add_gopt("max-spreading", "Maximum fragment length spreading",
                "MAX_SPREADING");
    p->add_gopt("min-identity",
                "Minimum block identity (only if alignment is known, "
                "columns without gaps as 1, columns with gaps as 0.5)",
                "MIN_IDENTITY");
    p->add_opt("max-identity", "Maximum block identity", D(1.0));
    p->add_opt("min-gaps",
               "Min gap columns percentage (only if alignment is known)",
               D(0.0));
    p->add_gopt("max-gaps", "Max gap columns percentage", "MAX_GAPS");
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
    p->set_opt_value("min-spreading", D(0.0));
    p->set_opt_value("max-spreading", D(999.9));
    p->set_opt_value("min-identity", D(0.0));
    p->set_opt_value("max-identity", D(1.0));
    p->set_opt_value("min-gaps", D(0.0));
    p->set_opt_value("max-gaps", D(1.0));
}

}

