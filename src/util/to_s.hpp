/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_TOS_HPP_
#define NPGE_TOS_HPP_

#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/cast.hpp>

#define L_CAST boost::lexical_cast
#define D_CAST boost::polymorphic_downcast

#define TO_S(x) L_CAST<std::string>(x)

#endif

