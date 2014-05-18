/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/current_function.hpp>

#include "to_s.hpp"
#include "config.hpp"

namespace bloomrepeats {
void assertion_failed_msg(char const* expr, char const* msg,
                          char const* function, char const* file, long line);
}

#ifdef BR_ASSERTS

#define ASSERT_MSG(expr, msg) ((expr) \
    ? ((void)0) \
    : ::bloomrepeats::assertion_failed_msg(#expr, msg,\
        BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))

#else

#define ASSERT_MSG(expr, msg) ((void)0)

#endif

#define ASSERT_TRUE(exp_a) \
    ASSERT_MSG(exp_a, TO_S(exp_a).c_str())

#define ASSERT_FALSE(exp_a) \
    ASSERT_MSG(!(exp_a), ("!" + TO_S(exp_a)).c_str())

#define ASSERT_EQ(exp_a, exp_b) \
    ASSERT_MSG((exp_a) == (exp_b), \
        (TO_S(exp_a) + " == " + TO_S(exp_b)).c_str())

#define ASSERT_NE(exp_a, exp_b) \
    ASSERT_MSG((exp_a) != (exp_b), \
        (TO_S(exp_a) + " != " + TO_S(exp_b)).c_str())

#define ASSERT_LT(exp_a, exp_b) \
    ASSERT_MSG((exp_a) < (exp_b), \
        (TO_S(exp_a) + " < " + TO_S(exp_b)).c_str())

#define ASSERT_LTE(exp_a, exp_b) \
    ASSERT_MSG((exp_a) <= (exp_b), \
        (TO_S(exp_a) + " <= " + TO_S(exp_b)).c_str())

#define ASSERT_GT(exp_a, exp_b) \
    ASSERT_MSG((exp_a) > (exp_b), \
        (TO_S(exp_a) + " > " + TO_S(exp_b)).c_str())

#define ASSERT_GTE(exp_a, exp_b) \
    ASSERT_MSG((exp_a) >= (exp_b), \
        (TO_S(exp_a) + " >= " + TO_S(exp_b)).c_str())

