/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#define BOOST_ENABLE_ASSERT_HANDLER

#include "boost-assert.hpp"
#include "to_s.hpp"

#define ASSERT_TRUE(exp_a) \
    BOOST_ASSERT_MSG(exp_a, TO_S(exp_a).c_str())

#define ASSERT_FALSE(exp_a) \
    BOOST_ASSERT_MSG(!(exp_a), ("!" + TO_S(exp_a)).c_str())

#define ASSERT_EQ(exp_a, exp_b) \
    BOOST_ASSERT_MSG((exp_a) == (exp_b), \
        (TO_S(exp_a) + " == " + TO_S(exp_b)).c_str())

#define ASSERT_NE(exp_a, exp_b) \
    BOOST_ASSERT_MSG((exp_a) != (exp_b), \
        (TO_S(exp_a) + " != " + TO_S(exp_b)).c_str())

#define ASSERT_LT(exp_a, exp_b) \
    BOOST_ASSERT_MSG((exp_a) < (exp_b), \
        (TO_S(exp_a) + " < " + TO_S(exp_b)).c_str())

#define ASSERT_LTE(exp_a, exp_b) \
    BOOST_ASSERT_MSG((exp_a) <= (exp_b), \
        (TO_S(exp_a) + " <= " + TO_S(exp_b)).c_str())

#define ASSERT_GT(exp_a, exp_b) \
    BOOST_ASSERT_MSG((exp_a) > (exp_b), \
        (TO_S(exp_a) + " > " + TO_S(exp_b)).c_str())

#define ASSERT_GTE(exp_a, exp_b) \
    BOOST_ASSERT_MSG((exp_a) >= (exp_b), \
        (TO_S(exp_a) + " >= " + TO_S(exp_b)).c_str())

