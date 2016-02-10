/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_DECIMAL_HPP_
#define NPGE_DECIMAL_HPP_

// TODO use Boost.Multiprecision

#include <stdint.h> // for int64_t
#include <iosfwd>
#include <cstdlib>
#include <cstdio>
#include <string>

#include "cast.hpp"

namespace npge {

template < typename Impl = int64_t,
         Impl SubPoint = 10000, int Digits = 4 >
class BasicDecimal {
public:
    /** Typede to self */
    typedef BasicDecimal<Impl, SubPoint> this_type;

    /** Base of number.
    For 4-digits number 1.2345, sub_point = 10000.
    */
    static const Impl sub_point = SubPoint;

    /** Number of subpoint digits */
    static const int digits = Digits;

    /** Constructor */
    BasicDecimal(int value = 0):
        impl_(Impl(value) * SubPoint) {
    }

    /** Constructor */
    BasicDecimal(const this_type& value):
        impl_(value.impl_) {
    }

    /** Constructor */
    BasicDecimal(const std::string& value) {
        int point = value.find('.');
        if (point != std::string::npos) {
            std::string c = value.substr(0, point);
            Impl i = c.empty() ? 0 : L_CAST<Impl>(c);
            std::string frac = value.substr(point + 1);
            frac.resize(Digits, '0');
            Impl j = L_CAST<Impl>(frac);
            if (i < 0) {
                j = -j;
            }
            impl_ = i * SubPoint + j;
        } else {
            impl_ = L_CAST<Impl>(value) * SubPoint;
        }
    }

    /** Operator */
    this_type& operator+=(const this_type& o) {
        impl_ += o.impl_;
        return *this;
    }

    /** Operator */
    this_type& operator-=(const this_type& o) {
        impl_ -= o.impl_;
        return *this;
    }

    /** Operator */
    this_type& operator*=(const this_type& o) {
        impl_ *= o.impl_;
        impl_ /= SubPoint;
        return *this;
    }

    /** Operator */
    this_type& operator/=(const this_type& o) {
        impl_ *= SubPoint;
        impl_ /= o.impl_;
        return *this;
    }

    /** Operator */
    this_type operator+(const this_type& o) const {
        this_type result;
        result.impl_ = impl_ + o.impl_;
        return result;
    }

    /** Operator */
    this_type operator-(const this_type& o) const {
        this_type result;
        result.impl_ = impl_ - o.impl_;
        return result;
    }

    /** Operator */
    this_type operator-() const {
        this_type result;
        result.impl_ = -impl_;
        return result;
    }

    /** Operator */
    this_type operator*(const this_type& o) const {
        this_type result;
        result.impl_ = impl_ * o.impl_ / SubPoint;
        return result;
    }

    /** Operator */
    this_type operator/(const this_type& o) const {
        this_type result;
        result.impl_ = impl_ * SubPoint / o.impl_;
        return result;
    }

    /** Operator */
    bool operator==(const this_type& o) const {
        return impl_ == o.impl_;
    }

    /** Operator */
    bool operator!=(const this_type& o) const {
        return impl_ != o.impl_;
    }

    /** Operator */
    bool operator<(const this_type& o) const {
        return impl_ < o.impl_;
    }

    /** Operator */
    bool operator<=(const this_type& o) const {
        return impl_ <= o.impl_;
    }

    /** Operator */
    bool operator>(const this_type& o) const {
        return impl_ > o.impl_;
    }

    /** Operator */
    bool operator>=(const this_type& o) const {
        return impl_ >= o.impl_;
    }

    /** Conver to double */
    double to_d() const {
        return double(impl_) / double(SubPoint);
    }

    /** Get integral part (floor).
    For negative x, return -(floor(-x)).
    */
    Impl to_i() const {
        if (impl_ >= 0) {
            return impl_ / SubPoint;
        } else {
            return -(-impl_ / SubPoint);
        }
    }

    /** Get fractional part.
    Example: D(-1.3).fraction() == 3 * Decimal::sub_point / 10
    */
    Impl fraction() const {
        return abs(impl_) % SubPoint;
    }

    /** Round number */
    Impl round() const {
        if (fraction() < sub_point / 2) {
            return to_i();
        } else {
            Impl r = to_i();
            r += (r < 0) ? (-1) : (1);
            return r;
        }
    }

    /** Conver to string */
    std::string to_s() const {
        std::string result = TO_S(to_i());
        Impl f = fraction();
        if (f) {
            std::string d = TO_S(fraction());
            int missing_zeros = Digits - d.length();
            result += '.';
            if (missing_zeros > 0) {
                result += std::string(missing_zeros, '0');
            }
            while (!d.empty() && d[d.size() - 1] == '0') {
                d.resize(d.size() - 1);
            }
            result += d;
        }
        return result;
    }

    Impl impl_;
};

typedef BasicDecimal<> Decimal;

#define D(...) ::npge::Decimal(#__VA_ARGS__)

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Decimal& d);

}

#endif

