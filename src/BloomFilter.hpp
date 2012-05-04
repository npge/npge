/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOOM_FILTER_HPP_
#define BR_BLOOM_FILTER_HPP_

#include <cstddef>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

namespace bloomrepeats {

/** Bloom filter.

See http://en.wikipedia.org/wiki/Bloom_filter
*/
class BloomFilter {
public:
    /** Shared pointer */
    typedef boost::shared_ptr<BloomFilter> Ptr;

    /** Default constructor.
    Postconditions: bits() = 0, hashes() = 0.
    */
    BloomFilter();

    /** Constructor.
    \see set_members, set_optimal_hashes
    */
    BloomFilter(size_t members, float error_prob);

    /** Set optimal bits number.
    \see optimal_bits(), set_bits()
    */
    void set_members(size_t members, float error_prob);

    /** Get bits number */
    size_t bits() const;

    /** Set bits number.
    \warning This method clears all added members.
    */
    void set_bits(size_t bits);

    /** Set optimal hash functions number.
    \see optimal_hashes(), set_hashes()
    */
    void set_optimal_hashes(size_t members);

    /** Get hash functions number */
    size_t hashes() const;

    /** Set hash functions number.
    \warning This method removes all existing hash functions
        and invalidates all added members. If you have added members,
        you should not call this method, unless you clear bits too.
    */
    void set_hashes(size_t hashes);

    /** Return if the member is likely to be added and add it.
    It is an equivalent to:
    \code
    bool was_added = test(start, length);
    add(start, length);
    return was_added;
    \endcode
    */
    bool test_and_add(const char* start, size_t length);

    /** Return if the member is likely to be added and add it.
    Overloaded method.
    */
    bool test_and_add(const std::string& member);

    /** Add member.
    \note Bytes are added as is, i.e, case sensitive.
    */
    void add(const char* start, size_t length);

    /** Add member.
    Overloaded method.
    */
    void add(const std::string& member);

    /** Return if the member is likely to be added.
    If returns false, then the member was added;
    if returns true, then the member is likely to be added.

    \see set_members
    */
    bool test(const char* start, size_t length) const;

    /** Return if the member is likely to be added.
    Overloaded method.
    */
    bool test(const std::string& member) const;

    /** Return optimal bits number.
    Optimal bits number is based on
    expected members number and false positive probability:
    \f$ m=-n\frac{\ln p}{(\ln 2)^2} \f$,
    where m is bits number,
    n is expected members number and
    p is false positive probability.
    */
    static size_t optimal_bits(size_t members, float error_prob);

    /** Return optimal hash functions number.
    Optimal hash functions number is based on
    expected members number and bits number.
    \f$ k=\frac{m}{n} \ln 2 \f$,
    where k is optimal hashes number,
    m is bits number,
    n is expected members number.
    */
    static size_t optimal_hashes(size_t members, size_t bits);

private:
    std::vector<bool> bits_;
    std::vector<size_t> hash_mul_;

    size_t make_index(size_t hash, const char* start, size_t length) const;
};

}

#endif

