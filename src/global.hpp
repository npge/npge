/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_GLOBAL_HPP_
#define BR_GLOBAL_HPP_

#include <boost/shared_ptr.hpp>

/** Namespace for bloomrepeats */
namespace bloomrepeats {

class BloomFilter;
class Sequence;
class InMemorySequence;
class AnchorFinder;

typedef boost::shared_ptr<BloomFilter> BloomFilterPtr;
typedef boost::shared_ptr<Sequence> SequencePtr;
typedef boost::shared_ptr<AnchorFinder> AnchorFinderPtr;

}

#endif

