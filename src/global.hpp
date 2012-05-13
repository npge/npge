/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_GLOBAL_HPP_
#define BR_GLOBAL_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

/** Namespace for bloomrepeats */
namespace bloomrepeats {

class BloomFilter;
class Sequence;
class InMemorySequence;
class AnchorFinder;
class Fragment;
class Block;
class BlockSet;
class PairAligner;

typedef boost::shared_ptr<BloomFilter> BloomFilterPtr;
typedef boost::shared_ptr<Sequence> SequencePtr;
typedef boost::shared_ptr<AnchorFinder> AnchorFinderPtr;
typedef boost::shared_ptr<Fragment> FragmentPtr;
typedef boost::shared_ptr<Block> BlockPtr;
typedef boost::shared_ptr<BlockSet> BlockSetPtr;
typedef boost::shared_ptr<PairAligner> PairAlignerPtr;

}

#endif

