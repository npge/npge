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

namespace boost {
namespace program_options {

// Missing forward declarations
class options_description;
class variables_map;

}
}

/** Namespace for bloomrepeats */
namespace bloomrepeats {

/** Namespace alias for program_options */
namespace po = boost::program_options;

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

class Exception;

}

#endif

