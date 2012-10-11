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

// util
class FastaReader;

// model
class Sequence;
class InMemorySequence;
class CompactSequence;
class Fragment;
struct FragmentDiff;
class Block;
class BlockSet;
class Alignment;
class AlignmentRow;

// algo
class BloomFilter;
class PairAligner;
class ExpanderBase;
// processors
class Processor;
class AnchorFinder;
class Filter;
class Connector;
class OverlapsResolver;
class Joiner;
class Union;
class BlocksExpander;
class FragmentsExpander;
class Pipe;
class CleanUp;
class Output;
class UniqueNames;
class Rest;

/** Shared pointer to BloomFilter */
typedef boost::shared_ptr<BloomFilter> BloomFilterPtr;

/** Shared pointer to Sequence */
typedef boost::shared_ptr<Sequence> SequencePtr;

/** Shared pointer to BlockSet */
typedef boost::shared_ptr<BlockSet> BlockSetPtr;

/** Shared pointer to PairAligner */
typedef boost::shared_ptr<PairAligner> PairAlignerPtr;

/** Shared pointer to Processor */
typedef boost::shared_ptr<Processor> ProcessorPtr;

class Exception;

}

#endif

