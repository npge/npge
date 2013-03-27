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
class StringToArgv;
class FastaReader;
template<typename V> class SortedVector;
template<typename V> class Graph;

// model
class Sequence;
class InMemorySequence;
class CompactSequence;
class Fragment;
struct FragmentDiff;
struct AlignmentStat;
class Block;
class BlockSet;
class AlignmentRow;
class MapAlignmentRow;
class CompactAlignmentRow;
class BlockSetFastaReader;

// algo
class BloomFilter;
class PairAligner;
class ExpanderBase;
class SeqStorage;
class RowStorage;
class FileReader;
class FileWriter;
class OptionsPrefix;
class SizeLimits;
class Meta;
// processors
class Processor;
class AddSequences;
class SequencesFromOther;
class AddBlocks;
class BlastRunner;
class ImportBlastHits;
class AddBlastBlocks;
class ResolveBlast;
class ResolveAnchors;
class AnchorFinder;
class Filter;
class Connector;
class OriByMajority;
class StickBoundaries;
class OverlapsResolver;
class OverlapsResolver2;
class CheckNoOverlaps;
class Joiner;
class Union;
class BlocksExpander;
class FragmentsExpander;
class BlocksJobs;
class Pipe;
class CleanUp;
class AbstractOutput;
class Output;
class OutputPipe;
class PrintOverlaps;
class BlockInfo;
class Stats;
class Consensus;
class UniqueNames;
class Rest;
class ExternalAligner;
class ConSeq;
class DeConSeq;
class MoveGaps;
class CutGaps;

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

/** Type of Sequence */
enum SequenceType {
    NO_SEQUENCE, /**< No sequence */
    ASIS_SEQUENCE, /**< InMemorySequence */
    COMPACT_SEQUENCE /**< CompactSequence */
};

/** Type of AlignmentRow */
enum RowType {
    MAP_ROW, /**< MapAlignmentRow */
    COMPACT_ROW /**< CompactAlignmentRow */
};

/** Creat new BlockSet and return shared pointer to it */
BlockSetPtr new_bs();

}

#endif

