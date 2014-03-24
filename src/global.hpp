/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_GLOBAL_HPP_
#define BR_GLOBAL_HPP_

#include <vector>
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
class AnyAs;
template<typename V> class SortedVector;
template<typename V> class Graph;
class ThreadTask;
class ThreadWorker;
class ThreadGroup;
class TreeNode;
class LeafNode;
template<typename Contents> class GeneralAligner;

// model
class Sequence;
class InMemorySequence;
class CompactSequence;
class Fragment;
struct FragmentDiff;
class AlignmentStat;
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
class FileReader;
class FileWriter;
template<typename F, typename C> class FragmentCollection;
class Meta;
// processors
class Processor;
class SequencesFromOther;
class AddBlocks;
class BlastRunner;
class ImportBlastHits;
class AddBlastBlocks;
class BlastFinder;
class ResolveBlast;
class ResolveAnchors;
class AnchorFinder;
class Filter;
class Stem;
class SameChr;
class Connector;
class OriByMajority;
class StickBoundaries;
class OverlapsResolver;
class OverlapsResolver2;
class CheckNoOverlaps;
class SelfOverlapsResolver;
class Joiner;
class Union;
class Move;
class Clear;
class OverlaplessUnion;
class OneByOne;
class Partition;
class FindGeneGroups;
class FindGeneConversion;
class Subtract;
class BlocksExpander;
class FragmentsExpander;
class BlocksJobs;
class Hash;
class FileRemover;
class FileCopy;
class Pipe;
class CleanUp;
class AbstractOutput;
class Output;
class OutputPipe;
class FragmentDistance;
class PrintTree;
class ConsensusTree;
class BlockSetAlignment;
class ChrBlockSetAlignment;
class PrintBlockSetAlignment;
class InputBlockSetAlignment;
class FastaBlockSetAlignment;
class PrintOverlaps;
class PrintPartition;
class PrintGeneGroups;
class PrintMutations;
class MutationsSequences;
class BlockInfo;
class Stats;
class Info;
class IsPangenome;
class MakePrePangenome;
class MakePangenome;
class Consensus;
class UniqueNames;
class Rest;
class ExternalAligner;
class RemoveAlignment;
class RemoveNames;
class MarkNonWeak;
class LinkEqualFragments;
class ConSeq;
class DeConSeq;
class MoveGaps;
class CutGaps;
class LiteAlign;
class Align;
class MetaProcessor;
class TrySmth;
class AddGenes;

/** Shared pointer to BloomFilter */
typedef boost::shared_ptr<BloomFilter> BloomFilterPtr;

/** Shared pointer to Sequence */
typedef boost::shared_ptr<Sequence> SequencePtr;

/** Shared pointer to BlockSet */
typedef boost::shared_ptr<BlockSet> BlockSetPtr;

/** Shared pointer to PairAligner */
typedef boost::shared_ptr<PairAligner> PairAlignerPtr;

/** Shared pointer to Processor */
typedef boost::shared_ptr<Processor> SharedProcessor;

/** Shared pointer to Pipe */
typedef boost::shared_ptr<Pipe> SharedPipe;

class Exception;

/** Type of Sequence */
enum SequenceType {
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

typedef std::vector<Fragment*> Fragments;
typedef std::vector<Block*> Blocks;

}

#endif

