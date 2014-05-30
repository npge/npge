/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GLOBAL_HPP_
#define NPGE_GLOBAL_HPP_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace boost {
namespace program_options {

// Missing forward declarations
class options_description;
class variables_map;

}
}

/** Namespace for NPG-explorer */
namespace npge {

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

class Processor;
class Pipe;
class BlocksJobs;

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
typedef std::vector<std::string> Strings;

/** Pair from alignment */
typedef std::pair<int, int> AlignmentPair;

/** Pair alignment, array of pairs of positions, -1 as gaps */
typedef std::vector<AlignmentPair> PairAlignment;

}

#endif

