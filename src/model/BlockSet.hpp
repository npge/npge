/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_SET_HPP_
#define BR_BLOCK_SET_HPP_

#include <iosfwd>
#include <set>
#include <vector>
#include <map>
#include <boost/pool/pool_alloc.hpp>
#include <boost/enable_shared_from_this.hpp> // FIXME

#include "global.hpp"

namespace bloomrepeats {

/** Container of blocks.
*/
class BlockSet :
    public boost::enable_shared_from_this<BlockSet> { // FIXME
public:
    /** Type of implementation container.
    Do not rely on ths type!
    */
    typedef std::set < Block*, std::less<Block*>,
            boost::fast_pool_allocator<void> > Impl;

    /** Iterator */
    typedef Impl::iterator iterator;

    /** Constant iterator */
    typedef Impl::const_iterator const_iterator;

    /** Destructor.
    Call clear().
    */
    ~BlockSet();

    /** Return a copy of this block set.
    Fragments and blocks are copied, sequences are not copied,
    sequence list is copied.
    Connections between the fragments
    (\ref Fragment::prev() "prev", \ref Fragment::next() "next")
    are rebuild with connect_fragments().
    \todo Preserve fragment connections from source block set.
    \see Block::clone()
    */
    BlockSetPtr clone() const;

    /** Add sequence.
    All the sequences, used by blocks, must be added.
    The sequence shared pointer is guaranteed to be kept until the set exists.
    */
    void add_sequence(SequencePtr seq);

    /** Add all sequence from the vector.
    \see add_sequence
    */
    void add_sequences(const std::vector<SequencePtr>& sequences);

    /** Get sequence list */
    const std::vector<SequencePtr>& seqs() const {
        return seqs_;
    }

    /** Return sequence with given name or null */
    SequencePtr seq_from_name(const std::string& name) const;

    /** Create a fragment based on its ID.
    Returns valid fragment such that fragment->id() == id.
    Sequence must be pre-added.
    On error return 0.
    \warning This function will create new instance even if
        fragment is already in block set.
    */
    Fragment* fragment_from_id(const std::string& id) const;

    /** Extract block name from description of sequnce in FASTA file.
    This method can be used to get block name from description of FASTA
    files, created by streaming BlockSet to file.

    On error return "".
    */
    static std::string block_from_description(const std::string& description);

    /** Add block.
    The same block can't be added twice.
    */
    void insert(Block* block);

    /** Remove fragment.
    The block is deleted.
    */
    void erase(Block* block);

    /** Return the number of blocks */
    size_t size() const;

    /** Return if there is no blocks */
    bool empty() const;

    /** Return if has the block */
    bool has(const Block* block) const;

    /** Remove all blocks.
    Removed blocks are deleted.
    */
    void clear();

    /** Get some block if any or an empty pointer */
    Block* front() const;

    /** Return iterator to beginning */
    iterator begin();

    /** Return constant iterator to beginning */
    const_iterator begin() const;

    /** Return iterator to end */
    iterator end();

    /** Return constant iterator to end */
    const_iterator end() const;

    /** Connect all the fragments (prev-next) */
    void connect_fragments();

    /** Expand all blocks (starting from blocks of large number of fragments).
    \see Block::expand()
    */
    void expand_blocks(PairAligner* aligner = 0, int batch = 100, int ori = 0,
                       int max_overlap = 0);

    /** Return if there are blocks which have overlaping fragments.
    \warning Fragments must be \ref BlockSet::connect_fragments "connected"
    */
    bool overlaps() const;

    /** Resolve overlaping fragments.
    If some blocks from the block set have overlaping fragments
    these two blocks are replaced with one higher (and narrower)
    block and several remainder blocks.

    Anyway, applying this method guarantees that no blocks of the block set
    have overlaping fragments.

    Since resolve_overlaps() can split blocks,
    applying \ref join "join(0)" is recommended.

    \verbatim
    Input:
        Block 1:
            seq1: ---xxxx----
            seq2: ---xxxx----
        Block 2
            seq2: -----xxxx--
            seq3: -----xxxx--
            seq4: -----xxxx--

    Output of resolve_overlaps:
        Block 1:
            seq1: ---xx------
            seq2: ---xx------
        Block 2
            seq2: -------xx--
            seq3: -------xx--
            seq4: -------xx--
        Block 3:
            seq1: -----xx----
            seq2: -----xx----
            seq3: -----xx----
            seq4: -----xx----
    \endverbatim
    \warning Fragments must be \ref BlockSet::connect_fragments "connected"
       for this to work correctly.
    */
    void resolve_overlaps();

    /** Expand all blocks by fragments.
    Return \p true is something was added.
    \see Block::expand_by_fragments().
    */
    bool expand_blocks_by_fragments(PairAligner* aligner = 0, int batch = 100);

    /** Return new block set of blocks of nucleotides, not included in this set.
    From sequences involved in this block set, nucleotides are selected,
    not included in this block set. They are grouped into fragments.
    Each fragment is inserted into one block.
    These blocks are inserted into resulting block set.
    \warning Fragments must be \ref BlockSet::connect_fragments "connected"
       for this to work correctly.
    */
    BlockSetPtr rest() const;

    /** Add options to options description */
    static void add_pangenome_options(po::options_description& desc);

    /** Apply several methods to convert the block set to pangenome.
    This method could be applied to block set produced by AnchorFinder.
    */
    void make_pangenome(const po::variables_map& vm);

    /** Add options to options description */
    static void add_output_options(po::options_description& desc);

    /** Write all blocks to file(s) or std::cout (depends on vm) */
    void make_output(const po::variables_map& vm);

    /** Set unique names to all blocks of this block set.
    Firstly, Block::set_name_from_fragments() is used, if name is null.
    Then Block::set_random_name() is called untill the name is unique.
    */
    void set_unique_block_names();

private:
    typedef std::map<std::string, SequencePtr> Name2Seq;

    Impl blocks_;
    std::vector<SequencePtr> seqs_;
    mutable Name2Seq name2seq_;

    friend class BlockSetFastaReader;
};

/** Streaming operator.
\note Sequence list must be pre-added using BlockSet::add_sequence.
*/
std::istream& operator>>(std::istream& i, BlockSet& block_set);

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const BlockSet& block_set);

}

#endif

