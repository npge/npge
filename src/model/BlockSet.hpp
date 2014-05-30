/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_BLOCK_SET_HPP_
#define NPGE_BLOCK_SET_HPP_

#include <iosfwd>
#include <set>
#include <vector>
#include <map>
#include "boost-xtime.hpp"
#include <boost/pool/pool_alloc.hpp>
#include <boost/utility.hpp>

#include "global.hpp"
#include "block_set_alignment.hpp"

namespace npge {

/** Container of blocks.
*/
class BlockSet : boost::noncopyable {
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

    /** Constructor */
    BlockSet();

    /** Destructor.
    Call clear().
    */
    ~BlockSet();

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
    std::vector<SequencePtr> seqs() const;

    /** Remove sequence */
    void remove_sequence(const SequencePtr& seq);

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

    /** Remove fragment.
    The block is not deleted.
    */
    void detach(Block* block);

    /** Return the number of blocks */
    size_t size() const;

    /** Return if there is no blocks */
    bool empty() const;

    /** Return if has the block */
    bool has(const Block* block) const;

    /** Remove all blocks and sequences.
    \see clear_blocks(), clear_seqs(), clear_bsas()
    */
    void clear();

    /** Remove all blocks.
    Removed blocks are deleted.
    */
    void clear_blocks();

    /** Remove all sequences.
    Removed sequences may remain, as they are handled by shared pointers.
    */
    void clear_seqs();

    /** Access block set alignment by its name.
    If name not present, new bsa is created.
    */
    BSA& bsa(const std::string& bsa_name);

    /** Access block set alignment by its name.
    If name not present, Exception is thrown.
    */
    const BSA& bsa(const std::string& bsa_name) const;

    /** List of names of block set alignments */
    Strings bsas() const;

    /** Return if block set has alignment with this name */
    bool has_bsa(const std::string& bsa_name) const;

    /** Remove block set alignment by its name.
    If name not present, do nothing.
    */
    void remove_bsa(const std::string& bsa_name);

    /** Remove all block set alignments */
    void clear_bsas();

    /** Exchange values of two objects */
    void swap(BlockSet& other);

    /** Copy the block set */
    BlockSetPtr clone() const;

    /** Copy blocks to other block set */
    void copy(BlockSet& target) const;

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

private:
    struct I;

    I* impl_;
};

/** Streaming operator.
\note Sequence list must be pre-added using BlockSet::add_sequence.

Alignment is not stored. COMPACT_ROW is used as row_type.
*/
std::istream& operator>>(std::istream& i, BlockSet& block_set);

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const BlockSet& block_set);

}

#endif

