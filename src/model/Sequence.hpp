/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SEQUENCE_HPP_
#define NPGE_SEQUENCE_HPP_

#include <iosfwd>
#include <string>
#include <vector>

#include "global.hpp"
#include "boundaries.hpp"

namespace npge {

/** Config = Sequence */
class Sequence {
public:
    static const int FIRST_ORI = -1;

    Sequence();

    static SequencePtr new_sequence(SequenceType seq_type);

    virtual ~Sequence();

    virtual void read_from_file(std::istream& input) = 0;

    virtual void read_from_string(const std::string& data) = 0;

    void map_from_string(const std::string& data, size_t min_pos);

    void push_back(const std::string& data);

    size_t size() const {
        return size_;
    }

    void set_size(size_t size) {
        size_ = size;
    }

    void print_header(std::ostream& o) const;

    void print_contents(std::ostream& o, int line = 60) const;

    std::string contents() const;

    void make_first_fragment(Fragment& fragment, size_t fragment_size,
                             int only_ori = 1) const;

    bool next_fragment(Fragment& fragment) const;

    bool next_fragment_keeping_ori(Fragment& fragment) const;

    /** Create a fragment based on its ID.
    Returns valid fragment such that fragment->id() == id.
    On error return 0.
    */
    Fragment* fragment_from_id(const std::string& id);

    /** Upper-case, s/[^atgcn]// */
    static void to_atgcn(std::string& data);

    const std::string& name() const {
        return name_;
    }

    void set_name(const std::string& name) {
        name_ = name;
    }

    const std::string& description() const {
        return description_;
    }

    void set_description(const std::string& description) {
        description_ = description;
    }

    /** Return name of genome, if can be deduced from name().
    Format: genome&chromosome&circular.
    Empty string is returned if name format is wrong.
    */
    std::string genome() const;

    /** Return name of chromosome, if can be deduced from name().
    Format: genome&chromosome&circular.
    Empty string is returned if name format is not accepted.
    */
    std::string chromosome() const;

    /** Return if the contig is circular, if can be deduced from name().
    Format: genome&chromosome&circular.
    Circular = c (for circular) or l (for linear).
    Exception std::logic_error is thrown if format is wrong of circular
    is not 'c' or 'l'.
    */
    bool circular() const;

    /** Return accession number of the sequence.
    Accession number is set if description contains "ac=XXXX".
    If no accession number was set, empty string is returned.
    */
    std::string ac() const;

    /** Return block, whose consensus is stored in this sequence */
    const Block* block() const {
        return block_;
    }

    /** Set block and store its consensus in this sequence.
    \param block Target block
    \param set_consensus Whether to set block's consensus
        as Sequence's contents.
    \warning If set_consensus = true,
        sequence must be empty (size() == 0).
    \note Block should be aligned to provide
        reliable consensus.
    \note Block's name is used as Sequence's name.
        See UniqueNames.
    */
    void set_block(const Block* block,
                   bool set_consensus = true);

    char char_at(size_t index) const;

    // TODO std::string slice()

protected:
    virtual char char_at_impl(size_t index) const = 0;

    virtual void map_from_string_impl(const std::string& data,
                                      size_t min_pos) = 0;

private:
    size_t size_;
    std::string name_;
    std::string description_;
    const Block* block_;

    friend class Fragment;
};

class InMemorySequence : public Sequence {
public:
    InMemorySequence();

    // reads first sequence
    InMemorySequence(const std::string& filename, int);

    InMemorySequence(std::istream& input);

    InMemorySequence(const std::string& data);

    void read_from_string(const std::string& data);

protected:
    char char_at_impl(size_t index) const;

    void map_from_string_impl(const std::string& data, size_t min_pos);

private:
    std::string data_;

    void read_from_file(std::istream& input);
};

class CompactSequence : public Sequence {
public:
    CompactSequence();

    CompactSequence(std::istream& input);

    CompactSequence(const std::string& data);

    void read_from_string(const std::string& data);

protected:
    char char_at_impl(size_t index) const;

    void map_from_string_impl(const std::string& data, size_t min_pos);

private:
    std::string data_;

    void read_from_file(std::istream& input);

    void add_hunk(const std::string& hunk);

    void set_item(size_t index, char value);

    size_t chunk_index(size_t index) const;

    size_t n_index(size_t index) const;

    size_t contents_index(size_t index) const;

    size_t index_in_chunk(size_t index) const;

    size_t index_in_contents(size_t index) const;
};

class CompactLowNSequence : public Sequence {
public:
    CompactLowNSequence();

    CompactLowNSequence(const std::string& data);

    void read_from_string(const std::string& data);

protected:
    char char_at_impl(size_t index) const;

    void map_from_string_impl(const std::string& data, size_t min_pos);

private:
    std::string data_;
    Boundaries ns_;

    void read_from_file(std::istream& input);

    void add_hunk(const std::string& hunk);

    void set_item(size_t index, char value);

    size_t byte_index(size_t index) const;

    size_t shift(size_t index) const;
};

/** Sequence returning the one letter for each position.
This utility sequence can be used to use in place of long
sequences without large memory allocations.
*/
class DummySequence : public Sequence {
public:
    DummySequence(char letter = 'N', int size = 0);

    char letter() const;

    void set_letter(char letter);

    void read_from_string(const std::string& data);

protected:
    char char_at_impl(size_t index) const;

    void map_from_string_impl(const std::string& data,
                              size_t min_pos);

private:
    char letter_;

    void read_from_file(std::istream& input);

    void add_hunk(const std::string& hunk);
};

class FragmentSequence : public Sequence {
public:
    FragmentSequence(Fragment* fragment = 0);

    Fragment* fragment() const;

    void set_fragment(Fragment* fragment);

    void read_from_string(const std::string& data);

protected:
    char char_at_impl(size_t index) const;

    void map_from_string_impl(const std::string& data,
                              size_t min_pos);

private:
    Fragment* fragment_;

    void read_from_file(std::istream& input);
};

/** Streaming operator.
\see Sequence::print_header
\see Sequence::print_contents
*/
std::ostream& operator<<(std::ostream& o, const Sequence& sequence);

}

#endif

