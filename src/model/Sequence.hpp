/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SEQUENCE_HPP_
#define BR_SEQUENCE_HPP_

#include <iosfwd>
#include <string>
#include <vector>

#include "global.hpp"

namespace bloomrepeats {

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

    size_t size() const {
        return size_;
    }

    void print_contents(std::ostream& o) const;

    std::string contents() const;

    void make_first_fragment(Fragment& fragment, size_t fragment_size,
                             int only_ori = 1) const;

    bool next_fragment(Fragment& fragment) const;

    bool next_fragment_keeping_ori(Fragment& fragment) const;

    static void to_atgc(std::string& data);

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

    /** Return block, whose consensus is stored in this sequence */
    const Block* block() const {
        return block_;
    }

    /** Set block and store its consensus in this sequence.
    \warning Sequence must be empty (size() == 0).
    \note Block should be aligned to provide reliable consensus.
    \note Block's name is used as Sequence's name. See UniqueNames.
    */
    void set_block(const Block* block);

protected:
    virtual char char_at(size_t index) const = 0;

    void set_size(size_t size) {
        size_ = size;
    }

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
    virtual char char_at(size_t index) const;

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
    virtual char char_at(size_t index) const;

    void map_from_string_impl(const std::string& data, size_t min_pos);

private:
    std::string data_;

    void read_from_file(std::istream& input);

    void add_hunk(const std::string& hunk);

    void set_item(size_t index, char value);

    size_t byte_index(size_t index) const;

    size_t shift(size_t index) const;
};

}

#endif

