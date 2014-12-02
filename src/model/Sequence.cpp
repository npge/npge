/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <sstream>
#include <vector>
#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/utility/binary.hpp>

#include "Sequence.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "FastaReader.hpp"
#include "block_stat.hpp"
#include "complement.hpp"
#include "char_to_size.hpp"
#include "make_hash.hpp"
#include "read_block_set.hpp"
#include "name_to_stream.hpp"
#include "key_value.hpp"
#include "po.hpp"
#include "cast.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "global.hpp"

namespace npge {

Sequence::Sequence():
    size_(0), block_(0) {
}

SequencePtr Sequence::new_sequence(SequenceType seq_type) {
    if (seq_type == ASIS_SEQUENCE) {
        return boost::make_shared<InMemorySequence>();
    } else if (seq_type == COMPACT_LOW_N_SEQUENCE) {
        return boost::make_shared<CompactLowNSequence>();
    } else {
        return boost::make_shared<CompactSequence>();
    }
}

Sequence::~Sequence() {
}

void Sequence::map_from_string(const std::string& data,
                               pos_t min_pos) {
    if (!data.empty()) {
        pos_t new_size = min_pos + data.size();
        if (new_size > size()) {
            set_size(new_size);
        }
        map_from_string_impl(data, min_pos);
    }
}

void Sequence::push_back(const std::string& data) {
    read_from_string(data);
}

void Sequence::print_header(std::ostream& o) const {
    o << name();
    if (!description().empty()) {
        o << ' ' << description();
    }
}

void Sequence::print_contents(std::ostream& o, int line) const {
    for (int pos = 0; pos < size(); pos++) {
        if (line != 0 && pos % line == 0 && pos > 0) {
            o << '\n';
        }
        o << char_at(pos);
    }
}

std::string Sequence::contents() const {
    return substr(0, size(), 1);
}

void Sequence::make_first_fragment(Fragment& f,
                                   pos_t fragment_size,
                                   int only_ori) const {
    f.set_min_pos(0 - 1);
    f.set_max_pos(fragment_size - 1 - 1);
    f.set_ori(only_ori ? : 1);
}

bool Sequence::next_fragment(Fragment& f) const {
    f.set_ori(-f.ori());
    if (f.ori() == -1) {
        f.set_min_pos(f.min_pos() + 1);
        f.set_max_pos(f.max_pos() + 1);
    }
    return f.max_pos() < size();
}

bool Sequence::next_fragment_keeping_ori(Fragment& f) const {
    f.set_min_pos(f.min_pos() + 1);
    f.set_max_pos(f.max_pos() + 1);
    return f.max_pos() < size();
}

Fragment* Sequence::fragment_from_id(const std::string& id) {
    if (id.empty()) {
        return 0;
    }
    size_t u1 = id.find('_');
    if (u1 == std::string::npos) {
        return 0;
    }
    std::string seq_name = id.substr(0, u1);
    ASSERT_EQ(seq_name, name());
    size_t u2 = id.find('_', u1 + 1);
    if (u2 == std::string::npos) {
        return 0;
    }
    std::string begin_pos_str = id.substr(u1 + 1, u2 - u1 - 1);
    int begin_pos = L_CAST<int>(begin_pos_str);
    std::string last_pos_str = id.substr(u2 + 1);
    int last_pos = L_CAST<int>(last_pos_str);
    int ori = begin_pos <= last_pos ? 1 : -1;
    if (last_pos == -1) {
        last_pos = begin_pos;
        ori = -1;
    }
    Fragment* f = new Fragment(this);
    f->set_ori(ori);
    f->set_begin_pos(begin_pos);
    f->set_last_pos(last_pos);
    return f;
}

void Sequence::to_atgcn(std::string& data) {
    int s = data.size();
    char* d = &(data[0]);
    char* f = d;
    char* e = d + s;
    while (f != e) {
        char c = toupper(*f);
        if (c == 'A' || c == 'T' || c == 'G' || c == 'C' ||
                c == 'N') {
            *d = c;
            d++;
        }
        f++;
    }
    int removed_count = e - d;
    data.resize(s - removed_count);
}

void Sequence::set_name(const std::string& name) {
    if (name.find(' ') != std::string::npos) {
        throw Exception("Sequence name must not "
                        "include space: " + name);
    }
    if (is_fragment_name(name)) {
        throw Exception("Sequence name must not look like "
                        "fragment name: " + name);
    }
    name_ = name;
}

std::string Sequence::genome() const {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, name(), is_any_of("&"));
    if (parts.size() == 3 && (parts[2] == "c" || parts[2] == "l")) {
        return parts[0];
    } else {
        return "";
    }
}

std::string Sequence::chromosome() const {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, name(), is_any_of("&"));
    if (parts.size() == 3 && (parts[2] == "c" || parts[2] == "l")) {
        return parts[1];
    } else {
        return "";
    }
}

bool Sequence::circular() const {
    using namespace boost::algorithm;
    Strings parts;
    split(parts, name(), is_any_of("&"));
    if (parts.size() == 3 && (parts[2] == "c" || parts[2] == "l")) {
        return parts[2] == "c";
    } else {
        return false;
    }
}

std::string Sequence::ac() const {
    return extract_value(description(), "ac");
}

char Sequence::char_at(pos_t index) const {
    ASSERT_GTE(index, 0);
    ASSERT_MSG(index < size(), ("Index out of sequence:"
                                " name=" + name() +
                                " index=" + TO_S(index) +
                                ", size=" + TO_S(size())).c_str());
    return char_at_impl(index);
}

std::string Sequence::substr(pos_t index, pos_t length,
                             int ori) const {
    if (length == 0) {
        return "";
    }
    return substr_impl(index, length, ori);
}

std::string Sequence::substr_impl(pos_t index, pos_t length,
                                  int ori) const {
    ASSERT_LT(index, size());
    ASSERT_LT(index + (length - 1) * ori, size());
    std::string result;
    result.reserve(length);
    for (pos_t i = 0; i < length; i += 1) {
        char c = char_at_impl(index);
        if (ori == -1) {
            c = complement(c);
        }
        result += c;
        index += ori;
    }
    return result;
}

template<int ori>
class SChar {
public:
    pos_t index_;
    const Sequence* seq_;

    SChar(pos_t index, const Sequence* seq):
        index_(index), seq_(seq) {
    }

    char operator()(pos_t pos) {
        if (ori == 1) {
            return seq_->char_at_impl(index_ + pos);
        } else {
            return seq_->char_at_impl(index_ - pos);
        }
    }
};

hash_t Sequence::hash(pos_t index, pos_t length,
                      int ori) const {
    return hash_impl(index, length, ori);
}

hash_t Sequence::hash_impl(pos_t index, pos_t length,
                           int ori) const {
    ASSERT_LT(index, size());
    ASSERT_LT(index + (length - 1) * ori, size());
    if (ori == 1) {
        typedef SChar<1> F;
        return make_hash_base<F, 1>(F(index, this), length);
    } else {
        typedef SChar < -1 > F;
        return make_hash_base < F, -1 > (F(index, this), length);
    }
}

void Sequence::set_block(const Block* block,
                         bool set_consensus) {
    if (set_consensus) {
        ASSERT_EQ(size(), 0);
        ASSERT_EQ(block_, 0);
    }
    block_ = block;
    if (set_consensus) {
        std::string name_value = name();
        std::string description_value = description();
        std::stringstream cons;
        cons << ">dummy\n";
        block->consensus(cons);
        read_from_file(cons);
        set_name(name_value);
        set_description(description_value);
    }
    if (name().empty()) {
        set_name(block->name());
    }
    if (description().empty()) {
        AlignmentStat stat;
        make_stat(stat, block);
        std::string d;
        d += "fragments=";
        BOOST_FOREACH (Fragment* f, *block) {
            d += f->id();
            d += ',';
        }
        d.resize(d.size() - 1); // cut last comma
        d += " columns=" + TO_S(block->alignment_length());
        d += " identity=" + TO_S(block_identity(stat));
        set_description(d);
    }
}

InMemorySequence::InMemorySequence() {
}

InMemorySequence::InMemorySequence(const std::string& filename, int) {
    boost::shared_ptr<std::istream> file = name_to_istream(filename);
    read_from_file(*file);
}

InMemorySequence::InMemorySequence(std::istream& input) {
    read_from_file(input);
}

InMemorySequence::InMemorySequence(const std::string& data):
    data_(data) {
    to_atgcn(data_);
    set_size(data_.size());
}

char InMemorySequence::char_at_impl(pos_t index) const {
    return data_[index];
}

template <typename F>
class SequenceFastaReader : public FastaReader {
public:
    SequenceFastaReader(Sequence& seq, std::istream& input, const F& f):
        FastaReader(input), seq_(seq), f_(f) {
    }

    void new_sequence(const std::string& name, const std::string& description) {
        seq_.set_name(name);
        seq_.set_description(description);
    }

    void grow_sequence(const std::string& data) {
        std::string line(data);
        Sequence::to_atgcn(line);
        f_(line);
    }

private:
    Sequence& seq_;
    const F& f_;
};

template <typename F>
static void read_fasta(Sequence& seq, std::istream& input, const F& f) {
    SequenceFastaReader<F> reader(seq, input, f);
    reader.read_one_sequence();
}

void InMemorySequence::read_from_file(std::istream& input) {
    typedef std::string& (std::string::* StringMethod)(const std::string&);
    StringMethod append = &std::string::append;
    read_fasta(*this, input, boost::bind(append, &data_, _1));
    set_size(data_.size());
}

void InMemorySequence::read_from_string(const std::string& data) {
    data_ += data;
    to_atgcn(data_);
    set_size(data_.size());
}

void InMemorySequence::map_from_string_impl(const std::string& data,
        pos_t min_pos) {
    pos_t new_size = min_pos + data.size();
    if (new_size > data_.size()) {
        data_.resize(new_size);
    }
    for (pos_t i = 0; i < data.size(); i++) {
        data_[min_pos + i] = data[i];
    }
}

CompactSequence::CompactSequence() {
}

CompactSequence::CompactSequence(std::istream& input) {
    read_from_file(input);
}

CompactSequence::CompactSequence(const std::string& data) {
    std::string data_copy(data);
    to_atgcn(data_copy);
    add_hunk(data_copy);
}

const size_t LAST_2_BITS = BOOST_BINARY(11);
const size_t LAST_BIT = BOOST_BINARY(1);
const size_t SEQ_CHUNK_LETTERS = 8;
const size_t SEQ_CHUNK_BYTES = 3;
const size_t SEQ_BITS_PER_LETTER = 2;
const size_t SEQ_BITS_IN_BYTE = 8;

char CompactSequence::char_at_impl(pos_t index) const {
    pos_t n_i = n_index(index);
    if (n_i >= data_.size()) {
        return '\0';
    }
    size_t index_in_ch = index_in_chunk(index);
    if ((data_[n_i] >> index_in_ch) & LAST_BIT) {
        return 'N';
    }
    size_t contents_i = contents_index(index);
    size_t index_in_c = index_in_contents(index);
    size_t s = (data_[contents_i] >> index_in_c) & LAST_2_BITS;
    return size_to_char(s);
}

void CompactSequence::read_from_file(std::istream& input) {
    read_fasta(*this, input,
               boost::bind(&CompactSequence::add_hunk, this, _1));
}

void CompactSequence::read_from_string(const std::string& data) {
    std::string data_copy(data);
    to_atgcn(data_copy);
    add_hunk(data_copy);
}

void CompactSequence::map_from_string_impl(const std::string& data,
        pos_t min_pos) {
    if (data.empty()) {
        return;
    }
    size_t new_size = min_pos + data.size();
    size_t chunk_needed = chunk_index(new_size - 1);
    size_t chunks_exist = data_.size() / SEQ_CHUNK_BYTES;
    if (chunk_needed >= chunks_exist) {
        data_.resize(SEQ_CHUNK_BYTES * (chunk_needed + 1));
    }
    for (size_t i = 0; i < data.size(); i++) {
        set_item(min_pos + i, data[i]);
    }
}

void CompactSequence::add_hunk(const std::string& hunk) {
    pos_t new_size = size() + hunk.size();
    map_from_string_impl(hunk, size());
    set_size(new_size);
}

void CompactSequence::set_item(size_t index, char value) {
    if (value == 'N') {
        size_t n_i = n_index(index);
        if (n_i >= data_.size()) {
            return;
        }
        size_t index_in_ch = index_in_chunk(index);
        data_[n_i] |= 1 << index_in_ch;
    } else {
        size_t contents_i = contents_index(index);
        if (contents_i >= data_.size()) {
            return;
        }
        size_t index_in_c = index_in_contents(index);
        data_[contents_i] |= char_to_size(value) << index_in_c;
    }
}

size_t CompactSequence::chunk_index(size_t index) const {
    return index / SEQ_CHUNK_LETTERS;
}

size_t CompactSequence::n_index(size_t index) const {
    return chunk_index(index) * SEQ_CHUNK_BYTES;
}

size_t CompactSequence::contents_index(size_t index) const {
    size_t index_in_ch = index_in_chunk(index);
    int byte_index = index_in_ch  / (SEQ_CHUNK_LETTERS / 2); // 0 or 1
    return chunk_index(index) * SEQ_CHUNK_BYTES
           + 1 // first byte - with Ns
           + byte_index;
}

size_t CompactSequence::index_in_chunk(size_t index) const {
    return index % SEQ_CHUNK_LETTERS;
}

size_t CompactSequence::index_in_contents(size_t index) const {
    return (SEQ_BITS_PER_LETTER * index_in_chunk(index)) %
           SEQ_BITS_IN_BYTE;
}

CompactLowNSequence::CompactLowNSequence()
{ }

CompactLowNSequence::CompactLowNSequence(
    const std::string& data) {
    read_from_string(data);
}

char CompactLowNSequence::char_at_impl(pos_t index) const {
    if (ns_.has_elem(index)) {
        return 'N';
    }
    size_t s = (data_[byte_index(index)] >> shift(index)) &
               LAST_2_BITS;
    return size_to_char(s);
}

void CompactLowNSequence::read_from_file(std::istream& input) {
    read_fasta(*this, input,
               boost::bind(&CompactLowNSequence::add_hunk,
                           this, _1));
}

void CompactLowNSequence::read_from_string(
    const std::string& data) {
    std::string data_copy(data);
    to_atgcn(data_copy);
    add_hunk(data_copy);
}

void CompactLowNSequence::map_from_string_impl(
    const std::string&, pos_t) {
    throw Exception("CompactLowNSequence::map_from_string "
                    "not implemented");
}

void CompactLowNSequence::add_hunk(const std::string& hunk) {
    if (hunk.empty()) {
        return;
    }
    size_t n_pos = -1;
    while (true) {
        n_pos = hunk.find('N', n_pos + 1);
        if (n_pos == std::string::npos) {
            break;
        }
        ns_.push_back(size() + n_pos);
    }
    size_t new_size = size() + hunk.size();
    if (byte_index(new_size - 1) >= data_.size()) {
        data_.resize(byte_index(new_size - 1) + 1);
    }
    for (size_t i = 0; i < hunk.size(); i++) {
        set_item(size() + i, hunk[i]);
    }
    set_size(new_size);
}

void CompactLowNSequence::set_item(size_t index, char value) {
    if (value != 'N') {
        data_[byte_index(index)] |=
            char_to_size(value) << shift(index);
    }
}

size_t CompactLowNSequence::byte_index(size_t index) const {
    return index / 4;
}

size_t CompactLowNSequence::shift(size_t index) const {
    return 2 * (index % 4);
}

DummySequence::DummySequence(char letter, int size) {
    set_letter(letter);
    set_size(size);
}

char DummySequence::letter() const {
    return letter_;
}

void DummySequence::set_letter(char letter) {
    letter_ = letter;
}

void DummySequence::read_from_string(const std::string& data) {
    set_size(data.size());
}

char DummySequence::char_at_impl(pos_t) const {
    return letter_;
}

void DummySequence::map_from_string_impl(const std::string&,
        pos_t) {
}

void DummySequence::read_from_file(std::istream& input) {
    read_fasta(*this, input,
               boost::bind(&DummySequence::add_hunk,
                           this, _1));
}

void DummySequence::add_hunk(const std::string& hunk) {
    size_t new_size = size() + hunk.size();
    set_size(new_size);
}

FragmentSequence::FragmentSequence(Fragment* fragment) {
    set_fragment(fragment);
}

Fragment* FragmentSequence::fragment() const {
    return fragment_;
}

void FragmentSequence::set_fragment(Fragment* fragment) {
    fragment_ = fragment;
    if (fragment) {
        set_size(fragment->length());
    }
}

void FragmentSequence::read_from_string(const std::string&) {
    throw Exception("Trying to modify const FragmentSequence");
}

char FragmentSequence::char_at_impl(pos_t index) const {
    ASSERT_TRUE(fragment_);
    ASSERT_LT(index, fragment_->length());
    return fragment_->raw_at(index);
}

void FragmentSequence::map_from_string_impl(const std::string&,
        pos_t min_pos) {
    throw Exception("Trying to modify const FragmentSequence");
}

void FragmentSequence::read_from_file(std::istream& input) {
    throw Exception("Trying to modify const FragmentSequence");
}

std::ostream& operator<<(std::ostream& o, const Sequence& s) {
    o << '>';
    s.print_header(o);
    o << '\n';
    s.print_contents(o);
    o << '\n';
    return o;
}

}

