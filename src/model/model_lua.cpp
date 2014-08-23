/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <luabind/luabind.hpp>
#include <luabind/operator.hpp>
#include <luabind/iterator_policy.hpp>
#include <luabind/out_value_policy.hpp>
#include <luabind/adopt_policy.hpp>

#include "model_lua.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "block_set_alignment.hpp"
#include "block_stat.hpp"
#include "block_hash.hpp"
#include "convert_position.hpp"

namespace npge {

static SequencePtr new_sequence0(SequenceType type) {
    return Sequence::new_sequence(type);
}

static SequencePtr new_sequence1() {
    return new_sequence0(COMPACT_SEQUENCE);
}

static std::string sequence_to_atgcn(const std::string& text) {
    std::string copy = text;
    Sequence::to_atgcn(copy);
    return copy;
}

static Fragment* new_fragment0() {
    return new Fragment;
}

static Fragment* new_fragment1(Sequence* seq) {
    return new Fragment(seq);
}

static Fragment* new_fragment2(Sequence* seq,
                               size_t min_pos, size_t max_pos) {
    return new Fragment(seq, min_pos, max_pos);
}

static Fragment* new_fragment3(Sequence* seq,
                               size_t min_pos, size_t max_pos, int ori) {
    return new Fragment(seq, min_pos, max_pos, ori);
}

static void delete_fragment(Fragment* f) {
    delete f;
}

static void fragment_inverse(Fragment* f) {
    f->inverse();
}

static void fragment_set_ori(Fragment* f, int ori) {
    f->set_ori(ori);
}

static std::string fragment_str(Fragment* f) {
    return f->str();
}

static std::string fragment_header(Fragment* f) {
    std::stringstream ss;
    f->print_header(ss);
    return ss.str();
}

static std::string fragment_header2(Fragment* f,
                                    const Block* block) {
    std::stringstream ss;
    f->print_header(ss, block);
    return ss.str();
}

static std::string fragment_contents(Fragment* f) {
    std::stringstream ss;
    f->print_contents(ss);
    return ss.str();
}

static std::string fragment_contents2(Fragment* f, int line) {
    std::stringstream ss;
    f->print_contents(ss, '-', line);
    return ss.str();
}

static AlignmentRow* alignmentrow_new_row() {
    return new CompactAlignmentRow;
}

static void alignmentrow_delete_row(AlignmentRow* row) {
    delete row;
}

static Block* new_block0() {
    return new Block;
}

static Block* new_block1(const std::string& name) {
    return new Block(name);
}

static void delete_block(Block* block) {
    delete block;
}

static const Block& block_fragments(const Block* block) {
    return *block;
}

static char block_consensus_char(const Block* block, int pos) {
    return block->consensus_char(pos);
}

static std::string block_consensus_string(const Block* block) {
    return block->consensus_string();
}

static void block_inverse(Block* block) {
    block->inverse();
}

static int alignmentstat_letter_count(const AlignmentStat& stat,
                                      const std::string& letter) {
    ASSERT_EQ(letter.size(), 1);
    return stat.letter_count(letter[0]);
}

static void make_stat0(AlignmentStat& stat,
                       const Block* block) {
    make_stat(stat, block, 0, -1);
}

static Decimal alignmentstat_identity(
    const AlignmentStat& stat) {
    return block_identity(stat);
}

static Decimal block_identity0(int ident_nogap, int ident_gap,
                               int noident_nogap, int noident_gap) {
    return block_identity(ident_nogap, ident_gap,
                          noident_nogap, noident_gap);
}

static void lite_test_column(const Block* block, int column,
                             bool& ident, bool& gap) {
    test_column(block, column, ident, gap);
}

static std::string block_hash_str(const Block* block) {
    return TO_S(block_hash(block));
}

static const BlockSet& blockset_blocks(const BlockSet& bs) {
    return bs;
}

static void blockset_read(BlockSet& bs,
                          const std::string& text) {
    std::stringstream ss(text);
    ss >> bs;
}

static std::string blockset_hash1(BlockSet& bs) {
    return TO_S(blockset_hash(bs));
}

static BSA& blockset_bsa(BlockSet& bs, const std::string& k) {
    return bs.bsa(k);
}

typedef std::vector<SequencePtr> Sequences;

static const Sequences& sequences_iter(const Sequences& seqs) {
    return seqs;
}

static const Fragments& fragments_iter(const Fragments& ff) {
    return ff;
}

static Fragment* fragments_at(const Fragments& v, int index) {
    return v.at(index);
}

static const Blocks& blocks_iter(const Blocks& bb) {
    return bb;
}

static Block* blocks_at(const Blocks& v, int index) {
    return v.at(index);
}

static bool bsa_has(const BSA& bsa, Sequence* seq) {
    return bsa.find(seq) != bsa.end();
}

static BSRow& bsa_get(BSA& bsa, Sequence* seq) {
    return bsa[seq];
}

static void bsa_erase(BSA& bsa, Sequence* seq) {
    bsa.erase(seq);
}

}

extern "C" int init_model_lua(lua_State* L) {
    using namespace luabind;
    using namespace npge;
    open(L);
    module(L) [
        class_<Sequence, SequencePtr>("Sequence")
        .enum_("Type") [
            value("ASIS_SEQUENCE", ASIS_SEQUENCE),
            value("COMPACT_SEQUENCE", COMPACT_SEQUENCE),
            value("COMPACT_LOW_N_SEQUENCE",
                  COMPACT_LOW_N_SEQUENCE)
        ]
        .scope [
            def("new", &new_sequence0),
            def("new", &new_sequence1),
            def("to_atgcn", &sequence_to_atgcn)
        ]
        .def("push_back", &Sequence::push_back)
        .def("size", &Sequence::size)
        .def("set_size", &Sequence::set_size)
        .def("contents", &Sequence::contents)
        .def("name", &Sequence::name)
        .def("set_name", &Sequence::set_name)
        .def("description", &Sequence::description)
        .def("set_description", &Sequence::set_description)
        .def("genome", &Sequence::genome)
        .def("chromosome", &Sequence::chromosome)
        .def("circular", &Sequence::chromosome)
        .def("ac", &Sequence::ac)
        .def("char_at", &Sequence::char_at)
        .def("substr", &Sequence::substr)
        .def("hash", &Sequence::hash)
        .def(tostring(self)),
        class_<DummySequence, Sequence, SequencePtr>
        ("DummySequence")
        .def(constructor<>())
        .def(constructor<char>())
        .def(constructor<char, int>())
        .def("letter", &DummySequence::letter)
        .def("set_letter", &DummySequence::set_letter),
        class_<FragmentSequence, Sequence, SequencePtr>
        ("FragmentSequence")
        .def(constructor<>())
        .def(constructor<Fragment*>())
        .def("fragment", &FragmentSequence::fragment)
        .def("set_fragment",
             &FragmentSequence::set_fragment),
        class_<Fragment>("Fragment")
        .scope [
            def("new", &new_fragment0),
            def("new", &new_fragment1),
            def("new", &new_fragment2),
            def("new", &new_fragment3),
            def("delete", &delete_fragment)
        ]
        .def("seq", &Fragment::seq)
        .def("block", &Fragment::block)
        .def("min_pos", &Fragment::min_pos)
        .def("set_min_pos", &Fragment::set_min_pos)
        .def("max_pos", &Fragment::max_pos)
        .def("set_max_pos", &Fragment::set_max_pos)
        .def("ori", &Fragment::ori)
        .def("set_ori", &Fragment::set_ori)
        .def("set_ori", &fragment_set_ori)
        .def("length", &Fragment::length)
        .def("alignment_length",
             &Fragment::alignment_length)
        .def("inverse", &Fragment::inverse)
        .def("inverse", &fragment_inverse)
        .def("begin_pos", &Fragment::begin_pos)
        .def("set_begin_pos", &Fragment::set_begin_pos)
        .def("last_pos", &Fragment::last_pos)
        .def("set_last_pos", &Fragment::set_last_pos)
        .def("set_begin_last", &Fragment::set_begin_last)
        .def("end_pos", &Fragment::end_pos)
        .def("str", &Fragment::str)
        .def("str", &fragment_str)
        .def("substr", &Fragment::substr)
        .def("subfragment", &Fragment::subfragment)
        .def("clone", &Fragment::clone)
        .def("id", &Fragment::id)
        .def("hash", &Fragment::hash)
        .def("valid", &Fragment::valid)
        .def("has", &Fragment::has)
        .def("raw_at", &Fragment::raw_at)
        .def("at", &Fragment::at)
        .def("alignment_at", &Fragment::alignment_at)
        .def("common_positions",
             &Fragment::common_positions)
        .def("dist_to", &Fragment::dist_to)
        .def("is_subfragment_of",
             &Fragment::is_subfragment_of)
        .def("is_internal_subfragment_of",
             &Fragment::is_internal_subfragment_of)
        .def("row", &Fragment::row)
        .def("detach_row", &Fragment::detach_row)
        .def("set_row", &Fragment::set_row)
        .def("header", &fragment_header)
        .def("header", &fragment_header2)
        .def("contents", &fragment_contents)
        .def("contents", &fragment_contents2)
        .def(tostring(self))
        .def(const_self == const_self)
        .def(const_self < const_self)
        .def("block_pos", &block_pos)
        .def("fragment_pos", &block_pos)
        .def("frag_to_seq", &frag_to_seq)
        .def("seq_to_frag", &seq_to_frag),
        class_<AlignmentRow>("AlignmentRow")
        .def("clear", &AlignmentRow::clear)
        .def("grow", &AlignmentRow::grow)
        .def("bind", &AlignmentRow::bind)
        .def("map_to_alignment",
             &AlignmentRow::map_to_alignment)
        .def("map_to_fragment",
             &AlignmentRow::map_to_fragment)
        .def("length", &AlignmentRow::length)
        .def("set_length", &AlignmentRow::set_length)
        .def("fragment", &AlignmentRow::fragment)
        .def("nearest_in_fragment",
             &AlignmentRow::nearest_in_fragment)
        .def("clone", &AlignmentRow::clone)
        .def("slice", &AlignmentRow::slice)
        .def("type", &AlignmentRow::type)
        .enum_("Type") [
            value("MAP_ROW", MAP_ROW),
            value("COMPACT_ROW", COMPACT_ROW)
        ]
        .scope [
            def("new", &AlignmentRow::new_row),
            def("new", &alignmentrow_new_row),
            def("delete", &alignmentrow_delete_row)
        ],
        class_<Block>("Block")
        .scope [
            def("new", &new_block0),
            def("new", &new_block1),
            def("delete", &delete_block)
        ]
        .def("insert", &Block::insert)
        .def("erase", &Block::erase)
        .def("detach", &Block::detach)
        .def("size", &Block::size)
        .def("empty", &Block::empty)
        .def("has", &Block::has)
        .def("clear", &Block::clear)
        .def("swap", &Block::swap)
        .def("front", &Block::front)
        .def("fragments", &block_fragments,
             return_stl_iterator)
        .def("alignment_length", &Block::alignment_length)
        .def("identity", &Block::identity)
        .def("consensus_char", &block_consensus_char)
        .def("consensus_string", &block_consensus_string)
        .def("match", &Block::match)
        .def("inverse", &Block::inverse)
        .def("inverse", &block_inverse)
        .def("slice", &Block::slice)
        .def("clone", &Block::clone)
        .def("remove_alignment", &Block::remove_alignment)
        .def("max_shift_end", &Block::max_shift_end)
        .def("common_positions", &Block::common_positions)
        .def("merge", &Block::merge)
        .def("name", &Block::name)
        .def("set_name", &Block::set_name)
        .def("set_canonical_name", &set_canonical_name)
        .def("set_random_name", &Block::set_random_name)
        .def("set_name_from_fragments",
             &Block::set_name_from_fragments)
        .def("weak", &Block::weak)
        .def("set_weak", &Block::set_weak)
        .def(tostring(self))
        .def("is_ident_nogap", &is_ident_nogap)
        .def("test_column", &lite_test_column,
             pure_out_value(_3) + pure_out_value(_4))
        .def("hash", &block_hash_str)
        .def("id", &block_id)
        .def("has_repeats", &has_repeats)
        .def("is_exact_stem", &is_exact_stem)
        .def("make_name", &block_name)
        .def("has_alignment", &has_alignment)
        .def("test_block", &test_block)
        .def("find_slice", &find_slice,
             pure_out_value(_1) + pure_out_value(_2)),
        class_<AlignmentStat>("AlignmentStat")
        .def(constructor<>())
        .def("ident_nogap", &AlignmentStat::ident_nogap)
        .def("ident_gap", &AlignmentStat::ident_gap)
        .def("noident_nogap", &AlignmentStat::noident_nogap)
        .def("noident_gap", &AlignmentStat::noident_gap)
        .def("pure_gap", &AlignmentStat::pure_gap)
        .def("total", &AlignmentStat::total)
        .def("spreading", &AlignmentStat::spreading)
        .def("alignment_rows",
             &AlignmentStat::alignment_rows)
        .def("min_fragment_length",
             &AlignmentStat::min_fragment_length)
        .def("overlapping_fragments",
             &AlignmentStat::overlapping_fragments)
        .def("letter_count", &alignmentstat_letter_count)
        .def("gc", &AlignmentStat::gc)
        .def("make_stat", &make_stat)
        .def("make_stat", &make_stat0)
        .def("block_identity", &alignmentstat_identity),
        def("block_identity", &block_identity0),
        def("strict_block_identity", &strict_block_identity),
        class_<Sequences>("Sequences")
        .def("iter", &sequences_iter, return_stl_iterator),
        class_<BlockSet, BlockSetPtr>("BlockSet")
        .def(constructor<>())
        .def("add_sequence", &BlockSet::add_sequence)
        .def("seqs", &BlockSet::seqs)
        .def("remove_sequence", &BlockSet::remove_sequence)
        .def("insert", &BlockSet::insert)
        .def("erase", &BlockSet::erase)
        .def("detach", &BlockSet::detach)
        .def("size", &BlockSet::size)
        .def("empty", &BlockSet::empty)
        .def("has", &BlockSet::has)
        .def("clear_blocks", &BlockSet::clear_blocks)
        .def("clear_seqs", &BlockSet::clear_seqs)
        .def("clear", &BlockSet::clear)
        .def("swap", &BlockSet::swap)
        .def("clone", &BlockSet::clone)
        .def("blocks", &blockset_blocks,
             return_stl_iterator)
        .def(tostring(self))
        .def("read", &blockset_read)
        .def("genomes_number", &genomes_number)
        .def("hash", &blockset_hash1)
        .def("bsa", &blockset_bsa)
        .def("bsas", &BlockSet::bsas)
        .def("has_bsa", &BlockSet::has_bsa)
        .def("remove_bsa", &BlockSet::remove_bsa)
        .def("clear_bsas", &BlockSet::clear_bsas),
        class_<Fragments>("Fragments")
        .def(constructor<>())
        .def("iter", &fragments_iter, return_stl_iterator)
        .def("empty", &Fragments::empty)
        .def("clear", &Fragments::clear)
        .def("size", &Fragments::size)
        .def("push_back", &Fragments::push_back)
        .def("at", &fragments_at),
        class_<Blocks>("Blocks")
        .def(constructor<>())
        .def("iter", &blocks_iter, return_stl_iterator)
        .def("empty", &Blocks::empty)
        .def("clear", &Blocks::clear)
        .def("size", &Blocks::size)
        .def("push_back", &Blocks::push_back)
        .def("at", &blocks_at),
        class_<BSRow>("BSRow")
        .def(constructor<>())
        .def_readwrite("ori", &BSRow::ori)
        .def_readwrite("fragments", &BSRow::fragments),
        class_<BSA>("BSA")
        .def(constructor<>())
        .def("empty", &BSA::empty)
        .def("clear", &BSA::clear)
        .def("size", &BSA::size)
        .def("has", &bsa_has)
        .def("get", &bsa_get)
        .def("erase", &bsa_erase)
        .def("length", &bsa_length)
        .def("is_circular", &bsa_is_circular)
    ];
    return 0;
}

