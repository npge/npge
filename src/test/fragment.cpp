/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "read_block_set.hpp"

BOOST_AUTO_TEST_CASE (Fragment_main) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 9, 1);
    BOOST_REQUIRE(f1.length() == 10);
    BOOST_CHECK(f1.str() == "TGGTCCGAGA");
    Fragment f2(s1, 0, 9, -1);
    BOOST_REQUIRE(f2.length() == 10);
    BOOST_CHECK(f2.str() == "TCTCGGACCA");
    BOOST_CHECK(f1.substr(1, 1) == "G");
    BOOST_CHECK(f1.substr(1, 2) == "GG");
    BOOST_CHECK(f1.substr(1, -1) == "GGTCCGAGA");
    BOOST_CHECK(f1.substr(-2, -1) == "GA");
    BOOST_CHECK(f2.substr(1, 1) == "C");
    BOOST_CHECK(f2.substr(1, 2) == "CT");
    BOOST_CHECK(f2.substr(1, -1) == "CTCGGACCA");
    BOOST_CHECK(f2.substr(-2, -1) == "CA");
}

BOOST_AUTO_TEST_CASE (Fragment_begin_last) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 9, 1);
    BOOST_REQUIRE(f1.length() == 10);
    f1.set_begin_pos(5);
    BOOST_CHECK(f1 == Fragment(s1, 5, 9));
    f1.set_begin_pos(9);
    BOOST_CHECK(f1 == Fragment(s1, 9, 9));
    f1.set_last_pos(9);
    BOOST_CHECK(f1 == Fragment(s1, 9, 9));
    f1.set_last_pos(11);
    BOOST_CHECK(f1 == Fragment(s1, 9, 11));
    f1.inverse();
    BOOST_CHECK(f1 == Fragment(s1, 9, 11, -1));
    f1.set_begin_pos(10);
    BOOST_CHECK(f1 == Fragment(s1, 9, 10, -1));
    f1.set_begin_pos(15);
    BOOST_CHECK(f1 == Fragment(s1, 9, 15, -1));
    f1.set_last_pos(0);
    BOOST_CHECK(f1 == Fragment(s1, 0, 15, -1));
    f1.set_last_pos(15);
    BOOST_CHECK(f1 == Fragment(s1, 15, 15, -1));
    f1.set_begin_last(0, 0);
    BOOST_CHECK(f1 == Fragment(s1, 0, 0, 1));
    f1.set_begin_last(0, 1);
    BOOST_CHECK(f1 == Fragment(s1, 0, 1, 1));
    f1.set_begin_last(5, 1);
    BOOST_CHECK(f1 == Fragment(s1, 1, 5, -1));
    f1.set_begin_last(100, 0);
    BOOST_CHECK(f1 == Fragment(s1, 0, 100, -1));
}

BOOST_AUTO_TEST_CASE (Fragment_subfragment) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    // ----------------------------------------------------0123456789
    Fragment f1(s1, 0, 9, 1);
    Fragment f2(s1, 0, 9, -1);
    Fragment* tmp;
    //
    tmp = f1.subfragment(1, 1);
    BOOST_CHECK(tmp->str() == "G");
    delete tmp;
    //
    tmp = f1.subfragment(0, 5);
    BOOST_CHECK(tmp->str() == "TGGTCC");
    delete tmp;
    //
    tmp = f1.subfragment(5, 0);
    BOOST_CHECK(tmp->str() == "GGACCA");
    delete tmp;
    //
    tmp = f2.subfragment(1, 1);
    BOOST_CHECK(tmp->str() == "C");
    delete tmp;
    //
    tmp = f2.subfragment(0, 5);
    BOOST_CHECK(tmp->str() == "TCTCGG");
    delete tmp;
    //
    tmp = f2.subfragment(5, 0);
    BOOST_CHECK(tmp->str() == "CCGAGA");
    delete tmp;
}

BOOST_AUTO_TEST_CASE (Fragment_assign) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 9, 1);
    Fragment f2(f1);
    f2 = f2;
    Fragment f3;
    f3 = f2;
    BOOST_CHECK(f1 == Fragment(s1, 0, 9, 1));
    BOOST_CHECK(f2 == Fragment(s1, 0, 9, 1));
    BOOST_CHECK(f3 == Fragment(s1, 0, 9, 1));
}

BOOST_AUTO_TEST_CASE (Fragment_equal) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 9, 1);
    Fragment f2(s1, 0, 9, 1);
    Fragment f3(s1, 0, 9, -1);
    BOOST_CHECK(f1 == f1);
    BOOST_CHECK(f1 == f2);
    BOOST_CHECK(f1 != f3);
}

BOOST_AUTO_TEST_CASE (Fragment_less) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    BOOST_CHECK(!(Fragment(s1, 0, 9, 1) < Fragment(s1, 0, 9, 1)));
    BOOST_CHECK(Fragment(s1, 0, 9, 1) < Fragment(s1, 2, 9, 1));
    BOOST_CHECK(Fragment(s1, 0, 9, 1) < Fragment(s1, 0, 10, 1));
    BOOST_CHECK(Fragment(s1, 0, 9, -1) < Fragment(s1, 0, 9, 1));
    BOOST_CHECK(Fragment(s1, 0, 9, -1) < Fragment(s2, 0, 9, 1) ||
                Fragment(s2, 0, 9, -1) < Fragment(s1, 0, 9, 1));
}

BOOST_AUTO_TEST_CASE (Fragment_raw_at) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment f1(s1, 5, 10, 1);
    Fragment f2(s1, 5, 10, -1);
    BOOST_CHECK(f1.raw_at(0) == 'C');
    BOOST_CHECK(f1.raw_at(1) == 'G');
    BOOST_CHECK(f1.raw_at(-1) == 'C');
    BOOST_CHECK(f1.raw_at(-2) == 'T');
    BOOST_CHECK(f2.raw_at(0) == 'A');
    BOOST_CHECK(f2.raw_at(1) == 'T');
    BOOST_CHECK(f2.raw_at(-1) == 'C');
    BOOST_CHECK(f2.raw_at(-2) == 'G');
}

BOOST_AUTO_TEST_CASE (Fragment_at) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment f1(s1, 5, 10, 1);
    Fragment f2(s1, 5, 10, -1);
    BOOST_CHECK(f1.at(0) == 'C');
    BOOST_CHECK(f1.at(1) == 'G');
    BOOST_CHECK(f1.at(-1) == 'T');
    BOOST_CHECK(f1.at(-2) == 'A');
    BOOST_CHECK(f2.at(0) == 'A');
    BOOST_CHECK(f2.at(1) == 'T');
    BOOST_CHECK(f2.at(-1) == 'G');
    BOOST_CHECK(f2.at(-2) == 'C');
}

BOOST_AUTO_TEST_CASE (Fragment_alignment_at) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTC");
    Fragment f(s1, 0, 4);
    BOOST_CHECK(f.alignment_at(-2) == 0);
    BOOST_CHECK(f.alignment_at(-1) == 0);
    BOOST_CHECK(f.alignment_at(0) == 'T');
    BOOST_CHECK(f.alignment_at(1) == 'G');
    BOOST_CHECK(f.alignment_at(4) == 'C');
    BOOST_CHECK(f.alignment_at(5) == 0);
    f.set_row(new MapAlignmentRow("T--GGT-C"));
    BOOST_CHECK(f.alignment_at(-2) == 0);
    BOOST_CHECK(f.alignment_at(-1) == 0);
    BOOST_CHECK(f.alignment_at(0) == 'T');
    BOOST_CHECK(f.alignment_at(1) == 0);
    BOOST_CHECK(f.alignment_at(4) == 'G');
    BOOST_CHECK(f.alignment_at(5) == 'T');
    BOOST_CHECK(f.alignment_at(6) == 0);
    BOOST_CHECK(f.alignment_at(7) == 'C');
    BOOST_CHECK(f.alignment_at(8) == 0);
}

BOOST_AUTO_TEST_CASE (Fragment_next) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Fragment::connect(f1, f2);
    Fragment::connect(f2, f3);
    Fragment::connect(f3, f1);
    BOOST_REQUIRE(f1->next() == f2);
    BOOST_REQUIRE(f2->next() == f3);
    BOOST_REQUIRE(f3->next() == f1);
    BOOST_REQUIRE(f1->prev() == f3);
    BOOST_REQUIRE(f2->prev() == f1);
    BOOST_REQUIRE(f3->prev() == f2);
    f2->disconnect(/* connect_neighbors */ true);
    BOOST_CHECK(f1->next() == f3);
    BOOST_CHECK(f1->prev() == f3);
    BOOST_CHECK(f3->prev() == f1);
    BOOST_CHECK(f3->next() == f1);
    f1->disconnect(/* connect_neighbors */ false);
    BOOST_CHECK(!f1->prev());
    BOOST_CHECK(!f1->next());
    BOOST_CHECK(!f3->prev());
    BOOST_CHECK(!f3->next());
    delete f1;
    delete f2;
    delete f3;
}

BOOST_AUTO_TEST_CASE (Fragment_dtor) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Fragment::connect(f1, f2);
    Fragment::connect(f2, f3);
    Fragment::connect(f3, f1);
    BOOST_REQUIRE(f1->next() == f2);
    BOOST_REQUIRE(f2->next() == f3);
    BOOST_REQUIRE(f3->next() == f1);
    BOOST_REQUIRE(f1->prev() == f3);
    BOOST_REQUIRE(f2->prev() == f1);
    BOOST_REQUIRE(f3->prev() == f2);
    delete f2; // (= disconnect(true) )
    BOOST_CHECK(f1->next() == f3);
    BOOST_CHECK(f1->prev() == f3);
    BOOST_CHECK(f3->prev() == f1);
    BOOST_CHECK(f3->next() == f1);
    f1->disconnect(/* connect_neighbors */ false);
    BOOST_CHECK(!f1->prev());
    BOOST_CHECK(!f1->next());
    BOOST_CHECK(!f3->prev());
    BOOST_CHECK(!f3->next());
    delete f1;
    delete f3;
}

BOOST_AUTO_TEST_CASE (Fragment_connect_ori) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment::connect(f1, f2);
    BOOST_REQUIRE(f1->next() == f2);
    BOOST_REQUIRE(!f2->next());
    BOOST_REQUIRE(!f1->prev());
    BOOST_REQUIRE(f2->prev() == f1);
    Fragment::connect(f1, f2, -1); // cycle
    BOOST_REQUIRE(f1->next() == f2);
    BOOST_REQUIRE(f2->next() == f1);
    BOOST_REQUIRE(f1->prev() == f2);
    BOOST_REQUIRE(f2->prev() == f1);
    delete f1;
    delete f2;
}

BOOST_AUTO_TEST_CASE (Fragment_rearrange_with) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtcCGAGatgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Fragment::connect(f1, f2, -1); // wrong order
    Fragment::connect(f2, f3, -1); // wrong order
    f1->rearrange_with(f3);
    BOOST_CHECK(f1->next() == f2);
    BOOST_CHECK(f2->next() == f3);
    BOOST_CHECK(!f3->next());
    BOOST_CHECK(!f1->prev());
    BOOST_CHECK(f2->prev() == f1);
    BOOST_CHECK(f3->prev() == f2);
    f1->rearrange_with(f2);
    BOOST_CHECK(f1->next() == f3);
    BOOST_CHECK(f2->next() == f1);
    BOOST_CHECK(!f3->next());
    BOOST_CHECK(f1->prev() == f2);
    BOOST_CHECK(!f2->prev());
    BOOST_CHECK(f3->prev() == f1);
    delete f1;
    delete f2;
    delete f3;
}

BOOST_AUTO_TEST_CASE (Fragment_find_place) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtcCGAGatgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Fragment::connect(f1, f2, -1); // wrong order
    Fragment::connect(f2, f3, -1); // wrong order
    f1->find_place();
    f2->find_place();
    f3->find_place();
    BOOST_CHECK(f1->next() == f2);
    BOOST_CHECK(f2->next() == f3);
    BOOST_CHECK(!f3->next());
    BOOST_CHECK(!f1->prev());
    BOOST_CHECK(f2->prev() == f1);
    BOOST_CHECK(f3->prev() == f2);
    delete f1;
    delete f2;
    delete f3;
}

BOOST_AUTO_TEST_CASE (Fragment_find_place_f) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtcCGAGatgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Fragment::connect(f1, f3);
    f2->find_place(f1);
    BOOST_CHECK(f1->next() == f2);
    BOOST_CHECK(f2->next() == f3);
    BOOST_CHECK(!f3->next());
    BOOST_CHECK(!f1->prev());
    BOOST_CHECK(f2->prev() == f1);
    BOOST_CHECK(f3->prev() == f2);
    f1->find_place(f2);
    BOOST_CHECK(f1->next() == f2);
    BOOST_CHECK(f2->next() == f3);
    BOOST_CHECK(!f3->next());
    BOOST_CHECK(!f1->prev());
    BOOST_CHECK(f2->prev() == f1);
    BOOST_CHECK(f3->prev() == f2);
    delete f1;
    delete f2;
    delete f3;
}

BOOST_AUTO_TEST_CASE (Fragment_neighbor) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Fragment::connect(f1, f2);
    Fragment::connect(f2, f3);
    Fragment::connect(f3, f1);
    BOOST_CHECK(f1->neighbor(1) == f2);
    BOOST_CHECK(f1->neighbor(-1) == f3);
    BOOST_CHECK(f2->neighbor(1) == f3);
    BOOST_CHECK(f2->neighbor(-1) == f1);
    BOOST_CHECK(f3->neighbor(1) == f1);
    BOOST_CHECK(f3->neighbor(-1) == f2);
    //
    BOOST_CHECK(f1->logical_neighbor(1) == f2);
    BOOST_CHECK(f1->logical_neighbor(-1) == f3);
    BOOST_CHECK(f2->logical_neighbor(1) == f1);
    BOOST_CHECK(f2->logical_neighbor(-1) == f3);
    BOOST_CHECK(f3->logical_neighbor(1) == f1);
    BOOST_CHECK(f3->logical_neighbor(-1) == f2);
    //
    BOOST_CHECK(f1->is_neighbor(*f2));
    BOOST_CHECK(f2->is_neighbor(*f1));
    //
    BOOST_CHECK(f2->another_neighbor(*f1) == f3);
    BOOST_CHECK(f3->another_neighbor(*f2) == f1);
    delete f1;
    delete f2;
    delete f3;
}

BOOST_AUTO_TEST_CASE (Fragment_common_positions) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 5, 1);
    Fragment f2(s1, 5, 10, -1);
    Fragment f3(s1, 6, 8, -1);
    Fragment f4(s2, 6, 8, -1);
    BOOST_CHECK(f1.common_positions(f2) == 1);
    BOOST_CHECK(f2.common_positions(f1) == 1);
    BOOST_CHECK(f2.common_positions(f3) == 3);
    BOOST_CHECK(f3.common_positions(f2) == 3);
    BOOST_CHECK(f1.common_positions(f3) == 0);
    BOOST_CHECK(f3.common_positions(f1) == 0);
    BOOST_CHECK(f3.common_positions(f4) == 0);
}

BOOST_AUTO_TEST_CASE (Fragment_dist_to) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 5, 1);
    Fragment f2(s1, 5, 6, -1);
    Fragment f3(s1, 7, 8, -1);
    BOOST_CHECK(f1.dist_to(f2) == 0);
    BOOST_CHECK(f2.dist_to(f1) == 0);
    BOOST_CHECK(f1.dist_to(f3) == 1);
    BOOST_CHECK(f3.dist_to(f1) == 1);
    BOOST_CHECK(f2.dist_to(f3) == 0);
    BOOST_CHECK(f3.dist_to(f2) == 0);
}

BOOST_AUTO_TEST_CASE (Fragment_common_fragment) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 5, 1);
    Fragment f2(s1, 5, 10, -1);
    Fragment f3(s1, 6, 8, -1);
    Fragment f4(s2, 6, 8, -1);
    BOOST_CHECK(f1.common_fragment(f2) == Fragment(s1, 5, 5));
    BOOST_CHECK(f2.common_fragment(f1) == Fragment(s1, 5, 5, -1));
    BOOST_CHECK(f2.common_fragment(f3) == Fragment(s1, 6, 8, -1));
    BOOST_CHECK(f3.common_fragment(f2) == Fragment(s1, 6, 8, -1));
    BOOST_CHECK(!f1.common_fragment(f3).valid());
    BOOST_CHECK(!f3.common_fragment(f1).valid());
    BOOST_CHECK(!f3.common_fragment(f4).valid());
}

BOOST_AUTO_TEST_CASE (Fragment_is_subfragment) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 5, 1);
    Fragment f2(s1, 5, 10, -1);
    Fragment f3(s1, 6, 8, -1);
    Fragment f3a(s1, 5, 8, -1);
    Fragment f4(s2, 6, 8, -1);
    BOOST_CHECK(!f1.is_subfragment_of(f2));
    BOOST_CHECK(!f1.is_internal_subfragment_of(f2));
    BOOST_CHECK(f2.is_subfragment_of(f2));
    BOOST_CHECK(!f2.is_internal_subfragment_of(f2));
    BOOST_CHECK(!f2.is_subfragment_of(f1));
    BOOST_CHECK(!f2.is_internal_subfragment_of(f1));
    BOOST_CHECK(!f3.is_subfragment_of(f4));
    BOOST_CHECK(!f3.is_internal_subfragment_of(f4));
    BOOST_CHECK(f3.is_subfragment_of(f3a));
    BOOST_CHECK(!f3.is_internal_subfragment_of(f3a));
    BOOST_CHECK(f3.is_subfragment_of(f2));
    BOOST_CHECK(f3.is_internal_subfragment_of(f2));
}

BOOST_AUTO_TEST_CASE (Fragment_id) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccga|tggtccga");
    BOOST_CHECK(Fragment(s1, 1, 2).id() == "_1_2");
    BOOST_CHECK(Fragment(s1, 1, 2, -1).id() == "_2_1");
    s1->set_name("seq");
    BOOST_CHECK(Fragment(s1, 1, 2, -1).id() == "seq_2_1");
}

BOOST_AUTO_TEST_CASE (Fragment_print_contents) {
    using namespace npge;
    std::string seq(10, 'A');
    SequencePtr s1 = boost::make_shared<InMemorySequence>(seq);
    Fragment f(s1, 0, 9);
    {
        std::stringstream ss;
        f.print_contents(ss);
        BOOST_CHECK(ss.str() == seq);
    }
    {
        std::stringstream ss;
        f.print_contents(ss, '-', 5);
        BOOST_CHECK(ss.str() == "AAAAA\nAAAAA");
    }
    {
        std::stringstream ss;
        f.print_contents(ss, '-', 3);
        BOOST_CHECK(ss.str() == "AAA\nAAA\nAAA\nA");
    }
    {
        f.set_row(new MapAlignmentRow("AAAAAAAAAA"));
        std::stringstream ss;
        f.print_contents(ss, '-', 3);
        BOOST_CHECK(ss.str() == "AAA\nAAA\nAAA\nA");
    }
    {
        f.set_row(new MapAlignmentRow("A-A-A-A-A-A-A-A-A-A-"));
        std::stringstream ss;
        f.print_contents(ss, '-', 10);
        BOOST_CHECK(ss.str() == "A-A-A-A-A-\nA-A-A-A-A-");
    }
}

BOOST_AUTO_TEST_CASE (Fragment_reverse_of_length_1) {
    using namespace npge;
    std::string seq(10, 'A');
    SequencePtr s1;
    s1 = boost::make_shared<CompactLowNSequence>(seq);
    s1->set_name("a");
    Fragment f(s1, 0, 0);
    BOOST_CHECK(f.id() == "a_0_0");
    f.inverse();
    BOOST_CHECK(f.id() == "a_0_-1");
    typedef boost::scoped_ptr<Fragment> FragmentSc;
    FragmentSc f2((s1->fragment_from_id("a_0_0")));
    BOOST_CHECK(f2->ori() == 1);
    FragmentSc f3((s1->fragment_from_id("a_0_-1")));
    BOOST_CHECK(f3->ori() == -1);
}

BOOST_AUTO_TEST_CASE (Fragment_is_fragment_name) {
    using namespace npge;
    BOOST_CHECK(is_fragment_name("A_1_2"));
    BOOST_CHECK(is_fragment_name("A&d&c_1_2"));
    BOOST_CHECK(is_fragment_name("A&d&c_1_-1"));
    BOOST_CHECK(is_fragment_name("A&d&c_1000_500"));
    BOOST_CHECK(is_fragment_name("A&d&c_500_500"));
    BOOST_CHECK(is_fragment_name("A&d&c_500_100000"));
    BOOST_CHECK(!is_fragment_name("_500_100000"));
    BOOST_CHECK(!is_fragment_name("A_1.2_2"));
    BOOST_CHECK(!is_fragment_name("A_1_2.2"));
    BOOST_CHECK(!is_fragment_name("A__2"));
    BOOST_CHECK(!is_fragment_name("A___2"));
    BOOST_CHECK(!is_fragment_name("A___2_"));
    BOOST_CHECK(!is_fragment_name("A_1_2_"));
    BOOST_CHECK(!is_fragment_name("_aa_3_5"));
    BOOST_CHECK(!is_fragment_name("___5"));
    BOOST_CHECK(!is_fragment_name("___"));
    BOOST_CHECK(!is_fragment_name("a__2_3"));
}

