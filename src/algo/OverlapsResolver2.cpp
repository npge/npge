/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <ostream>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "OverlapsResolver2.hpp"
#include "Filter.hpp"
#include "Connector.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Exception.hpp"
#include "SortedVector.hpp"
#include "Graph.hpp"
#include "boundaries.hpp"
#include "convert_position.hpp"
#include "stick_impl.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

OverlapsResolver2::OverlapsResolver2(int min_distance):
    min_distance_(min_distance)
{ }

void OverlapsResolver2::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("min-distance", po::value<int>()->default_value(min_distance()),
     "Min distance between fragment boundaries")
   ;
}

void OverlapsResolver2::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("min-distance")) {
        int min_distance = vm["min-distance"].as<int>();
        if (min_distance < 0) {
            throw Exception("'min-distance' must be >= 0");
        }
        set_min_distance(min_distance);
    }
}

typedef std::pair<Sequence*, size_t> Point;
typedef std::pair<Point, Point> PointsPair;
typedef Graph<Point> PointsGraph;

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Point& point) {
    o << point.first->name() << " " << point.second;
    return o;
}

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const PointsPair& pair) {
    o << pair.first << " - " << pair.second;
    return o;
}

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const PointsGraph& graph) {
    BOOST_FOREACH (const PointsPair& points, graph) {
        o << points << std::endl;
    }
    return o;
}

/** Sum of size of all boundaries of Seq2Boundaries */
static size_t sb_size(const Seq2Boundaries& sb) {
    size_t result = 0;
    BOOST_FOREACH (const Seq2Boundaries::value_type& s_and_b, sb) {
        const Boundaries& b = s_and_b.second;
        result += b.size();
    }
    return result;
}

static void cat_boundaries(Seq2Boundaries& dest_sb,
                           const Seq2Boundaries& src_sb) {
    BOOST_FOREACH (const Seq2Boundaries::value_type& s_and_b, src_sb) {
        Sequence* seq = s_and_b.first;
        const Boundaries& src_b = s_and_b.second;
        Boundaries& dest_b = dest_sb[seq];
        dest_b.extend(src_b);
    }
}

static void sort_unique(Seq2Boundaries& sb) {
    BOOST_FOREACH (Seq2Boundaries::value_type& s_and_b, sb) {
        Boundaries& b = s_and_b.second;
        b.sort_unique();
    }
}

static bool is_sorted_unique(const Seq2Boundaries& sb) {
    BOOST_FOREACH (const Seq2Boundaries::value_type& s_and_b, sb) {
        const Boundaries& b = s_and_b.second;
        if (!b.is_sorted_unique()) {
            return false;
        }
    }
    return true;
}

static void stick_point(Point& point, const Seq2Boundaries& sb) {
    Sequence* seq = point.first;
    Seq2Boundaries::const_iterator it = sb.find(seq);
    BOOST_ASSERT(it != sb.end());
    const Boundaries& b = it->second;
    point.second = nearest_element(b, point.second);
}

// only for second points
static void stick_point_graph(PointsGraph& graph,
                              const Seq2Boundaries& boundaries) {
    BOOST_FOREACH (PointsPair& pair, graph) {
        stick_point(pair.second, boundaries);
    }
}

class HasNearest {
public:
    HasNearest(const Boundaries& boundaries, int min_distance):
        boundaries_(boundaries), min_distance_(min_distance)
    { }

    bool operator()(size_t boundary) {
        size_t new_pos = nearest_element(boundaries_, boundary);
        return std::abs(int(new_pos) - int(boundary)) < min_distance_;
    }

private:
    const Boundaries& boundaries_;
    size_t min_distance_;
};

static void filter_new_boundaries(Seq2Boundaries& new_sb,
                                  const Seq2Boundaries& old_sb,
                                  int min_distance) {
    BOOST_FOREACH (Seq2Boundaries::value_type& s_and_b, new_sb) {
        Sequence* seq = s_and_b.first;
        Boundaries& new_b = s_and_b.second;
        const Boundaries& old_b = old_sb.find(seq)->second;
        new_b.erase(std::remove_if(new_b.begin(), new_b.end(),
                                   HasNearest(old_b, min_distance)),
                    new_b.end());
    }
}

// /** Add new boundaries and graph edges.
// \param graph Graph where new edges are added to.
// \param new_b Empty vector where new points are added to.
// \param expand_b Points which are tried to be added.
// \param all_sb All added points, includes expand_b.
// \param bs BlockSet.
// \param min_distance Distance between points caused them to be sticked.
// */

int boundary_to_frag(const Fragment* f, size_t seq_pos) {
    int fr_pos = seq_to_frag(f, seq_pos);
    if (f->ori() == -1) {
        fr_pos += 1;
    }
    return fr_pos;
}

int frag_to_boundary(const Fragment* f, size_t fr_pos) {
    int seq_pos = frag_to_seq(f, fr_pos);
    if (f->ori() == -1) {
        seq_pos += 1;
    }
    return seq_pos;
}

/** Map selected point to each other fragment in the block, add to graph */
static void add_edges(PointsGraph& graph, const Block& block, int block_length,
                      const Fragment* from, size_t from_seq_pos) {
    int from_fr_pos = boundary_to_frag(from, from_seq_pos);
    BOOST_ASSERT(from_fr_pos >= 0 && from_fr_pos <= from->length());
    Point from_point(from->seq(), from_seq_pos);
    int block_p = block_pos(from, from_fr_pos, block_length);
    BOOST_FOREACH (const Fragment* to, block) {
        BOOST_ASSERT(to->length() > 0);
        if (to != from || block.size() == 1) { // for 1-blocks self-loops
            Sequence* to_seq = to->seq();
            int to_fr_pos = fragment_pos(to, block_p, block_length);
            size_t to_seq_pos = frag_to_boundary(to, to_fr_pos);
            Point to_point(to_seq, to_seq_pos);
            PointsPair pair(from_point, to_point);
            graph.push_back(pair);
        }
    }
}

/** Select appropriate points inside each fragment (from) */
static void add_edges(PointsGraph& graph, const Seq2Boundaries& expand_b,
                      const Block& block) {
    int block_length = block.alignment_length();
    BOOST_ASSERT(block_length > 0);
    BOOST_FOREACH (const Fragment* from, block) {
        BOOST_ASSERT(from->length() > 0);
        BOOST_ASSERT(from->length() <= block_length);
        Sequence* from_seq = from->seq();
        Seq2Boundaries::const_iterator expand_it = expand_b.find(from_seq);
        if (expand_it == expand_b.end()) {
            continue; // no new boundaries inside this fragment
        }
        const Boundaries& e_b = expand_it->second;
        Boundaries::const_iterator begin = e_b.lower_bound(from->min_pos());
        Boundaries::const_iterator end = e_b.upper_bound(from->max_pos() + 1);
        for (Boundaries::const_iterator it = begin; it != end; ++it) {
            size_t from_seq_pos = *it;
            add_edges(graph, block, block_length, from, from_seq_pos);
        }
    }
}

/** Add new boundaries.
\param graph Graph where new edges are added to.
\param expand_b Points which are tried to be added.
\param bs BlockSet.
*/
static void add_edges(PointsGraph& graph,
                      const Seq2Boundaries& expand_b,
                      const BlockSet& bs) {
    BOOST_FOREACH (const Block* block, bs) {
        add_edges(graph, expand_b, *block);
    }
}

// only for second points
static void extract_boundaries(Seq2Boundaries& boundaries,
                               const PointsGraph& graph) {
    BOOST_FOREACH (const PointsPair& pair, graph) {
        const Point& point = pair.second;
        Sequence* seq = point.first;
        size_t pos = point.second;
        boundaries[seq].push_back(pos);
    }
}

bool has_point(const Seq2Boundaries& sb, const Point& point) {
    Sequence* seq = point.first;
    Seq2Boundaries::const_iterator it = sb.find(seq);
    if (it == sb.end()) {
        return false;
    }
    const Boundaries& boundaries = it->second;
    size_t pos = point.second;
    return boundaries.has_elem(pos);
}

/** Return if all vertices from graph belong to sb
sb must be sorted.
*/
bool has_all_vertices(const PointsGraph& graph, const Seq2Boundaries& sb) {
    BOOST_FOREACH (const PointsPair& pair, graph) {
        if (!has_point(sb, pair.first) || !has_point(sb, pair.second)) {
            return false;
        }
    }
    return true;
}

static void build_point_graph(PointsGraph& graph, Seq2Boundaries& all_sb,
                              const BlockSetPtr& other, int min_distance) {
    BlockSet& bs = *other;
    Seq2Boundaries new_sb;
    bs_to_sb(new_sb, bs);
    stick_boundaries(new_sb, min_distance);
    stick_fragments(bs, new_sb, min_distance);
    remove_extra_sb(new_sb, bs);
    BOOST_ASSERT(sb_match_bs(new_sb, bs));
    Filter filter;
    filter.allow_everything();
    filter.set_min_fragment_length(1);
    filter.set_min_block_size(1);
    filter.apply(other);
    cat_boundaries(all_sb, new_sb); // new_sb is part of all_sb
    while (sb_size(new_sb) > 0) {
        PointsGraph new_g; // new edges of graph
        add_edges(new_g, new_sb, bs); // all mappings from new_sb
        Seq2Boundaries next_sb; // new points
        extract_boundaries(next_sb, new_g); // copy all destinations from new_g
        filter_new_boundaries(next_sb, all_sb, min_distance); // only new
        stick_boundaries(next_sb, min_distance); // only new
        BOOST_ASSERT(is_sorted_unique(next_sb));
        cat_boundaries(all_sb, next_sb); // append new points to all_sb
        size_t s1 = sb_size(all_sb);
        sort_unique(all_sb); // reorder all_sb
        BOOST_ASSERT(sb_size(all_sb) == s1);
        stick_point_graph(new_g, all_sb); // fix destinations in new_g
        BOOST_ASSERT(has_all_vertices(new_g, all_sb));
        graph.extend(new_g); // append new_g to graph
        new_sb.swap(next_sb); // new_sb = next_sb
    }
    graph.sort_unique();
    BOOST_ASSERT(has_all_vertices(graph, all_sb));
}

static Point neighbour_point(const Point& point, int ori,
                             const Seq2Boundaries& all_sb) {
    Sequence* seq = point.first;
    size_t pos = point.second;
    const Boundaries& b = all_sb.find(seq)->second;
    Boundaries::const_iterator it = b.lower_bound(pos);
    BOOST_ASSERT(it != b.end() && *it == pos);
    if (ori == 1) {
        ++it;
        return it == b.end() ? Point(0, 0) : Point(seq, *it);
    } else if (it == b.begin()) {
        return Point(0, 0); // no point
    } else {
        --it;
        return Point(seq, *it);
    }
}

struct MarkedFragment {
    Sequence* seq;
    size_t min_pos;
    size_t max_pos;

    /** Edge information.
    In edge.first flag is used for confirmation by input fragment.
    In edge.second flag is used as ori.
    */
    char flag;

    MarkedFragment(Sequence* s, size_t min, size_t max, char f = 0):
        seq(s), min_pos(min), max_pos(max), flag(f)
    { }

    bool operator<(const MarkedFragment& other) const {
        typedef boost::tuple<Sequence*, size_t> Tie;
        return Tie(seq, min_pos) < Tie(other.seq, other.min_pos);
    }

    bool operator==(const MarkedFragment& other) const {
        bool result = seq == other.seq && min_pos == other.min_pos;
        BOOST_ASSERT(!result || max_pos == other.max_pos);
        return result;
    }

    bool is_subfragment_of(const Fragment& f) const {
        return Fragment(seq, min_pos, max_pos).is_subfragment_of(f);
    }
};

static bool mgf_min(const MarkedFragment& a, const Fragment& b) {
    typedef boost::tuple<Sequence*, size_t> Tie;
    return Tie(a.seq, a.min_pos) < Tie(b.seq(), b.min_pos());
}

static bool mgf_min(const MarkedFragment& a, const MarkedFragment& b) {
    typedef boost::tuple<Sequence*, size_t> Tie;
    return Tie(a.seq, a.min_pos) < Tie(b.seq, b.min_pos);
}

static bool mgf_max(const Fragment& a, const MarkedFragment& b) {
    typedef boost::tuple<Sequence*, size_t> Tie;
    return Tie(a.seq(), a.max_pos()) < Tie(b.seq, b.max_pos);
}

typedef Graph<MarkedFragment> FragmentGraph;

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const MarkedFragment& mf) {
    o << mf.seq->name() << '_' << mf.min_pos << '_' << mf.max_pos;
    o << '(' << mf.flag << ')';
    return o;
}

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const FragmentGraph::Edge& edge) {
    o << edge.first << " - " << edge.second;
    return o;
}

/** Add edges to the graph of fragments */
static void build_fragment_graph(FragmentGraph& fg,
                                 const Seq2Boundaries& sb,
                                 const PointsGraph& pg) {
    BOOST_FOREACH (const Seq2Boundaries::value_type& s_and_b, sb) {
        Sequence* seq = s_and_b.first;
        const Boundaries& b = s_and_b.second;
        Boundaries::const_iterator min_pos_it, max_pos_it;
        min_pos_it = b.begin();
        max_pos_it = min_pos_it;
        ++max_pos_it;
        for (; max_pos_it < b.end(); ++min_pos_it, ++max_pos_it) {
            size_t min_pos = *min_pos_it;
            size_t max_pos = *max_pos_it;
            Point min_pos_point(seq, min_pos);
            Point max_pos_point(seq, max_pos);
            BOOST_ASSERT(min_pos < max_pos);
            MarkedFragment mf(seq, min_pos, max_pos - 1, 0);
            PointsGraph::Vertices min_friends, max_friends;
            pg.connected_with(min_friends, min_pos_point);
            pg.connected_with(max_friends, max_pos_point);
            BOOST_FOREACH (const Point& min_friend, min_friends) {
                Sequence* seq2 = min_friend.first;
                for (int ori = -1; ori <= 1; ori += 2) {
                    Point neighbour = neighbour_point(min_friend, ori, sb);
                    if (neighbour.first != 0) {
                        if (max_friends.has_elem(neighbour)) {
                            size_t f2_min_pos = (ori == 1) ?
                                                min_friend.second :
                                                neighbour.second;
                            size_t f2_max_pos = (ori == -1) ?
                                                min_friend.second :
                                                neighbour.second;
                            BOOST_ASSERT(f2_min_pos < f2_max_pos);
                            MarkedFragment mf2(seq2, f2_min_pos,
                                               f2_max_pos - 1, ori);
                            FragmentGraph::Edge fe(mf, mf2);
                            fg.push_back(fe);
                        }
                    }
                }
            }
        }
    }
}

/** Remove conflict edges (equal fragments, but different ori) */
void remove_conflict_edges(FragmentGraph& fg) {
    fg.sort();
    fg.remove_multiple();
}

static void add_block(BlockSet& bs,
                      const FragmentGraph::Vertices& marked_fragments,
                      const FragmentGraph& edges) {
    BOOST_ASSERT(!marked_fragments.empty());
    BOOST_ASSERT(edges.size() == marked_fragments.size() - 1);
    std::map<MarkedFragment, int> oris;
    const MarkedFragment* main = &edges.front().first;
    if (edges.empty()) {
        main = &marked_fragments.front();
    } else {
        main = &edges.front().first;
    }
    oris[*main] = 1;
    BOOST_FOREACH (const FragmentGraph::Edge& edge, edges) {
        const MarkedFragment& mf1 = edge.first;
        const MarkedFragment& mf2 = edge.second;
        BOOST_ASSERT(oris.find(mf1) != oris.end());
        BOOST_ASSERT(oris.find(mf2) == oris.end());
        int ori = mf2.flag;
        oris[mf2] = oris[mf1] * ori;
    }
    BOOST_ASSERT(oris.size() == marked_fragments.size());
    Block* block = new Block;
    BOOST_FOREACH (const MarkedFragment& mf, marked_fragments) {
        Fragment* new_f = new Fragment(mf.seq, mf.min_pos, mf.max_pos);
        BOOST_ASSERT(oris[mf]);
        new_f->set_ori(oris[mf]);
        block->insert(new_f);
    }
    bs.insert(block);
}

static void add_blocks(BlockSet& bs, const FragmentGraph& fg) {
    fg.connected_components(boost::bind(add_block, boost::ref(bs), _1, _2));
}

// TODO test all_sb min_distance

static void mark(FragmentGraph::Edge& edge) {
    edge.first.flag = 1;
}

static bool is_marked(const FragmentGraph::Edge& edge) {
    return edge.first.flag == 1;
}

typedef FragmentGraph::iterator FgIt;

struct CompareFirstBegin {
    bool operator()(const FragmentGraph::Edge& e,
                    const Fragment& src_f1) const {
        return mgf_min(e.first, src_f1);
    }

    bool operator()(const MarkedFragment& f1,
                    const FragmentGraph::Edge& e) const {
        return mgf_min(f1, e.first);
    }
};

struct CompareFirstEnd {
    bool operator()(const Fragment& src_f1,
                    const FragmentGraph::Edge& e) const {
        return mgf_max(src_f1, e.first);
    }
};

void find_internal_first(FgIt& begin, FgIt& end, FragmentGraph& g,
                         const Fragment& src_f1) {
    begin = std::lower_bound(g.begin(), g.end(), src_f1, CompareFirstBegin());
    end = std::upper_bound(g.begin(), g.end(), src_f1, CompareFirstEnd());
}

struct CompareSecondBegin {
    bool operator()(const FragmentGraph::Edge& e,
                    const Fragment& src_f2) const {
        return mgf_min(e.second, src_f2);
    }
};

struct CompareSecondEnd {
    bool operator()(const Fragment& src_f2,
                    const FragmentGraph::Edge& e) const {
        return mgf_max(src_f2, e.second);
    }
};

void find_internal_second(FgIt& begin2, FgIt& end2, FgIt begin, FgIt end,
                          const Fragment& src_f2) {
    begin2 = std::lower_bound(begin, end, src_f2, CompareSecondBegin());
    end2 = std::upper_bound(begin, end, src_f2, CompareSecondEnd());
}

static void mark_edges(FgIt begin, FgIt end, const Fragment& src_f2) {
    FgIt begin2, end2;
    find_internal_second(begin2, end2, begin, end, src_f2);
    for (FgIt it = begin2; it != end2; ++it) {
        FragmentGraph::Edge& edge = *it;
        BOOST_ASSERT(edge.second.is_subfragment_of(src_f2));
        mark(edge);
    }
}

static bool is_bad_edge(const FragmentGraph::Edge& edge) {
    return !is_marked(edge);
}

static void filter_fragment_graph(FragmentGraph& g, const BlockSet& bs) {
    BOOST_FOREACH (const Block* block, bs) {
        BOOST_FOREACH (const Fragment* f1, *block) {
            FgIt begin, end;
            find_internal_first(begin, end, g, *f1);
            for (FgIt it = begin; it != end; ++it) {
                FragmentGraph::Edge& edge = *it;
                BOOST_ASSERT(edge.first.is_subfragment_of(*f1));
            }
            FgIt begin1, end1; // with same edge.first
            for (begin1 = begin; begin1 != end; begin1 = end1) {
                end1 = std::upper_bound(begin1, end, begin1->first,
                                        CompareFirstBegin());
                BOOST_FOREACH (const Fragment* f2, *block) {
                    if (f1 != f2 || block->size() == 1 /* self-loops */) {
                        mark_edges(begin1, end1, *f2);
                    }
                }
            }
        }
    }
    g.erase(std::remove_if(g.begin(), g.end(), is_bad_edge), g.end());
}

bool OverlapsResolver2::run_impl() const {
    PointsGraph points_graph;
    Seq2Boundaries all_sb;
    build_point_graph(points_graph, all_sb, other(), min_distance());
    points_graph.add_for_symmetric();
    FragmentGraph fragment_graph;
    build_fragment_graph(fragment_graph, all_sb, points_graph);
    remove_conflict_edges(fragment_graph);
    BOOST_ASSERT(fragment_graph.is_symmetric());
    points_graph.clear();
    all_sb.clear();
    filter_fragment_graph(fragment_graph, *other());
    BOOST_ASSERT(fragment_graph.is_symmetric());
    block_set()->clear();
    add_blocks(*block_set(), fragment_graph);
#ifndef NDEBUG
    BOOST_ASSERT(!overlaps());
    Connector connector;
    connector.apply(block_set());
    BOOST_ASSERT(!overlaps());
#endif
    return true;
}

const char* OverlapsResolver2::name_impl() const {
    return "Resolve overlaping fragments (version 2)";
}

}

