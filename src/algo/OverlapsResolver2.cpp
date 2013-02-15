/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "OverlapsResolver2.hpp"
#include "Connector.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Exception.hpp"
#include "SortedVector.hpp"
#include "Graph.hpp"
#include "boundaries.hpp"
#include "convert_position.hpp"
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
    int min_distance = vm["min-distance"].as<int>();
    if (min_distance < 0) {
        throw Exception("'min-distance' must be >= 0");
    }
    set_min_distance(min_distance);
}

typedef std::pair<Sequence*, size_t> Point;
typedef SortedVector<Point> Points;
typedef std::pair<Point, Point> PointsPair;
typedef Graph<Point> PointsGraph;
typedef std::map<Sequence*, Boundaries> Seq2Boundaries;

static void cat_boundaries(Seq2Boundaries& dest_sb,
                           const Seq2Boundaries& src_sb) {
    BOOST_FOREACH (const Seq2Boundaries::value_type& s_and_b, src_sb) {
        Sequence* seq = s_and_b.first;
        const Boundaries& src_b = s_and_b.second;
        Boundaries& dest_b = dest_sb[seq];
        dest_b.insert(dest_b.end(), src_b.begin(), src_b.end());
    }
}

static void cat_point_graphs(PointsGraph& dest, const PointsGraph& src) {
    dest.insert(dest.end(), src.begin(), src.end());
}

static void stick_boundaries(Seq2Boundaries& sb, int min_distance) {
    BOOST_FOREACH (Seq2Boundaries::value_type& s_and_b, sb) {
        Boundaries& b = s_and_b.second;
        select_boundaries(b, min_distance);
    }
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
    graph.sort_unique();
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
                                  const Seq2Boundaries& old_sb, int min_distance) {
    BOOST_FOREACH (Seq2Boundaries::value_type& s_and_b, new_sb) {
        Sequence* seq = s_and_b.first;
        Boundaries& new_b = s_and_b.second;
        const Boundaries& old_b = old_sb.find(seq)->second;
        new_b.erase(std::remove_if(new_b.begin(), new_b.end(),
                                   HasNearest(old_b, min_distance)), new_b.end());
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

/** Map selected point to each other fragment in the block, add to graph */
static void add_edges(PointsGraph& graph, const Block& block, int block_length,
                      const Fragment* from, size_t from_seq_pos) {
    int from_fr_pos = seq_to_frag(from, from_seq_pos);
    if (from_fr_pos >= 0 && from_fr_pos <= from->length()) {
        Point from_point(from->seq(), from_seq_pos);
        int block_p = block_pos(from, from_fr_pos, block_length);
        BOOST_FOREACH (const Fragment* to, block) {
            if (to != from) {
                Sequence* to_seq = from->seq();
                int to_fr_pos = fragment_pos(to, block_p, block_length);
                size_t to_seq_pos = frag_to_seq(to, to_fr_pos);
                Point to_point(to_seq, to_seq_pos);
                PointsPair pair(from_point, to_point);
                graph.push_back(pair);
            }
        }
    }
}

/** Select appropriate points inside each fragment (from) */
static void add_edges(PointsGraph& graph, const Seq2Boundaries& expand_b,
                      const Block& block) {
    int block_length = block.alignment_length();
    BOOST_FOREACH (const Fragment* from, block) {
        Sequence* from_seq = from->seq();
        const Boundaries& e_b = expand_b.find(from_seq)->second;
        Boundaries::const_iterator begin = lower_bound(e_b, from->min_pos());
        Boundaries::const_iterator end = upper_bound(e_b, from->max_pos() + 1);
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

static void build_point_graph(PointsGraph& graph, Seq2Boundaries& all_sb,
                              const BlockSet& bs, int min_distance) {
    Seq2Boundaries new_sb;
    BOOST_FOREACH (const Block* block, bs) {
        BOOST_FOREACH (const Fragment* f, *block) {
            // new_sb contains all boundaries of initial fragments
            new_sb[f->seq()].push_back(f->min_pos());
            new_sb[f->seq()].push_back(f->max_pos() + 1);
        }
    }
    stick_boundaries(new_sb, min_distance);
    cat_boundaries(all_sb, new_sb); // new_sb is part of all_sb
    while (!new_sb.empty()) {
        PointsGraph new_g; // new edges of graph
        add_edges(new_g, new_sb, bs); // all mappings from new_sb
        Seq2Boundaries next_sb; // new points
        extract_boundaries(next_sb, new_g); // copy all destinations from new_g
        filter_new_boundaries(next_sb, all_sb, min_distance); // only new
        cat_boundaries(all_sb, next_sb); // append new points to all_sb
        stick_boundaries(all_sb, min_distance); // reorder all_sb
        stick_point_graph(new_g, all_sb); // fix destinations in new_g
        cat_point_graphs(graph, new_g); // append new_g to graph
        new_sb.swap(next_sb); // new_sb = next_sb
    }
    graph.sort_unique();
}

// TODO assert is_symmetric

static Point neighbour_point(const Point& point, int ori,
                             const Seq2Boundaries& all_sb) {
    Sequence* seq = point.first;
    size_t pos = point.second;
    const Boundaries& b = all_sb.find(seq)->second;
    Boundaries::const_iterator it = lower_bound(b, pos);
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

void mark(Fragment* f) {
    f->prev_ = reinterpret_cast<Fragment*>(-1);
}

bool is_marked(const Fragment* f) {
    return f->prev_ == reinterpret_cast<Fragment*>(-1);
}

typedef Graph<Fragment> FragmentGraph;

/** Add edges to the graph of fragments */
static void build_fragment_graph(FragmentGraph& fg,
                                 const Seq2Boundaries& sb, const PointsGraph& pg) {
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
            Fragment f(seq, min_pos, max_pos - 1);
            PointsGraph::Vertices min_friends, max_friends;
            pg.connected_with(min_friends, min_pos_point);
            pg.connected_with(max_friends, max_pos_point);
            BOOST_FOREACH (const Point& min_friend, min_friends) {
                Sequence* seq2 = min_friend.first;
                for (int ori = -1; ori <= 1; ori += 2) {
                    Point neighbour = neighbour_point(min_friend, ori, sb);
                    if (neighbour.first != 0) {
                        if (max_friends.has_elem(neighbour)) {
                            size_t f2_min_pos = ori == 1 ? min_friend.second : neighbour.second;
                            size_t f2_max_pos = ori == -1 ? min_friend.second : neighbour.second;
                            BOOST_ASSERT(f2_min_pos < f2_max_pos);
                            Fragment f2(seq2, f2_min_pos, f2_max_pos);
                            if (ori == -1) {
                                mark(&f2);
                            }
                            FragmentGraph::Edge fe(f, f2);
                            fg.push_back(fe);
                        }
                    }
                }
            }
        }
    }
    std::sort(fg.begin(), fg.end());
}

static void add_block(BlockSet& bs, const std::vector<Fragment>& fragments,
                      const FragmentGraph& edges) {
    std::map<Fragment, int> oris;
    const Fragment& main = edges.front().first;
    oris[main] = 1;
    BOOST_FOREACH (const FragmentGraph::Edge& edge, edges) {
        const Fragment& f1 = edge.first;
        const Fragment& f2 = edge.second;
        BOOST_ASSERT(oris.find(f1) != oris.end());
        BOOST_ASSERT(oris.find(f2) == oris.end());
        int ori = is_marked(&f2) ? -1 : 1;
        oris[f2] = oris[f1] * ori;
    }
    BOOST_ASSERT(oris.size() == fragments.size());
    Block* block = new Block;
    BOOST_FOREACH (const Fragment& f, fragments) {
        block->insert(new Fragment(f));
    }
    bs.insert(block);
}

static void add_blocks(BlockSet& bs, const FragmentGraph& fg) {
    fg.connected_components(boost::bind(add_block, boost::ref(bs), _1, _2));
}

// TODO preserve blocks of one fragment
// TODO test all_sb min_distance

bool OverlapsResolver2::run_impl() const {
    PointsGraph points_graph;
    Seq2Boundaries all_sb;
    build_point_graph(points_graph, all_sb, *other(), min_distance());
    FragmentGraph fragment_graph;
    build_fragment_graph(fragment_graph, all_sb, points_graph);
    points_graph.clear();
    all_sb.clear();
    add_blocks(*block_set(), fragment_graph);
#ifndef NDEBUG
    BOOST_ASSERT(!overlaps());
    Connector connector;
    connector.apply(block_set());
    BOOST_ASSERT(!overlaps());
#endif
    return true;
}

}

