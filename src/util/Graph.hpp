/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_GRAPH_HPP_
#define BR_GRAPH_HPP_

#include <set>
#include <algorithm>
#include <boost/foreach.hpp>

#include "global.hpp"
#include "SortedVector.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

template<typename Edge, typename Vertex>
struct EdgeCompare {
    bool operator()(const Edge& edge, const Vertex& vertex) const {
        return edge.first < vertex;
    }
};

template<typename V>
class GraphExtender {
public:
    GraphExtender(Graph<V>& graph):
        graph_(graph)
    { }

    void operator()(const SortedVector<V>& vertices,
                    const Graph<V>& edges) const;

private:
    Graph<V>& graph_;
};

/** A graph represented as sorted array of pairs */
template <typename V>
class Graph : public SortedVector<std::pair<V, V> > {
public:
    /** Type of edge (pair of vertices) */
    typedef typename std::pair<V, V> Edge;

    /** Typedef for vector */
    typedef SortedVector<Edge> BaseVector;

    /** Sorted vector of vertices */
    typedef SortedVector<V> Vertices;

    /** Constant iterator */
    typedef typename BaseVector::const_iterator const_iterator;

    /** Iterator */
    typedef typename BaseVector::iterator iterator;

    using BaseVector::begin;
    using BaseVector::end;
    using BaseVector::front;
    using BaseVector::back;
    using BaseVector::erase;
    using BaseVector::size;
    using BaseVector::push_back;
    using BaseVector::pop_back;
    using BaseVector::is_sorted_unique;
    using BaseVector::has_elem;
    using BaseVector::sort;
    using BaseVector::unique;
    using BaseVector::sort_unique;
    using BaseVector::extend;

    /** Return if the graph is symmetric */
    bool is_symmetric() const {
        BOOST_FOREACH (const Edge& e, *this) {
            Edge e_1(e.second, e.first);
            if (!has_elem(e_1)) {
                return false;
            }
        }
        return true;
    }

    /** Add lacking reverse edges.
    Input graph must be sorted, output graph is sorted and unique.
    */
    void add_for_symmetric() {
        Graph<V> new_edges;
        BOOST_FOREACH (const Edge& e, *this) {
            Edge e_1(e.second, e.first);
            if (!has_elem(e_1)) {
                new_edges.push_back(e_1);
            }
        }
        extend(new_edges);
        sort_unique();
        BOOST_ASSERT(is_symmetric());
    }

    /** Return iterator to first edge matching given vertex if it exists */
    const_iterator first_point_pair(const V& vertex) const {
        return std::lower_bound(begin(), end(), vertex,
                                EdgeCompare<Edge, V>());
    }

    /** Add all vertices, to which given vertex is connected */
    void connected_with(Vertices& vertices, const V& vertex) const {
        for (const_iterator it = first_point_pair(vertex);
                it != end() && it->first == vertex; ++it) {
            vertices.push_back(it->second);
        }
    }

    /** Add all vertices */
    void all_vertices(Vertices& vertices) const {
        BOOST_ASSERT(is_sorted_unique());
        BOOST_FOREACH (const Edge& edge, *this) {
            const V& v = edge.first;
            if (vertices.empty() || !(vertices.back() == v)) {
                vertices.push_back(v);
            }
        }
    }

    /** Apply a function to each connected component and vector of edges.
    Connected components may be not sorted.
    Order of edges is the order of vertex discovery.
    */
    template<typename F>
    void connected_components(const F& f) const {
        Vertices vertices;
        all_vertices(vertices);
        typedef std::set<V> VerticesSet;
        VerticesSet visited;
        BOOST_FOREACH (const V& vertex, vertices) {
            if (visited.find(vertex) == visited.end()) {
                visited.insert(vertex);
                Vertices component;
                Graph edges;
                Vertices new_vertices;
                new_vertices.push_back(vertex);
                while (!new_vertices.empty()) {
                    V v = new_vertices.back();
                    new_vertices.pop_back();
                    component.push_back(v);
                    Vertices connected;
                    connected_with(connected, v);
                    BOOST_FOREACH (const V& v2, connected) {
                        if (visited.find(v2) == visited.end()) {
                            visited.insert(v2);
                            edges.push_back(Edge(v, v2));
                            new_vertices.push_back(v2);
                        }
                    }
                }
                f(component, edges);
            }
        }
    }

    /** Remove edges preserving connected components and symmetric graph */
    void remove_extra_edges() {
        Graph<V> new_graph;
        connected_components(GraphExtender<V>(new_graph));
        sort();
        BOOST_ASSERT(is_sorted_unique());
        add_for_symmetric();
    }
};

template<typename V>
void GraphExtender<V>::operator()(const SortedVector<V>& /* vertices */,
                                  const Graph<V>& edges) const {
    graph_.extend(edges);
}

}

#endif

