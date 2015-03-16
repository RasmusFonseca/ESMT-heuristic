#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <tr1/unordered_map>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"

/**
 * Namespace Utils.
 * Contains utility functions for the heuristic ESMT algorithm
 */
namespace Utils {
  
  /**
   * Performs Kruskal's MST algorithm on the given graph.
   * Returns a list of the edges in the MST, which is a subset
   * of the edges in the given graph.
   *
   * @param graph   The graph to find the MST of.
   *
   * @return        The MST of graph.
   */
  Graph MSTKruskal(Graph &graph);

  /**
   * Gets a random number between 0 and 1
   *
   * @return    A random number between 0 and 1
   */
  double frand();

  /**
   * Constructs a bottleneck-graph in b_tree. Only edges from del will appear
   * in b_tree. Bottlenck-graph will be constructed from mst.
   *
   * @param b_graph  A hashmap, which will be filled with the bottleneck edges from
   *                 del.
   * @param del      A graph containing the edges to be added to the b-graph
   *                 (e.g. Delaunay). mst must be a subset of del.
   * @param mst      The MST for the gived del-graph.
   */
  void constructBottleneckGraph(std::tr1::unordered_map<long, double> &b_graph,
				Graph &del, Graph &mst);
  
  /**
   * Validates a SteinerTree
   *
   * @param st   The SteinerTree to be validated (|S|<|N|-1, connected, etc.).
   *
   * @return     True if valid, otherwise false (error will be printed
   *             to std::cerr).
   */
  bool validate(SteinerTree &st);
}

#endif // UTILS_H
