#ifndef SUBGRAPH_HEURISTIC_H
#define SUBGRAPH_HEURISTIC_H

#include <vector>

#include "steiner/steiner_tree.hpp"
#include "steiner/utils/point.hpp"

/**
 * @interface SubgraphHeuristic
 *
 * SubgraphHeuristic works as an interface for implementing subgraph
 * heuristics, thus giving a method independent API for finding Steiner
 * points in subgraphs.
 */
class SubgraphHeuristic {
public:
  /**
   * Constructor
   */
  SubgraphHeuristic() {}
  /**
   * Destructor
   */
  virtual ~SubgraphHeuristic() {}

  /**
   * Finds an approximated ESMT for the given graph, which must be
   * a MST.
   *
   * @param subgraph   The input subgraph. Must be an MST for the given
   *                   point set.
   */
  virtual SteinerTree *findSteinerPoints(Graph &subgraph) = 0;

protected:
private:
};

#endif
