#ifndef STEINER_FINDER_H
#define STEINER_FINDER_H

#include <vector>

#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/utils/point.hpp"

/**
 * @class SteinerFinder
 *
 * This class uses the concept of Simplex Partitioning. Fermat-Torricelli is used
 * for solving 3-point instances.
 */
class SteinerFinder : public SubgraphHeuristic {
public:
  /**
   * Constructor
   */
  SteinerFinder(SubgraphHeuristic* sh);
  /**
   * Destructor
   */
  ~SteinerFinder();

  /**
   * Finds an approximated ESMT for the given graph.
   *
   * @param subgraph   Input subgraph.
   *
   * @return           Approximated ESMT for the given point set.
   */
  SteinerTree* findSteinerPoints(Graph &subgraph);
  
protected:
private:
  SubgraphHeuristic* sh;
  /**
   * dimension
   */
  int dim;

  double best_length;
  double best_ratio;
  SteinerTree *best_tree;
  double mst_length;

  void simplexPartitionRec(std::vector<Utils::Point> &points,
			   std::vector<int> &curPart, unsigned int index);

  SteinerTree *merge(std::vector<Utils::Point> &points, std::vector<int> &partOne);

  /**
   * Find the centroid of a point set
   * 
   * @param subgraph   The point set
   *
   * @return           The centroid of the point set
   */
  Utils::Point centroid(Graph& subgraph);
};

#endif // STEINER_FINDER_H
