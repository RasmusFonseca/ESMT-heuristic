#ifndef ITERATIVE_H
#define ITERATIVE_H

#include <vector>

#include "steiner/graph.hpp"
#include "steiner/utils/point.hpp"

/**
 * @class Iterative
 *
 * Implements the optimisation method given by Smith in the paper:
 *
 *   W. D. Smith, How to find steiner minimal trees in euclidean d-space,
 *         Algorithmica 7 (1992), 137–177
 */
class Iterative {
public:
  /**
   * Constructor
   */
  Iterative(int dim, int n);
  
  /**
   * Destructor
   */
  virtual ~Iterative();

protected:
  void buildTree(int k, int *topo_vec);
  
  /**
   * Optimises the tree by solving equations.
   *
   * @param tol   A small number
   */
  void optimise(double tol);
  
  /**
   * Returns the length of the current tree and fills the EL attribute.
   *
   * @return   The length of the current tree.
   */
  double length();
  
  /**
   * Returns the error figure for the current tree
   *
   * @return   The error figure for the current tree (based on angles).
   */
  double error();
  
  /**
   * Calculates the euclidean distance between two points in the tree.
   *
   * @param i1   Index of first point.
   * @param i2   Index of second point.
   *
   * @return     The euclidean distance between points[i1] and points[i2]
   */
  double dist(int i1, int i2);

  /**
   * Cleans up the tree, e.g. removes SPs coincidents with terminals.
   * Fills points and edges.
   *
   * @param points    The vector to fill with points.
   * @param edges     The vector to fill with edges.
   * @param doCleanUp If set to true, Steiner points close to terminals will be
   *                  removed.
   */
  void cleanUp(std::vector<Utils::Point> *points, std::vector<Utils::Edge> *edges,
	       bool doCleanUp = true);

  std::vector<Utils::Point> P; // Points -> [0..N-1] is terminals, [N..S-1] steiner
  int **adj;                   // Adjacent for steiner points [0..S-1][3]
  int **edge;                  // Edges. Indexes into points  [2*N][2]
  double **EL;                 // Contains length of edges to SPs [N][3]
  unsigned int N;              // No of terminals
  unsigned int S;              // No of steiner points
  unsigned int max_dim;        // Max-dimension (allocated)
  unsigned int dim;            // Dimension
  int alloc;                   // Places allocated


private:

  /*
   * Subprocedure prep.
   */
  void prep(int i, int a, int b, double c);

  /*
   * Update references to a point
   */
  void updateRef(unsigned int i, int oldp, int newp,
		 std::vector< std::vector<int> > &tadj);

  // Private attribute, for use in optimise, etc.
  double **B;  // [N][3]
  double **C;  // [N][3]
  int *val;    // [N]

  // Stacks
  std::vector<int> leafQ;
  std::vector<int> eqnstack;
};

#endif // ITERATIVE_H
