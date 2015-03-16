#ifndef FERMAT_H
#define FERMAT_H

#include <vector>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/utils/point.hpp"

/*
 * Namespace Utils.
 * Contains utility functions for the heuristic ESMT algorithm
 */
namespace Utils {
  
  /**
   * Finds the Fermat-Torricelli point of a triangle
   *
   * @param A        1st point of triangle
   * @param B        2nd point of triangle
   * @param C        3rd point of triangle
   * @param res      Point to store the result in.
   *
   * @return         -1 if a Steiner Point is added, otherwise the index of the
   *                 point, which the SP is incident with (A, B or C)
   *                 is returned. The Fermat-Torricelli point is stored in res.
   */
  int getFermatPoint(Point &A, Point &B, Point &C, Point &res);

  /**
   * Similar to getFermatPoint, but an allocated Steiner Tree is returned instead.
   *
   * @param triangle The triangle stored in a Graph (size must be 3)
   * @param add_sp   If true, the sp will be added even if the sp is coincident
   *                 with a terminal.
   *
   * @return         An allocated Steiner Tree.
   */
  SteinerTree* getFermatSMT(Graph &triangle, bool add_sp = false);
}

#endif // FERMAT_H
