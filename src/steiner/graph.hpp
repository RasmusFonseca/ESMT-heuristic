#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

#include "steiner/utils/point.hpp"

class Graph {
public:

  /**
   * Default constructor.
   */
  Graph();
  
  /**
   * Constructor. Creates a graph with no edges.
   *
   * @param points   The points for this graph.
   */
  Graph(std::vector<Utils::Point> &points);

  /**
   * Constructor.
   *
   * @param points   The points for this graph.
   * @param edges    The edges for this graph.
   */
  Graph(std::vector<Utils::Point> &points, std::vector<Utils::Edge> &edges);
  
  /**
   * Destructor.
   */
  ~Graph();

  /**
   * Getter for the points.
   *
   * @return  A copy of the set of points for this graph.
   */
  std::vector<Utils::Point> getPoints();

  /**
   * Getter for the edges.
   *
   * @return  A copy of the set of edges for this graph.
   */
  std::vector<Utils::Edge> getEdges();

  /**
   * Getter for a pointer to the points.
   *
   * @return  A pointer to the list of points for this set.
   */
  std::vector<Utils::Point> *getPointsPtr();

  /**
   * Getter for a pointer to the edges.
   *
   * @return  A pointer to the list of edges for this set.
   */
  std::vector<Utils::Edge> *getEdgesPtr();

  /**
   * Calculates the length of the MST for this graph.
   *
   * @return   The length of the MST for this graph.
   */
  double getMSTLength();

  /**
   * Setter for the MST length.
   *
   * Used to set the MST length to avoid overhead of
   * calculating it twice.
   *
   * @param l   The new MST length.
   */
  void setMSTLength(double l);

  /**
   * Calculates the length of this graph.
   *
   * @return   The length of this graph.
   */
  double getLength();

  /**
   * Returns the dimension of this graph
   *
   * @return    The dimension of this graph
   */
  unsigned int dimension();
  
protected:
  
  /** The points of this graph. */
  std::vector<Utils::Point> points;
  /** The edges of this graph. */
  std::vector<Utils::Edge>  edges;
  /** The length of the MST of this graph. */
  double                    mst_length;
  
private:
};

#endif // GRAPH_H
