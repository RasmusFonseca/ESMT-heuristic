#include <vector>
#include <iostream>

#include <steiner/graph.hpp>
#include <steiner/utils/utils.hpp>
#include <steiner/utils/point.hpp>

Graph::Graph() {
  this->mst_length = -1;
  this->points     = std::vector<Utils::Point>();
  this->edges      = std::vector<Utils::Edge>();
}

Graph::Graph(std::vector<Utils::Point> &points) {
  this->points     = points;
  this->edges      = std::vector<Utils::Edge>();
  this->mst_length = -1;
}

Graph::Graph(std::vector<Utils::Point> &points, std::vector<Utils::Edge> &edges) {
  this->mst_length = -1;
  this->points     = points;
  this->edges      = edges;
}

Graph::~Graph() {}

/*
 * Getter for the points
 */
std::vector<Utils::Point> Graph::getPoints() {
  return this->points;
}

/*
 * Getter for the edges
 */
std::vector<Utils::Edge> Graph::getEdges() {
  return this->edges;
}

/*
 * Getter for the terminals ptr
 */
std::vector<Utils::Point> *Graph::getPointsPtr() {
  return &this->points;
}

/*
 * Getter for the edges ptr
 */
std::vector<Utils::Edge> *Graph::getEdgesPtr() {
  return &this->edges;
}


double Graph::getMSTLength() {
  if(this->mst_length > 0)
    return this->mst_length;
  
  Graph mst = Utils::MSTKruskal(*this);
  this->mst_length = mst.getLength();
  
  return this->mst_length;
}

void Graph::setMSTLength(double l) {
  this->mst_length = l;
}

double Graph::getLength() {
  std::vector<Utils::Edge>::iterator it;
  double result = 0.0;
  for(it = this->edges.begin(); it != this->edges.end(); it++) {
    result += Utils::length(this->points[it->i0], this->points[it->i1]);
  }
  return result;
}

unsigned int Graph::dimension() {
  if(this->points.size() < 1)
    return 0;
  return this->points[0].dim();
}
