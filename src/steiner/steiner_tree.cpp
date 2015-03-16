#include <vector>
#include <iostream>

#include "steiner/steiner_tree.hpp"
#include "steiner/utils/point.hpp"

SteinerTree::SteinerTree()
  : Graph()
{
  this->steiner_ratio  = -1.0;
  this->smt_length     = -1.0;
}

SteinerTree::SteinerTree(std::vector<Utils::Point> &terminals)
  : Graph(terminals)
{
  this->steiner_ratio  = -1.0;
  this->smt_length     = -1.0;
}

SteinerTree::SteinerTree(std::vector<Utils::Point> &terminals,
			 std::vector<Utils::Edge> &edges)
  : Graph(terminals, edges)
{
  this->steiner_ratio  = -1.0;
  this->smt_length     = -1.0;
}

SteinerTree::~SteinerTree() {}

double SteinerTree::getSMTLength() {
  if(this->smt_length > 0.0)
    return this->smt_length;
  
  // Assume that this is a SMT. Simply get the length
  this->smt_length = this->Graph::getLength();
  
  return this->smt_length;
}

void SteinerTree::setSMTLength(double l) {
  this->smt_length = l;
}

double SteinerTree::getSteinerRatio() {
  if(this->steiner_ratio > 0.0)
    return this->steiner_ratio;

  this->steiner_ratio = this->getSMTLength() / this->getMSTLength();
  return this->steiner_ratio;
}

void SteinerTree::setSteinerRatio(double l) {
  this->steiner_ratio = l;
}

std::ostream& operator<<(std::ostream& os, SteinerTree &st) {
  unsigned int i;
  std::string indent("  ");
  std::vector<Utils::Point> points = st.getPoints();
  std::vector<Utils::Edge>  edges  = st.getEdges();
  os << "# Terminals" << std::endl;
  for(i = 0; i < points.size(); i++) {
    if(points[i].isSteiner())
      break;
    os << indent << i << " " << points[i] << std::endl;
  }
  os << std::endl << "# Steiner points" << std::endl;
  for(; i < points.size(); i++) {
    os << indent << i << " " << points[i] << std::endl;
  }
  os << std::endl << "# Edges" << std::endl;
  for(i = 0; i < edges.size(); i++) {
    os << indent << "(" << edges[i].i0
       << " " << edges[i].i1 << ")" << std::endl;
  }
  os << std::endl << "# |MST|: " << st.getMSTLength() << std::endl
     << "# |SMT|: " << st.getSMTLength() << std::endl
     << "# Steiner ratio: " << st.getSteinerRatio() << std::endl;
  return os;
}
