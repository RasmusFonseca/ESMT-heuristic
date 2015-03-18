#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <unordered_map>

#include "steiner/graph.hpp"
#include "steiner/utils/delaunay.hpp"
#include "steiner/utils/point.hpp"

typedef Utils::Edge                Edge;
typedef Utils::Point               Point;
typedef Utils::Delaunay::Simplex   Simplex;

/* Implementation of the constructor */
Utils::Delaunay::Delaunay(std::vector<Utils::Point> &points) {
  this->points = points;
  this->doDelaunayQHull(this->points[0].dim());
}

/* Implementation of the destructor */
Utils::Delaunay::~Delaunay() {}

/* Implementation of getSimplices() */
std::vector<Utils::Delaunay::Simplex> *Utils::Delaunay::getSimplices() {
  return &this->simplices;
}

/* Implementation of getPointHandles() */
std::vector<Utils::Delaunay::PointHandle> *Utils::Delaunay::getPointHandles() {
  return &this->point_handles;
}

/* Implementation of doDelaunayQHull */
void Utils::Delaunay::doDelaunayQHull(int d) {
  doPrepInputQHull(d);
  int res = system("./steiner/utils/qhull/qdelaunay i < steiner/utils/qhull/in.txt Fn Qt > steiner/utils/qhull/out.txt");
  if(res) {
    std::cerr << "QHull error: " << res << std::endl;
    exit(1);
  }
  doParseOutputQHull(d);
}

/* Implementation of doPrepInputQHull */
void Utils::Delaunay::doPrepInputQHull(int d) {
  std::ofstream infile;
  int n = this->points.size();
  infile.open("steiner/utils/qhull/in.txt");
  infile << d << std::endl << n << std::endl;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < d; j++)
      infile << this->points[i][j] << " ";
    infile << std::endl;
  }
  infile.close();
}

/* Implementation of doParseOutputQHull */
void Utils::Delaunay::doParseOutputQHull(int d) {
  this->simplices.clear();
  this->point_handles.clear();
  this->point_handles.reserve(this->points.size());
  for(unsigned int i = 0; i < this->points.size(); i++)
    this->point_handles.push_back(PointHandle());
  
  std::ifstream outfile("steiner/utils/qhull/out.txt");
  std::string buffer;
  
  if(outfile.is_open()) {
    // Allocate needed resources:
    // - edge in MST flags
    std::unordered_map<unsigned long, bool> flag;
    // 1. Read all simplices
    int no_of_simplices;
    getline(outfile, buffer);
    // Get number of simplices
    std::istringstream(buffer) >> no_of_simplices;
    this->simplices.reserve(no_of_simplices);
    for(int i = 0; i < no_of_simplices; i++) {
      getline(outfile, buffer);
      std::istringstream iss(buffer);
      std::vector<std::string> tokens;
      std::copy(std::istream_iterator<std::string>(iss),
		std::istream_iterator<std::string>(),
		std::back_inserter< std::vector<std::string> >(tokens));
      if((int)tokens.size() != d+1) {
	std::cerr << "Bad result size from qdelaunay" << std::endl;
	exit(1);
      }
      Simplex s(d);
      for(int j = 0; j < d+1; j++)
	std::istringstream(tokens[j]) >> s.map[j];
      // Add edges and update pointhandles
      for(int j = 0; j < d+1; j++) {
	this->point_handles[s.map[j]].simplices.push_back(i);
	for(int k = j+1; k < d+1; k++) {
	  int i0 = s.map[j];
	  int i1 = s.map[k];
	  unsigned long key = Edge::key(i0,i1);
	  if(flag.find(key) == flag.end()) {
	    flag[key] = true;
	    this->edges.push_back(Edge(i0,i1, Utils::length(this->points[i0], this->points[i1])));
	  }
	}
      }
      this->simplices.push_back(s);
    }
    // 2. Read all neighbours
    // Remove first line with number of simplices
    getline(outfile, buffer);
    for(int i = 0; i < no_of_simplices; i++) {
      getline(outfile, buffer);
      std::istringstream iss(buffer);
      std::vector<std::string> tokens;
      std::copy(std::istream_iterator<std::string>(iss),
		std::istream_iterator<std::string>(),
		std::back_inserter< std::vector<std::string> >(tokens));
      if((int)tokens.size() != d+2) {
	std::cerr << "Bad result size from qdelaunay" << std::endl;
	exit(1);
      }
      for(int s = 0; s < d+1; s++) {
	this->simplices[i].nextSimplex[s] = -1;
	this->simplices[i].nextVertex[s]  = -1;  
      }
      for(int s = 0; s < d+1; s++) {
	int next_s = 0, next_v = 0, opposite = 0;
	std::istringstream(tokens[s+1]) >> next_s;
	if(next_s < 0)
	  continue;
	// Find the extra vertex in the neighbour
	for(int j = 0; j < d+1; j++) {
	  bool found = false;
	  for(int k = 0; k < d+1; k++) {
	    if(this->simplices[next_s].map[j] == this->simplices[i].map[k]) {
	      found = true;
	      break;
	    }
	  }
	  if(!found) {
	    // Found neighbour point
	    next_v = j;
	    break;
	  }
	}
	// Find the opposite vertex
	for(int j = 0; j < d+1; j++) {
	  bool found = false;
	  for(int k = 0; k < d+1; k++) {
	    if(this->simplices[i].map[j] == this->simplices[next_s].map[k]) {
	      found = true;
	      break;
	    }
	  }
	  if(!found) {
	    // Found neighbour point
	    opposite = j;
	    break;
	  }
	}
	// Set values
	this->simplices[i].nextSimplex[opposite] = next_s;
	this->simplices[i].nextVertex[opposite]  = this->simplices[next_s].map[next_v];
      }
    }
    outfile.close();
  }
}

//////////////////////////////
// Struct functions

/* Constructor */
Simplex::Simplex(int d)  {
  this->n   = d+1;
  this->map = new int[d+1];
  this->nextSimplex = new int[d+1];
  this->nextVertex  = new int[d+1];
}

/* Copy constructor */
Simplex::Simplex(const Simplex &s)  {
  this->n           = s.n;
  this->map         = new int[s.n];
  this->nextSimplex = new int[s.n];
  this->nextVertex  = new int[s.n];
  std::copy(&s.map[0], &s.map[0]+s.n, &this->map[0]);
  std::copy(&s.nextSimplex[0], &s.nextSimplex[0]+s.n, &this->nextSimplex[0]);
  std::copy(&s.nextVertex[0], &s.nextVertex[0]+s.n, &this->nextVertex[0]);
}

/* Destructor */
Simplex::~Simplex() {
  delete this->map;
  delete this->nextSimplex;
  delete this->nextVertex;
}    
