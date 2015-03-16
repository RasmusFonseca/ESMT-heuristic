#include <vector>
#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <cmath>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/utils.hpp"

typedef Utils::Point Point;
typedef Utils::Edge  Edge;

/**
 * Constructor
 */
IterativeSmith::IterativeSmith(int dim, bool doCleanUp)
  : SubgraphHeuristic(), Iterative(dim, 12), do_clean_up(doCleanUp) {}

/**
 * Destructor
 */
IterativeSmith::~IterativeSmith() {}

/**
 * Finds the steiner points and edges.
 */
SteinerTree *IterativeSmith::findSteinerPoints(Graph &subgraph) {
  unsigned int i;
  // Create result
  std::vector<Edge> edgelist;
  SteinerTree *result = new SteinerTree(*subgraph.getPointsPtr(),
					edgelist);
  // Get pointers to result lists
  std::vector<Point> *points = result->getPointsPtr();
  std::vector<Edge>  *edges  = result->getEdgesPtr();

  this->dim = subgraph.dimension();
  assert(this->dim < this->max_dim);
  
  // Check that the input is valid
  assert(points->size() > 2);
  assert(points->size() < NMAX);

  // Generate the MST as an upper bound
  double mst_length = subgraph.getMSTLength();
  
  assert(mst_length > 0.0);

  this->N = points->size();

  // Topology vector
  int topo_vec[NMAX];
  
  // Copy points
  this->P = subgraph.getPoints();
  // Add SP-dummies
  this->P.reserve(2*this->N-2);
  for(i = 0; i < this->N-2; i++)
    this->P.push_back(Point(this->dim));
  
  double l, r;

  if(this->N == 3) {
    // Special case
    this->buildTree(0, topo_vec);
    l = this->length();
    r = this->error();
    do {
      this->optimise(0.0001*r/this->N);
      l = this->length();
      r = this->error();
    }
    while(r>l*0.0001);
  }
  else {
    // General case - generate topologies
    int best_vec[NMAX];
    
    // Init topology vector
    for(i = 0; i < this->N - 3; i++) {
      topo_vec[i] = 0;
    }
    
    int k = 0;
    double upper_bound = mst_length;
    
    while(1) {
      // Build tree for given topology
      this->buildTree(k+1, topo_vec);
      
      // Optimise
      l = this->length();
      r = this->error();
      do {
	this->optimise(0.0001*r/this->N);
	l = this->length();
	r = this->error();
      }
      while(r>l*0.0001);
      
      if(l*0.999 <= upper_bound) {
	// This is a candidate ...
	if(k == (int)(this->N - 4)) {
	  // ... with all SPs. Save!
	  upper_bound = l;
	  for(i = 0; i < this->N - 3; i++)
	    best_vec[i] = topo_vec[i];
	}
	else {
	  // Simply continue
	  k++;
	  continue;
	}
      }
      // If above fails, there is no reason to consider any topo. vector
      // with this prefix. Increment to next prefix
      
      // Move back through the vector, if max has been reached
      while(k >= 0 && topo_vec[k] == 2*k + 2) {
	topo_vec[k] = 0;
	k--;
      }
      
      // All topologies considered?
      if(k == -1)
	// Done
	break;
      
      // Increment current level
      topo_vec[k]++;
    }
    
    // Rebuild the best tree
    this->buildTree(this->N - 3, best_vec);
    l = this->length();
    r = this->error();
    do {
      this->optimise(0.0001*r/this->N);
      l = this->length();
      r = this->error();
    }
    while(r>l*0.0001);  
  }
  // The tree is done

  this->cleanUp(points, edges, this->do_clean_up);
  
  // Compute ratio
  result->setMSTLength(mst_length);
  result->getSMTLength();
  result->getSteinerRatio();
  return result;
}

void IterativeSmith::buildTree(unsigned int l, int *topo_vec) {
  unsigned int k, d;
  // Assume: points[0..N-1] contains the N terminals
  assert(this->N > 2);
  
  // First, build a 3-terminal tree.
  int si = this->N; // Steiner index
  int m  = 0;

  for(d = 0; d < this->dim; d++) {
    this->P[si][d] = (this->P[0][d]
		      + this->P[1][d]
		      + this->P[2][d]) / 3.0 + 0.001 * Utils::frand();
  }
  
  // Set adjs
  this->adj[m][0] = 0;
  this->adj[m][1] = 1;
  this->adj[m][2] = 2;

  // And edges:
  this->edge[0][0] = 0;
  this->edge[0][1] = si;
  this->edge[1][0] = 1;
  this->edge[1][1] = si;
  this->edge[2][0] = 2;
  this->edge[2][1] = si;

  this->S = 1;

  // Now build the rest of the topology
  for(k = 0; k < l; k++) {
    // Next vertex index
    int en = k + 3;
    // Next Steiner point index
    m++;
    // Next Steiner point index in table
    si++;
    // Edge to split index
    int e = topo_vec[k];
    
    // End points of edge
    unsigned int ea, eb;
    ea = this->edge[e][0];
    eb = this->edge[e][1];
    // Insert new connections in adj
    this->adj[m][0] = ea;
    this->adj[m][1] = eb;
    this->adj[m][2] = en;

    // Check if ea is a sp and update adj if it is
    if(ea >= this->N) {
      for(int ed = 0; ed < 3; ed++) {
	if(this->adj[ea-this->N][ed] == (int)eb) {
	  this->adj[ea-this->N][ed] = si;
	  break;
	}
      }
    }
    // Similar check for eb
    if(eb >= this->N) {
      for(int ed = 0; ed < 3; ed++) {
	if(this->adj[eb-this->N][ed] == (int)ea) {
	  this->adj[eb-this->N][ed] = si;
	  break;
	}
      }
    }
    
    // Set edge e end point
    this->edge[e][1] = si;

    // Add two new edges at 2en-4 = 2(k+3)-4 = 2k+2
    e = 2*en - 3;
    
    // Add first edge (en -- si)
    this->edge[e][0] = en;
    this->edge[e][1] = si;
    
    // Add second edge (eb - si)
    e++;
    this->edge[e][0] = eb;
    this->edge[e][1] = si;

    // Lastly, add the new steiner point point
    Point p1,p2,p3;
    p1 = this->P[ea];
    p2 = this->P[eb];
    p3 = this->P[en];
    for(d = 0; d < this->dim; d++) {
      this->P[si][d] = (p1[d] + p2[d] + p3[d]) / 3.0 + 0.0001 * Utils::frand();
    }
    this->S++;
  }
  return;
}
