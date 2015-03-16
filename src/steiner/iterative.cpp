#include <vector>
#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <cmath>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/iterative.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/utils.hpp"

typedef Utils::Point Point;
typedef Utils::Edge  Edge;

/**
 * Constructor
 */
Iterative::Iterative(int dim, int n) {
  this->alloc   = n;
  this->max_dim = dim;
  // Allocate resources
  this->adj    = new int*[n];
  this->edge   = new int*[2*n];
  this->EL     = new double*[n];
  this->B      = new double*[n];
  this->C      = new double*[n];
  this->val    = new int[n];
  for(int i = 0; i < n; i++) {
    this->adj[i] = new int[3];
    this->EL[i]  = new double[3];
    this->B[i]   = new double[3];
    this->C[i]   = new double[dim];
  }
  for(int i = 0; i < 2*n; i++) {
    this->edge[i] = new int[2];
  }
}

/**
 * Destructor
 */
Iterative::~Iterative() {
  // Free resources
  for(int i = 0; i < this->alloc; i++) {
    delete this->adj[i];
    delete this->EL[i];
    delete this->B[i];
    delete this->C[i];
  }
  for(int i = 0; i < 2*this->alloc; i++) {
    delete this->edge[i];
  }
  delete this->adj;
  delete this->edge;
  delete this->EL;
  delete this->B;
  delete this->C;
  delete this->val;
}

void Iterative::optimise(double tol) {
  this->leafQ.clear();
  this->eqnstack.clear();
  
  int n0, n1, n2, i, i2, j;
  unsigned int m, d;
  double q0, q1, q2, t;
  
  for(i = this->S - 1; i >= 0; i--) {
    n0 = this->adj[i][0];
    n1 = this->adj[i][1];
    n2 = this->adj[i][2];
    q0 = 1.0 / (this->EL[i][0] + tol);
    q1 = 1.0 / (this->EL[i][1] + tol);
    q2 = 1.0 / (this->EL[i][2] + tol);
    t = q0 + q1 + q2;
    q0 /= t;
    q1 /= t;
    q2 /= t;
    this->val[i] = 0;
    this->B[i][0] = 0.0;
    this->B[i][1] = 0.0;
    this->B[i][2] = 0.0;
    for(d = 0; d < this->dim; d++)
      this->C[i][d] = 0.0;
    this->prep(i, 0, n0, q0);
    this->prep(i, 1, n1, q1);
    this->prep(i, 2, n2, q2);
    
    if(this->val[i] <= 1)
      this->leafQ.push_back(i);
  }
  
  // Solve
  while(this->leafQ.size() > 1) {
    i = this->leafQ.back();
    this->leafQ.pop_back();
    this->val[i]--;
    i2 = i + this->N;
    this->eqnstack.push_back(i);
    for(j = 0; j < 3; j++) {
      if(this->B[i][j] != 0.0)
	break;
    }
    
    // Neighbour is #j
    q0 = this->B[i][j];
    j = this->adj[i][j] - this->N;
    this->val[j]--;
    if(this->val[j] == 1)
      this->leafQ.push_back(j);

    for(m = 0; m < 3; m++) {
      if(this->adj[j][m] == i2)
	break;
    }
    
    q1 = this->B[j][m];
    this->B[j][m] = 0.0;
    t = 1.0 / (1.0 - q1 * q0);
    
    for(m = 0; m < 3; m++)
      this->B[j][m] *= t;
    for(m = 0; m < this->dim; m++) {
      this->C[j][m] += q1 * this->C[i][m];
      this->C[j][m] *= t;
    }
  }
  
  // One leaf left
  i = this->leafQ.back();
  this->leafQ.pop_back();
  i2 = i + this->N;
  for(m = 0; m < this->dim; m++) {
    this->P[i2][m] = this->C[i][m];
  }
  
  // Back solve
  while(this->eqnstack.size() > 0) {
    i = this->eqnstack.back();
    this->eqnstack.pop_back();
    i2 = i + this->N;
    for(j = 0; j < 3; j++) {
      if(this->B[i][j] != 0.0)
	break;
    }
    q0 = this->B[i][j];
    j = this->adj[i][j];
    for(m = 0; m < this->dim; m++)
      this->P[i2][m] = this->C[i][m] + q0 * this->P[j][m];
  }
  return;
}

// Calculates valens and B/C arrays
void Iterative::prep(int i, int a, int b, double c) {
  if(b >= (int)this->N) {
    // Steiner Point
    this->val[i]++;
    this->B[i][a] = c;
  }
  else {
    // Normal point
    for(unsigned int d = 0; d < this->dim; d++) {
      this->C[i][d] += this->P[b][d]*c;
    }
  }
}

// Computes total length of tree
// Furthermore fills the EL array with lengths of adj-edges
double Iterative::length() {
  double total = 0.0;
  double t;
  unsigned i, j;
  int i2, n0, n1, n2;
  for(i = 0; i < this->S; i++) {
    // Index of SP in points
    i2 = i + this->N;
    n0 = this->adj[i][0];
    n1 = this->adj[i][1];
    n2 = this->adj[i][2];
    // Get distance from i to n0
    if(n0 < i2) {
      t = this->dist(n0,i2);
      total += t;
      this->EL[i][0] = t;
      n0 -= this->N;
      if(n0 >= 0) {
	// Update other SP (n0) EL entry too
	for(j = 0; j < 3; j++) {
	  if(this->adj[n0][j] == i2) {
	    this->EL[n0][j] = t;
	    break;
	  }
	}
      }
    }
    // Get distance from i to n1
    if(n1 < i2) {
      t = this->dist(n1,i2);
      total += t;
      this->EL[i][1] = t;
      n1 -= this->N;
      if(n1 >= 0) {
	// Update other SP (n1) EL entry too
	for(j = 0; j < 3; j++) {
	  if(this->adj[n1][j] == i2) {
	    this->EL[n1][j] = t;
	    break;
	  }
	}
      }
    }
    // Get distance from i to n2
    if(n2 < i2) {
      t = this->dist(n2,i2);
      total += t;
      this->EL[i][2] = t;
      n2 -= this->N;
      if(n2 >= 0) {
	// Update other SP (n0) EL entry too
	for(j = 0; j < 3; j++) {
	  if(this->adj[n2][j] == i2) {
	    this->EL[n2][j] = t;
	    break;
	  }
	}
      }
    }
  }
  return total;
}

/* Implementation of error */
double Iterative::error() {
  unsigned int i, j;
  double res = 0.0;
  for(i = 0; i < this->S; i++) {
    int i2 = i + this->N;
    int n0 = this->adj[i][0];
    int n1 = this->adj[i][1];
    int n2 = this->adj[i][2];

    double d01=0, d12=0, d02=0, d0=0, d1=0, d2=0, pi2j=0;
    for(j = 0; j < this->dim; j++) {
      pi2j = P[i2][j];
      d0 = P[n0][j]-pi2j;
      d1 = P[n1][j]-pi2j;
      d2 = P[n2][j]-pi2j;
      d01 += d0*d1;
      d12 += d1*d2;
      d02 += d0*d2;
    }

    double t;
    t = 2*d01 + this->EL[i][0]*this->EL[i][1];
    if(t > 0.0)
      res += t;
    t = 2*d12 + this->EL[i][1]*this->EL[i][2];
    if(t > 0.0)
      res += t;
    t = 2*d02 + this->EL[i][0]*this->EL[i][2];
    if(t > 0.0)
      res += t;
  }
  return std::sqrt(res);
}

// Compute distance between points
double Iterative::dist(int i1, int i2) {
  return Utils::length(this->P[i1], this->P[i2]);
}

void Iterative::cleanUp(std::vector<Point> *points, std::vector<Edge> *edges,
			bool doCleanUp) {
  // Now we have to remove any (near-) zero edges.
  unsigned int i, j, k, n, sp, g, a, b;
  double l;

  points->clear();
  edges->clear();
  points->insert(points->end(), this->P.begin(), this->P.begin()+this->N);
  if(!doCleanUp) {
    points->insert(points->end(), this->P.begin()+this->N, this->P.end());
    for(i = this->N; i < points->size(); i++)
      (*points)[i].setSteiner();
    for(i = 0; i < this->S; i++) {
      k = i + this->N;
      for(j = 0; j < 3; j++) {
	n = this->adj[i][j];
	if(n < k)
	  edges->push_back(Edge(k, n, this->dist(k, n)));
      }
    }
  }
  else {
    // Remove a SP if dist between SP and  terminal is less than
    // 0.01 % of the avg. length of an edge in the tree
    l = this->length();
    l /= this->N+this->S;
    l *= 0.0001;
    
    std::vector< std::vector<int> > tadj(this->N);
    std::vector<int> map(this->S);
    // Adj is in order
    for(i = 0; i < this->S; i++) {
      for(j = 0; j < 3; j++) {
	n = this->adj[i][j];
	if(n < this->N) {
	  tadj[n].push_back(i+this->N);
	}
      }
    }
    
    g = this->N + this->S + 1;
    for(i = 0; i < this->N; i++) {
      for(j = 0; j < tadj[i].size(); j++) {
	sp = tadj[i][j];
	if(sp >= this->N && this->dist(sp, i) < l) {
	  // Remove the connected SP.
	  n = sp-this->N;
	  for(k = 0; k < 3; k++)
	    if(this->adj[n][k] == (int)i)
	      break;
	  a = this->adj[n][(k+1)%3];
	  b = this->adj[n][(k+2)%3];
	  // Update terminal edges
	  tadj[i][j] = a;
	  tadj[i].push_back(b);
	  // Update neighbours
	  this->updateRef(a, sp, i, tadj);
	  this->updateRef(b, sp, i, tadj);
	  // Make sure not to add edges from sp
	  for(k = 0; k < 3; k++)
	    this->adj[n][k] = g;
	  j--; // Check new edges.
	}
      }
    }
    // Add sps
    for(i = 0; i < this->S; i++) {
      if(this->adj[i][0] != (int)g) {
	map[i] = points->size();
	points->push_back(this->P[i+this->N]);
	points->back().setSteiner();
      }
    }
    // Add terminal-terminal edges.
    for(i = 0; i < this->N; i++) {
      for(j = 0; j < tadj[i].size(); j++) {
	n = tadj[i][j];
	if(n < i)
	  edges->push_back(Edge(i, n, this->dist(i, n)));
      }
    }
    for(i = 0; i < this->S; i++) {
      k = i + this->N;
      for(j = 0; j < 3; j++) {
	n = this->adj[i][j];
	if(n < k) {
	  if(n >= this->N)
	    n = map[n-this->N];
	  edges->push_back(Edge(map[k-this->N], n, this->dist(map[k-this->N], n)));
	}
      }
    }
  }
}

void Iterative::updateRef(unsigned int i, int oldp, int newp,
			  std::vector< std::vector<int> > &tadj) {
  unsigned int j, k;
  if(i >= this->N) {
    k = i-this->N;
    for(j = 0; j < 3; j++)
      if(this->adj[k][j] == oldp)
	break;
    this->adj[k][j] = newp;
  }
  else {
    for(j = 0; j < tadj[i].size(); j++)
      if(tadj[i][j] == oldp)
	break;
    tadj[i][j] = newp;
  }
}
