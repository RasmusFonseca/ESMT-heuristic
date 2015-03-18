#include <vector>
#include <list>
#include <algorithm>
#include <iostream>

#include "steiner/graph.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/disjoint_set.hpp"
#include "steiner/utils/point.hpp"

typedef Utils::Edge         Edge;
typedef Utils::Point        Point;

/*
 * Compare function for the std::sort in Utils::MSTKruskal
 */
bool compareLength(const Edge &e1, const Edge &e2) {
  return e1.length < e2.length;
}


/*
 * Implementation of Kruskal's MST algorithm
 *
 * Runs in O(n lg n)
 */
Graph Utils::MSTKruskal(Graph &graph) {
  // Iterators
  std::vector<Point>::iterator pit;
  std::vector<Edge>::iterator eit;

  // Terminals and edges
  std::vector<Point> *points = graph.getPointsPtr();
  std::vector<Edge>  *edges  = graph.getEdgesPtr();
  
  Graph result = Graph(*points);
  // Result
  std::vector<Edge> *res_edges = result.getEdgesPtr();

  // List of disjoint sets
  std::vector< Utils::DisjointSet<Point*> > sets;

  // Create the disjoint sets
  for(pit = points->begin(); pit != points->end(); ++pit) {
    Utils::DisjointSet<Point*> ds(&(*pit));
    sets.push_back(ds);
  }
  
  // Sort the list of edges according to length
  std::sort(edges->begin(), edges->end(), compareLength);

  // Iterate through the edges
  for(eit = edges->begin(); eit != edges->end(); ++eit) {
    // Get the two sets from the list of sets
    Utils::DisjointSet<Point*> *ds1 = &sets[eit->i0];
    Utils::DisjointSet<Point*> *ds2 = &sets[eit->i1];

    // Check if we may safely concatenate
    if(ds1->findSet() != ds2->findSet()) {
      // Add edge to the resulting list and union
      res_edges->push_back(*eit);
      ds1->setUnion(ds2);
    }
  }

  return result;
}

double Utils::frand() {
  return (double)rand() / RAND_MAX;
}

void bDist(std::unordered_map<long, double> &b_tree,
	   std::vector<Edge> &edges, int org, int cur, int prev, double longest) {
  unsigned int i;
  unsigned long key = Edge::key(org,cur);
  if(org != cur && b_tree.find(key) != b_tree.end()) {
    if(longest > b_tree[key])
      b_tree[key] = longest;
  }
  for(i = 0; i < edges.size(); i++) {
    double l = edges[i].length > longest ? edges[i].length : longest;
    if(edges[i].i0 == cur && edges[i].i1 != prev)
      bDist(b_tree, edges, org, edges[i].i1, edges[i].i0, l);
    else if(edges[i].i1 == cur && edges[i].i0 != prev)
      bDist(b_tree, edges, org, edges[i].i0, edges[i].i1, l);
  }
}

void Utils::constructBottleneckGraph(std::unordered_map<long, double> &b_tree,
				     Graph &del, Graph &mst) {
  unsigned int i;

  for(i = 0; i < del.getEdgesPtr()->size(); i++)
    b_tree[(*del.getEdgesPtr())[i].key()] = -1.0;
  
  std::vector<Edge> *edges = mst.getEdgesPtr();
  for(i = 0; i < mst.getPointsPtr()->size(); i++)
    bDist(b_tree, *edges, i, i, -1, -1.0);
}

void validateRec(std::vector<Edge> &edges, std::vector<int> &pf,
		 std::vector<int> &ef, int cur, int prev) {
  unsigned int i;
  pf[cur]++;
  for(i = 0; i < edges.size(); i++) {
    if(edges[i].i0 == cur && edges[i].i1 != prev) {
      ef[i]++;
      if(ef[i] > 1)
	return;
      validateRec(edges, pf, ef, edges[i].i1, edges[i].i0);
    }
    else if(edges[i].i1 == cur && edges[i].i0 != prev) {
      ef[i]++;
      if(ef[i] > 1)
	return;
      validateRec(edges, pf, ef, edges[i].i0, edges[i].i1);
    }
  }
}

bool Utils::validate(SteinerTree &st) {
  unsigned int i, j, p = 0, sp = 0;
  std::vector<Point> *points = st.getPointsPtr();
  std::vector<Edge>  *edges  = st.getEdgesPtr();
  // No of edges check
  if(points->size() != edges->size()+1) {
    std::cerr << "Wrong number of edges!!!" << std::endl;
    return false;
  }

  for(i = 0; i < points->size(); i++) {
    if((*points)[i].isSteiner())
      sp++;
    else
      p++;
  }
  if(sp > p-2) {
    std::cerr << "Too many Steiner points: " << sp << " vs. " << p
	      << " terminals." << std::endl;
    return false;
  }
  
  std::vector<int> pf(points->size(), 0);
  std::vector<int> ef(edges->size(), 0);
  validateRec(*edges, pf, ef, 0, -1);

  for(i = 0; i < pf.size(); i++) {
    if(pf[i] == 0) {
      std::cerr << "Point not reached: " << i << " " << (*points)[i] << std::endl;
      return false;
    }
    else if(pf[i] > 1) {
      std::cerr << "Point reached more than once: "
		<< i << " " << (*points)[i] << std::endl;
      return false;
    }
  }
  for(i = 0; i < ef.size(); i++) {
    if(!ef[i]) {
      std::cerr << "Edge not reached: " << i << " (" << (*edges)[i].i0
		<< " : " << (*edges)[i].i1 << ")" << std::endl;
      return false;
    }
    else if(ef[i] > 1) {
      std::cerr << "Edge reached more than once: " << i << " (" << (*edges)[i].i0
		<< " : " << (*edges)[i].i1 << ")" << std::endl;
      return false;
    }
  }

  Graph g;
  for(i = 0; i < p; i++) {
    g.getPointsPtr()->push_back((*points)[i]);
    for(j = i+1; j < p; j++) {
      g.getEdgesPtr()->push_back(Edge(i,j,Utils::length((*points)[i],(*points)[j])));
    }
  }
  g = Utils::MSTKruskal(g);
  double mst_length = g.getLength();
  if(st.getMSTLength()-mst_length>0.0001) {
    std::cerr << "MST-length is not correct!" << std::endl
	      << "  Correct |MST|: " << mst_length << std::endl
	      << "  Actual |MST|:  " << st.getMSTLength() << std::endl;
    return false;
  }

  return true;
}
