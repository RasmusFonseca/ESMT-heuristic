#include <vector>
#include <algorithm>
#include <unordered_map>
#include <assert.h>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/esmt.hpp"
#include "steiner/iterative.hpp"
#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/heuristics/concat.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/heuristics/steiner_finder.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/fermat.hpp"
#include "steiner/utils/disjoint_set.hpp"
#include "steiner/utils/delaunay.hpp"

#define ESMT_COLLECT_STATS  1
#define DO_USE_B_GRAPH      0
#define DO_USE_B_GRAPH_FULL 0

/*
 * Utils type definitions
 */
typedef Utils::Point                 Point;
typedef Utils::Edge                  Edge;
typedef Utils::DisjointSet<int>      DisjointSet;
typedef Utils::Delaunay              Delaunay;
typedef Utils::Delaunay::Simplex     Simplex;
typedef Utils::Delaunay::PointHandle PointHandle;

/*
 * Constructor
 */
ESMT::ESMT(std::vector<Point> &points)
: SteinerTree(), Iterative(points[0].dim(), points.size()),
	iterCon(points[0].dim()), iterSmith(points[0].dim())
{
	this->points = points;
	SubgraphHeuristic *sh = new IterativeSmith(points[0].dim());
	this->findESMT(points, sh, true, true, false);
	delete sh;
}

/*
 * Constructor
 */
ESMT::ESMT(std::vector<Point> &points, SubgraphHeuristic *sh,
		bool concat_subgraphs, bool post_optimise, bool special_concat,
		bool verbose)
: SteinerTree(), Iterative(points[0].dim(), points.size()),
	iterCon(points[0].dim()), iterSmith(points[0].dim())
{
	if(verbose)
		std::cout << "*** ESMT Heuristic ***" << std::endl;
	this->findESMT(points, sh, concat_subgraphs, post_optimise, special_concat, verbose);
	if(verbose)
		std::cout << "**********************" << std::endl;
}

/*
 * Constructor
 */
ESMT::ESMT(Utils::Delaunay &del, SubgraphHeuristic *sh,
		bool concat_subgraphs, bool post_optimise, bool special_concat,
		bool verbose)
: SteinerTree(), Iterative(del.dimension(), del.getPointsPtr()->size()),
	iterCon(del.dimension()), iterSmith(del.dimension())
{
	if(verbose)
		std::cout << "*** ESMT Heuristic ***" << std::endl;

	this->findESMT(del, sh, concat_subgraphs, post_optimise, special_concat, verbose);
	/* Disabled old special concat. New version: added to org. heuristic,
	 * d-simplices only
	 if(special_concat)
	 this->findESMTSct(del, sh, concat_subgraphs, post_optimise, verbose);
	 else
	 this->findESMT(del, sh, concat_subgraphs, post_optimise, verbose);
	 */
	if(verbose)
		std::cout << "**********************" << std::endl;
}

ESMT::Stats *ESMT::getStats() {
	return &this->stats;
}

std::vector<ESMT::SubST> &ESMT::getComponents() {
	return this->smts;
}

/*
 * Implementation of findESMT(...)
 */
void ESMT::findESMT(std::vector<Point> &points,
		SubgraphHeuristic *sh,
		bool concat_subgraphs,
		bool post_optimise,
		bool special_concat,
		bool verbose) {
	// Do Delaunay
	Delaunay del(points);
	/* Disabled for now : special_concat implemented in normal findESMT
	   if(special_concat)
	   this->findESMTSct(del, sh, concat_subgraphs, post_optimise, verbose);
	   else*/
	this->findESMT(del, sh, concat_subgraphs, post_optimise, special_concat, verbose);
}

/*
 * Implementation of findESMT(...) with Delaunay as argument
 */
void ESMT::findESMT(Delaunay &del, 
		SubgraphHeuristic *sh,
		bool concat_subgraphs,
		bool post_optimise,
		bool special_concat,
		bool verbose) {
	unsigned int i,j,n,c;

	assert(sh);
	this->dim = del.dimension();
	this->points = *del.getPointsPtr();
	this->N = this->points.size();

	this->edges              = *del.getEdgesPtr();
	this->simplices          = del.getSimplices();
	this->is_covered_simplex = std::vector<bool>(this->simplices->size(), false);
	std::vector<PointHandle> *handles   = del.getPointHandles();

#if(ESMT_COLLECT_STATS)
	this->stats.no_of_simplices = this->simplices->size();
#endif

	if(verbose) {
		std::cout << "Delaunay done!" << std::endl
			<< "  Simplices:  " << this->simplices->size() << std::endl;
	}

	// Get the mst. This is simply MST Kruskal, only we mark each edge in in_MST
	// List of disjoint sets
	std::vector< DisjointSet > sets;
	// Create sets with indicies
	for(i = 0; i < this->N; i++)
		sets.push_back(DisjointSet(i));

	// Sort the list of edges according to length
	std::sort(this->edges.begin(), this->edges.end(), compareLength);

	std::vector<Edge> mst_edges;
	double mst_length = 0.0;
	c                 = 0;

	// Iterate through the edges
	std::vector<Edge>::iterator  eit;
	for(eit = this->edges.begin(); eit != this->edges.end(); ++eit) {
		// Get the two sets from the list of sets
		DisjointSet *ds1 = &sets[eit->i0];
		DisjointSet *ds2 = &sets[eit->i1];

		// Check if we may safely concatenate
		if(ds1->findSet() != ds2->findSet()) {
			// Add edge to the resulting list and union
			mst_edges.push_back(*eit);
			// Mark the edge as in MST
			in_MST[eit->key()] = true;
			// Length
			mst_length += eit->length;
			// Union the sets
			ds1->setUnion(ds2);
			c++;
			if(c == this->N - 1)
				// All points connected
				break;
		}
	}

	// Set edges to mst_edges
	this->edges = mst_edges;
	this->setMSTLength(mst_length);

	if(verbose) {
		std::cout << "MST done!" << std::endl
			<< "  Length:  " << mst_length << std::endl;
	}

	// Get components to concatenate
	std::vector< Subset > components;
	components.reserve(this->points.size());

	//////////////////
	// Preprocessing

	// 1. Build the connections list and sort all handles
	n = this->points.size();
	std::vector< std::vector<int> > connections;
	for(i = 0; i < n; i++) {
		connections.push_back(std::vector<int>());
		std::sort((*handles)[i].simplices.begin(),(*handles)[i].simplices.end());
	}
	int i0, i1;
	for(i = 0; i < this->edges.size(); i++) {
		i0 = this->edges[i].i0;
		i1 = this->edges[i].i1;
		connections[i0].push_back(i1);
		connections[i1].push_back(i0);
	}

	// 2. Generate faces
	for(i = 0; i < n; i++)
		this->findFaces(*handles, components, connections, i);

	if(verbose) {
		std::cout << "Preprocessing done!" << std::endl
			<< "  Number of covered faces: " << components.size() << std::endl;
#if(ESMT_COLLECT_STATS)
		for(i = 2; i < this->stats.covered_faces.size(); i++)
			std::cout << "   [" << i << "]: "
				<< this->stats.covered_faces[i] << std::endl;
#endif
	}

	// We now have a list of all faces covered by the MST.
	// Begin finding small SMTs.
	// Allocate vector now. Approx each element should be added.
	this->smts.reserve(components.size());

	std::vector< Subset >::iterator sit;
	Graph submst;
	for(sit = components.begin(); sit != components.end(); sit++) {
		if(sit->map.size() == 2) {
			SteinerTree *st = new SteinerTree();
			this->buildMST(*sit, *st);
			st->setMSTLength(st->getLength());
			st->getSteinerRatio();

			this->smts.push_back(SubST(st,sit->map.size(),&sit->map[0]));
		}
		else {
			this->buildMST(*sit, submst);
			SteinerTree *st;
			if(sit->map.size() == 3)
				st = Utils::getFermatSMT(submst);
			else
				st = sh->findSteinerPoints(submst);

			if(st->getSteinerRatio() >= 1.0) {
				delete st;
				continue;
			}
			SubST subst(st,sit->map.size(),&sit->map[0]);
			this->smts.push_back(subst);
			if(sit->map.size() == this->dim+1 && concat_subgraphs) {
				int simplex = sit->simplex_index;
				std::vector<int> map;
				for(i = 0; i < (*this->simplices)[simplex].n; i++)
					map.push_back((*this->simplices)[simplex].map[i]);
				for(i = 0; i < map.size(); i++) {
					// Swap first element with i
					int tmp = map[0];
					map[0]  = map[i];
					map[i]  = tmp;
					this->buildSausage(map, subst, 0, i, simplex, simplex, 1);
				}
			}
		}
	}

#if(ESMT_COLLECT_STATS)
	this->stats.sub_trees_in_queue = this->smts.size();
#endif

	if(verbose) {
		std::cout << "Sub-trees found!" << std::endl
			<< "  Number of sub-trees in queue: " << this->smts.size() << std::endl;
#if(ESMT_COLLECT_STATS)
		std::cout << "  Number of covered sausages: " << std::endl;
		for(i = this->dim+2; i < this->stats.covered_faces.size(); i++)
			std::cout << "   [" << i << "]: "
				<< this->stats.covered_faces[i] << std::endl;
#endif
	}

	// Sort priority queue of sub-graphs according to ratio
	std::sort(this->smts.begin(), this->smts.end(), ESMT::compareSteinerRatio);

	if(special_concat) {
		// Starting concat
		this->doConcatenate(false);
		std::vector<SubST> non_covered;
		// Preprocessing: find SMTs of all non_covered
		std::vector< Simplex >::iterator spit;
		for(i = 0, spit = this->simplices->begin();
				spit != this->simplices->end(); i++, spit++) {
			if(this->is_covered_simplex[i])
				// Simplex is covered. Continue
				continue;
			// All is of size = d+1 (since only simplices)
			Graph submst;
			for(i = 0; i < spit->n; i++) {
				submst.getPointsPtr()->push_back(this->points[spit->map[i]]);
				for(j = i+1; j < spit->n; j++) {
					int a = spit->map[i];
					int b = spit->map[j];
					submst.getEdgesPtr()->push_back(Edge(i,j,Utils::length(this->points[a],this->points[b])));
				}
			}

			submst = Utils::MSTKruskal(submst);
			submst.setMSTLength(submst.getLength());

			SteinerTree *st;
			if(spit->n == 3)
				st = Utils::getFermatSMT(submst);
			else
				st = sh->findSteinerPoints(submst);
			if(st->getSteinerRatio() >= 1.0) {
				delete st;
				continue;
			}
			SubST subst(st,spit->n,spit->map);
			non_covered.push_back(subst);
		}

		std::sort(non_covered.begin(), non_covered.end(), ESMT::compareSteinerRatio);

		// Now, redo concatenation
		std::vector<SubST> pre_add;
		double best_ratio = this->getSteinerRatio();
		unsigned int k = 1000; // Maximum of max(1000,non_covered.size()) tries
		if(verbose)
			std::cout << "Ratio after initial concatenation: "
				<< best_ratio << std::endl;
		for(i = 0; i < k && i < non_covered.size(); i++) {
			pre_add.push_back(non_covered[i]);
			// Now try to concatenate again.
			this->points.erase(this->points.begin()+this->N, this->points.end());
			this->doConcatenate(pre_add, false);
			this->setSMTLength(-1);
			this->setSteinerRatio(-1);
			if(this->getSteinerRatio() < best_ratio) {
				if(verbose)
					std::cout << "Better ratio achieved: " << this->getSteinerRatio()
						<< std::endl;
				best_ratio = this->getSteinerRatio();
			}
			else {
				// Remove this FST from pre_add again.
				pre_add.pop_back();
			}
		}

#if(ESMT_COLLECT_STATS)
		this->stats.added_sub_trees.clear();
		this->stats.add_sub_trees_total = 0;
#endif

		this->points.erase(this->points.begin()+this->N, this->points.end());
		this->doConcatenate(pre_add, verbose);
		this->setSMTLength(-1);
		this->setSteinerRatio(-1);
	}
	else {
		// Simple one-time concat
		this->doConcatenate(verbose);
	}

#if(ESMT_COLLECT_STATS)
	this->stats.no_of_sp = this->points.size()-this->N;
#endif

	if(post_optimise) {
		if(verbose)
			std::cout << "Starting post optimisation." << std::endl
				<< "  Current |SP| = " << this->points.size()-this->N << std::endl;
		// Setup for Smiths and add extra SPs. Then optimise, do clean-up and return.
		this->postOptimisation();

#if(ESMT_COLLECT_STATS)
		// Find overlapping SPs
		double l = this->length();
		l /= this->N+this->S;
		for(i = 0; i < this->N-2; i++) {
			for(j = 0; j < 3; j++) {
				unsigned int op = this->adj[i][j];
				if(op >= this->N && op-this->N > i) {
					if(this->EL[i][j] < 0.0001*l)
						this->stats.no_of_sp_overlapping++;
				}
			}
		}
#endif

		// Convert to ESMT structure
		this->points.erase(this->points.begin()+this->N, this->points.end());
		this->edges.clear();

		this->cleanUp(&this->points, &this->edges);

#if(ESMT_COLLECT_STATS)
		this->stats.no_of_sp_post_optimisation = this->points.size()-this->N;
#endif

		if(verbose) {
			std::cout << "Post optimisation done!" << std::endl
				<< "  |SP| = " << this->points.size()-this->N << std::endl;
#if(ESMT_COLLECT_STATS)
			std::cout << "  Number of short SP-SP edges: "
				<< this->stats.no_of_sp_overlapping << std::endl;
#endif
		}
	}
}

/*
 * Destructor
 */
ESMT::~ESMT() {
	// Do clean-up of components
	std::vector<SubST>::iterator stit;
	for(stit = this->smts.begin(); stit != this->smts.end(); stit++) {
		delete stit->st;
	}
}

// Private functions

/*
 * Performs the concatenation step of the algorithm.
 */
void ESMT::doConcatenate(bool verbose) {
	std::vector<SubST> dummy;
	this->doConcatenate(dummy, verbose);
}

/*
 * Performs the concatenation step of the algorithm.
 */
void ESMT::doConcatenate(std::vector<SubST> &pre_add_list, bool verbose) {
	unsigned int i, c, index;
	// Concatenation process.
	// This is really just MSTKruskal.

	// Clear edges from Delaunay
	this->edges.clear();

	if(verbose) {
		std::cout << "Starting concatenation process." << std::endl;
	}  

	c = 0;
	// First, create Disjoint sets. Each set is build over vertex indicies
	std::vector< DisjointSet > sets;
	for(i = 0; i < this->points.size(); i++) {
		sets.push_back(DisjointSet(i));
	}

	// Flags used to determine if two sets are disjoint
	bool *flags = new bool[this->N];
	for(i = 0; i < this->N; i++)
		flags[i] = false;

	std::vector<SubST>::iterator sit;
	for(sit = pre_add_list.begin(); sit != this->smts.end(); ++sit) {
		// If we are done with pre_add_list -> start on smts
		if(sit == pre_add_list.end())
			sit = this->smts.begin();
		// Check all terminals from sub-tree. If any of them has same set ->
		// Bad join

		// Indexes is a map of the graph vertices to the full graph vertices
		// indexes[j] is the index of vertice j in the big graph
		//std::vector<int> indexes = sit->map;
		// A list of all flags set. Used for fast clean-up
		std::vector<int> flags_set;
		for(i = 0; i < sit->n; i++) {
			// id is the vertex index in the full graph
			int id = sets[sit->map[i]].findSet()->getValue();
			if(flags[id])
				// Not good, has been encountered before = same set
				break;
			// Set flag
			flags[id] = true;
			flags_set.push_back(id);
		}
		if(i == sit->n) {	
			// No conflicts, this graph may be inserted without problem
			SteinerTree *st = sit->st;

			std::vector<Point> *points = st->getPointsPtr();
			std::vector<Edge>  *edges  = st->getEdgesPtr();
			std::vector<Point>::iterator pit;
			std::vector<Edge>::iterator  eit;
			unsigned int no_of_terminals = points->size();

			// Add points
			index = 0;
			std::vector<int> indexes;
			for(pit = points->begin(); pit != points->end(); pit++) {
				if(pit->isSteiner()) {
					no_of_terminals--;
					this->points.push_back(*pit);
					indexes.push_back(this->points.size()-1);
				}
				else
					indexes.push_back(sit->map[index++]);
			}
			// Add edges
			for(eit = edges->begin(); eit != edges->end(); eit++) {
				Edge e(indexes[eit->i0],indexes[eit->i1]);
				this->edges.push_back(e);
			}
			// Now union all
			// All points in flag_set should be unioned
			for(i = 1; i < flags_set.size(); i++) {
				sets[flags_set[0]].setUnion(&sets[flags_set[i]]);
			}

#if(ESMT_COLLECT_STATS)
			while(this->stats.added_sub_trees.size() < no_of_terminals+1)
				this->stats.added_sub_trees.push_back(0);
			this->stats.added_sub_trees[no_of_terminals]++;
			this->stats.add_sub_trees_total++;
#endif

			// Increment counter
			c += no_of_terminals-1;
			// Check if we have added N-1 mst-edges
			if(c == this->N-1)
				// Connected graph
				break;
		}

		// Reset flags
		for(i = 0; i < flags_set.size(); i++) {
			flags[flags_set[i]] = false;
		}
	}

	if(verbose) {
		std::cout << "Concatenation done." << std::endl;
#if(ESMT_COLLECT_STATS)
		std::cout << "  Sub-trees added: " << std::endl;
		for(i = 2; i < this->stats.added_sub_trees.size(); i++)
			std::cout << "   [" << i << "]: "
				<< this->stats.added_sub_trees[i] << std::endl;
		std::cout << "   Total: " << this->stats.add_sub_trees_total << std::endl;
#endif
	}

	delete flags;
}

/* Implementation of postOptimisation() */
void ESMT::postOptimisation() {
	unsigned int i, j, k, v;

	this->S = this->points.size() - this->N;
	this->P.clear();
	for(i = 0; i < this->N+this->S; i++)
		this->P.push_back(this->points[i]);
	for(i = this->N+this->S; i < 2*this->N-2; i++)
		this->P.push_back(Point(this->dim));

	for(i = 0; i < this->N; i++)
		for(j = 0; j < 3; j++)
			this->adj[i][j] = -1;

	std::vector< std::vector<int> > tadj(this->N);

	// Fill adj and tadj array
	unsigned int i0, i1, is0, is1;
	for(i = 0; i < this->edges.size(); i++) {
		i0 = this->edges[i].i0;
		i1 = this->edges[i].i1;
		is0 = i0-this->N;
		is1 = i1-this->N;
		if(i0 >= this->N)
			for(j = 0; j < 3; j++) {
				if(this->adj[is0][j] < 0) {
					this->adj[is0][j] = i1;
					break;
				}
			}
		else
			tadj[i0].push_back(i1);
		if(i1 >= this->N)
			for(j = 0; j < 3; j++) {
				if(this->adj[is1][j] < 0) {
					this->adj[is1][j] = i0;
					break;
				}
			}
		else
			tadj[i1].push_back(i0);
	}

	// Adj and Tadj are now up to date.
	// Now do the inserts of new SPs.
	for(i = 0; i < this->N; i++) {
		// Insert SPs if valens > 1
		unsigned int val;
		while((val = tadj[i].size()) > 1) {
			// This terminal has valens > 1
			float smallest_angle = Utils::angle(this->P[tadj[i][0]], this->P[i], this->P[tadj[i][1]]);
			float tmp_angle;
			int   a, b;
			a = 0;
			b = 1;
			for(j = 0; j < val; j++) {
				for(k = j+1; k < val; k++) {
					tmp_angle = Utils::angle(this->P[tadj[i][j]], this->P[i], this->P[tadj[i][k]]);
					if(tmp_angle > smallest_angle) {
						smallest_angle = tmp_angle;
						a = j;
						b = k;
					}
				}
			}

			// Best angle is now aib. Insert point.
			this->P[this->N+this->S] = (this->P[tadj[i][a]] + this->P[i] + this->P[tadj[i][b]]) * (1.0 / 3.0);

			// Update edges
			this->adj[this->S][0] = tadj[i][a];
			this->adj[this->S][1] = tadj[i][b];
			this->adj[this->S][2] = i;

			if(tadj[i][a] >= (int)this->N) {
				int ai = tadj[i][a]-this->N;
				for(j = 0; j < 3; j++)
					if(this->adj[ai][j] == (int)i) {
						this->adj[ai][j] = this->N+this->S;
						break;
					}
			}
			else {
				v = tadj[tadj[i][a]].size();
				for(j = 0; j < v; j++)
					if(tadj[tadj[i][a]][j] == (int)i) {
						tadj[tadj[i][a]][j] = this->N+this->S;
						break;
					}
			}
			if(tadj[i][b] >= (int)this->N) {
				int bi = tadj[i][b]-this->N;
				for(j = 0; j < 3; j++)
					if(this->adj[bi][j] == (int)i) {
						this->adj[bi][j] = this->N+this->S;
						break;
					}
			}
			else {
				v = tadj[tadj[i][b]].size();
				for(j = 0; j < v; j++)
					if(tadj[tadj[i][b]][j] == (int)i) {
						tadj[tadj[i][b]][j] = this->N+this->S;
						break;
					}
			}

			// Update a
			tadj[i][a] = this->N+this->S;
			// Erase b
			tadj[i].erase(tadj[i].begin()+b);
			this->S++;
		}
	}

	// Now all terminals should have val = 1.
	// Optimise
	double l = this->length();
	double r = this->error(); 
	do {
		this->optimise(0.0001*r/this->N);
		l = this->length();
		r = this->error();
	} while(r>l*0.0001);
}

bool ESMT::compareSteinerRatio(const SubST &st1, const SubST &st2) {
	double r1 = st1.st->getSteinerRatio();
	double r2 = st2.st->getSteinerRatio();
	if(st1.st->getPointsPtr()->size() == 2 && st2.st->getPointsPtr()->size() == 2)
		return st1.st->getMSTLength() < st2.st->getMSTLength();
	return r1 < r2;
}

bool ESMT::compareLength(Edge e1, Edge e2) {
	return e1.length < e2.length;
}

bool ESMT::isInMST(int i0, int i1) {
	return (this->in_MST.find(Edge::key(i0,i1)) != this->in_MST.end());
}

void ESMT::buildMST(Subset &subset, Graph &g) {
	std::vector<Point> *points = g.getPointsPtr();
	std::vector<Edge>  *edges  = g.getEdgesPtr();
	points->clear();
	edges->clear();

	unsigned int i, j, i0, i1;
	for(i = 0; i < subset.map.size(); i++) {
		points->push_back(this->points[subset.map[i]]);
		for(j = i+1; j < subset.map.size(); j++) {
			i0 = subset.map[i];
			i1 = subset.map[j];
			if(this->isInMST(i0, i1))
				edges->push_back(Edge(i, j,
							Utils::length(this->points[i0], this->points[i1])));
		}
	}
	g.setMSTLength(g.getLength());
}

void ESMT::findFaces(std::vector< PointHandle > &handles,
		std::vector< Subset > &components,
		std::vector< std::vector<int> > &connections,
		unsigned int point) {

	unsigned int n, i, r, next;
	n = this->points.size();
	int map[n];
	int mapmax;
	mapmax = 0;
	for(i = 0; i < n; i++)
		map[i] = -1;
	n = connections[point].size();
	// Add 1 point set.
	r = components.size();
	components.push_back(Subset());
	components[r].map.push_back(point);

	map[point] = mapmax++;
	for(i = 0; i < n; i++) {
		next = connections[point][i];
		if(next > point)
			this->findFacesRec(handles, components, connections,
					handles[point].simplices,
					map, next, &mapmax, point, r);
	}
	// Remove one point set.
	components.erase(components.begin()+r);
}

void ESMT::findFacesRec(std::vector< PointHandle > &handles,
		std::vector< Subset > &components,
		std::vector< std::vector<int> > &connections,
		std::vector< int > &currentSimplices, int* map,
		int cur, int *mapmax, int bound, int prevSet) {

	unsigned int c, n, i, j;
	int next, rank;

	// Check, if we may add this point.
	// Assume point handle lists and currentSimplices to be sorted.
	std::vector<int> intersection;
	std::set_intersection(currentSimplices.begin(),currentSimplices.end(),
			handles[cur].simplices.begin(), 
			handles[cur].simplices.end(),
			std::back_inserter(intersection));
	if(intersection.size() == 0)
		// We cannot add this point :(
		return;
	// Add this point, and add the new set
	components.push_back(components[prevSet]);
	components[components.size()-1].map.push_back(cur);
	prevSet = components.size()-1;
	if(intersection.size() == 1 && components[prevSet].map.size() == this->dim+1) {
		// This is a simplex. Add index to simplex list
		components[prevSet].simplex_index = intersection[0];
		this->is_covered_simplex[intersection[0]] = true;
	}

#if(ESMT_COLLECT_STATS)
	unsigned int size = components[prevSet].map.size();
	while(this->stats.covered_faces.size() < size+1)
		this->stats.covered_faces.push_back(0);
	this->stats.covered_faces[size]++;
#endif

	// For each connection, try to add new point.
	n    = components[prevSet].map.size();
	rank = map[components[prevSet].map[n-1]];
	for(j = 0; j < n; j++) {
		cur = components[prevSet].map[j];
		c = connections[cur].size();
		for(i = 0; i < c; i++) {
			next = connections[cur][i];
			if(next > bound && (map[next] == -1 || map[next] > rank)) {
				if(std::find(components[prevSet].map.begin(),
							components[prevSet].map.end(), next)
						!= components[prevSet].map.end())
					continue;
				if(map[next] == -1)
					map[next] = (*mapmax)++;
				this->findFacesRec(handles, components, connections,
						intersection, map, next, mapmax, bound, prevSet);
			}
		}
	}
}

void ESMT::buildSausage(std::vector<int> &prevSet,
		SubST &prevTree,
		unsigned int added,
		unsigned int next_simplex_index,
		int cur_simplex,
		int org_simplex,
		unsigned int c) {
	unsigned int i, j, z;
	double d_mst_length = 0.0;

	Simplex *cur  = &(*this->simplices)[cur_simplex];
	int next_vertex = cur->nextVertex[next_simplex_index];
	int next_simplex = cur->nextSimplex[next_simplex_index];

	// Check for neighbour?
	if(next_vertex < 0)
		return;

	Simplex *next = &(*this->simplices)[next_simplex];

	// Detect if we are connected - must have at least c connections
	for(i = 0; i < cur->n; i++)
		if(this->isInMST(cur->map[i], next_vertex)) {
			c--;
			d_mst_length += Utils::length(this->points[cur->map[i]], this->points[next_vertex]);
		}
	if(c > 0)
		// Stop - but we may have to look in other direction here!
		return;
	// Check if we are adding a duplet.
	for(i = 0; i < prevSet.size(); i++)
		if(prevSet[i] == next_vertex)
			return;
	// Check if the sharing facet has d edges in MST
	c = 0;
	for(i = 0; i < cur->n; i++)
		if(i != next_simplex_index && this->isInMST(cur->map[i],cur->map[next_simplex_index]))
			c++;
	if(c == 1 && next_simplex < cur_simplex)
		return;
	std::vector<int> curSet = prevSet;
	curSet.push_back(next_vertex);
	// Is connected. Calculate ST and add.
	SteinerTree *st = this->iterCon.insertTerminal(prevTree.st, this->points[next_vertex], prevTree.st->getMSTLength()+d_mst_length);
	if(st->getSteinerRatio() >= 1.0) {
		delete st;
		return;
	}

#if(ESMT_COLLECT_STATS)
	unsigned int size = prevTree.n+1;
	while(this->stats.covered_faces.size() < size+1)
		this->stats.covered_faces.push_back(0);
	this->stats.covered_faces[size]++;
#endif

	SubST subst(st, prevTree.n, prevTree.map, next_vertex);
	this->smts.push_back(subst);  
	added += 1;
	// We build from the largest index, so we have to choose the
	// simplex on the face with the last dim points from cur-set.
	// However, at first, we have different directions, with one
	// less for each added simplex. Thus, we check (opposite) points from
	// size-(this->dim+1) to size-(this->dim+1-max(added,dim)), where
	c = added >= this->dim ? 1 : this->dim - added;
	z = curSet.size()-(this->dim+1);
	for(i = z; i < z+c; i++) {
		// Swap first element with i
		int tmp   = curSet[z];
		curSet[z] = curSet[i];
		curSet[i] = tmp;
		// Get index (j) for next simplex
		for(j = 0; j < next->n; j++)
			if(next->map[j] == curSet[z])
				break;
		// Now, next vertex is in next->nextVertex[j]. Recursive call now
		// if next simplex is less than bound.
		this->buildSausage(curSet, subst, added, j, next_simplex, org_simplex, 1);
	}
	// We have to start the concatenation in the other direction.
	// prevSet[0..(added-1)] contains the set indices, which must be respected.
	// If(added==dim+1), direction is uniqly determined, otherwise, we must make
	// A permutation of prevSet[added..dim]. First exit will be prevSet[dim].
	int start = curSet[0];
	// Do not build if neighbour point is <= start
	for(i = this->dim; i >= added; i--) {
		if(curSet[i] <= start)
			continue;
		int tmp           = curSet[this->dim];
		curSet[this->dim] = curSet[i];
		curSet[i]         = tmp;
		// Get index (j) for next simplex
		for(j = 0; j < (*this->simplices)[org_simplex].n; j++)
			if((*this->simplices)[org_simplex].map[j] == curSet[this->dim])
				break;
		// j is now index of the next connection
		this->buildSausageReverse(curSet, subst, added, 0, j, org_simplex);
	}
}

void ESMT::buildSausageReverse(std::vector<int> &prevSet,
		SubST &prevTree,
		unsigned int prev_added,
		unsigned int added,
		unsigned int next_simplex_index,
		int cur_simplex) {
	unsigned int i, j, z, c;
	double d_mst_length = 0.0;
	Simplex *cur  = &(*this->simplices)[cur_simplex];
	int next_vertex = cur->nextVertex[next_simplex_index];
	int next_simplex = cur->nextSimplex[next_simplex_index];
	// Check for neighbour?
	if(next_vertex < 0)
		return;

	Simplex *next = &(*this->simplices)[next_simplex];

	// Detect if we are connected - must have at least 1 connection
	for(i = 0; i < cur->n; i++)
		if(this->isInMST(cur->map[i], next_vertex)) {
			d_mst_length = Utils::length(this->points[cur->map[i]], this->points[next_vertex]);
			break;
		}
	if(i >= cur->n)
		// Not connected
		return;
	// Check if we are adding a duplet.
	for(i = 0; i < prevSet.size(); i++)
		if(prevSet[i] == next_vertex)
			return;
	// Check if the sharing facet has d edges in MST
	c = 0;
	for(i = 0; i < cur->n; i++)
		if(i != next_simplex_index && this->isInMST(cur->map[i],cur->map[next_simplex_index]))
			c++;
	if(c == 1 && next_simplex < cur_simplex)
		return;
	// Is connected. Calculate ST and add.
	std::vector<int> curSet;
	// Push front because of reverse build
	curSet.reserve(prevSet.size()+1);
	curSet.push_back(next_vertex);
	curSet.insert(curSet.end(), prevSet.begin(), prevSet.end());
	SteinerTree *st = this->iterCon.insertTerminal(prevTree.st, this->points[next_vertex], prevTree.st->getMSTLength()+d_mst_length);

	if(st->getSteinerRatio() >= 1.0) {
		delete st;
		return;
	}

#if(ESMT_COLLECT_STATS)
	unsigned int size = prevTree.n+1;
	while(this->stats.covered_faces.size() < size+1)
		this->stats.covered_faces.push_back(0);
	this->stats.covered_faces[size]++;
#endif

	SubST subst(st, prevTree.n, prevTree.map, next_vertex);
	this->smts.push_back(subst);
	added += 1;
	// We now build in reverse. Given curSet of the form:
	// curSet[0..(1st index)..((d+1)th index)..n]
	// = curSet[0..added..(added+this->dim+1)..n]
	//   c = added+prev_added is in the middle.
	// Consider all between c and (added+this->dim+1) for next direction.
	c = added+prev_added > this->dim ? this->dim : added+prev_added;
	z = this->dim;
	for(i = this->dim; i >= c; i--) {
		// Swap first element with i
		int tmp   = curSet[z];
		curSet[z] = curSet[i];
		curSet[i] = tmp;
		// Get index (j) for next simplex
		for(j = 0; j < next->n; j++)
			if(next->map[j] == curSet[z])
				break;
		// Now, next vertex is in next->nextVertex[j]. Recursive call now
		// if next simplex is less than bound.
		this->buildSausageReverse(curSet, subst, prev_added, added, j, next_simplex);
	}
}

/////////////////////////////////////////
// Struct functions

/* Subset constructor */
ESMT::Subset::Subset() {
	this->simplex_index = -1;
}

/* Subset assignment operator */
ESMT::SubST &ESMT::SubST::operator=(const ESMT::SubST &other) {
	if(this != &other) {
		delete this->map;

		this->n = other.n;
		this->map = new int[n];
		for(unsigned int i = 0; i < n; i++)
			this->map[i] = other.map[i];
		this->st = other.st;
		this->bmst_length = other.bmst_length;
	}
	return *this;
}

/* SubST constructor */
ESMT::SubST::SubST(SteinerTree *st, int n, int *map) {
	this->n   = n;
	this->st  = st;
	this->map = new int[n];
	for(int i = 0; i < n; i++)
		this->map[i] = map[i];
}

/* SubST constructor */
ESMT::SubST::SubST(SteinerTree *st, int n, int *map, int nidx) {
	this->n   = n+1;
	this->st  = st;
	this->map = new int[n+1];
	for(int i = 0; i < n; i++)
		this->map[i] = map[i];
	this->map[n] = nidx;
}

/* SubST copy constructor */
ESMT::SubST::SubST(const ESMT::SubST &s) {
	this->n   = s.n;
	this->st  = s.st;
	this->bmst_length = s.bmst_length;
	this->map = new int[s.n];
	for(unsigned int i = 0; i < n; i++)
		this->map[i] = s.map[i];
}

/* SubST destructor */
ESMT::SubST::~SubST() {
	delete this->map;
}

/* Stats constructor */
ESMT::Stats::Stats() {
	this->no_of_simplices            = 0;
	this->sub_trees_in_queue         = 0;
	this->add_sub_trees_total        = 0;
	this->no_of_sp                   = 0;
	this->no_of_sp_post_optimisation = 0;
	this->no_of_sp_overlapping       = 0;
}
