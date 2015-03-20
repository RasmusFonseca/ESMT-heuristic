#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "test/test.hpp"
#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/esmt.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/heuristics/concat.hpp"
#include "steiner/heuristics/steiner_finder.hpp"
#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/delaunay.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point_set_generator.hpp"

typedef Utils::Point      Point;
typedef Utils::Edge       Edge;
typedef Utils::Generator  Generator;
typedef Test::Result      Result;

Test::Test()
  : seed(0), max(100), min(-100),
    loop_time(2), concat(false), po(true), sct(false),
    collect_stats(false), do_delaunay(true),
    sh(NULL), mean_ratio(0.0), mean_time(0.0), no_of_bad_trees(0),
    no_of_suboptimal_trees(0)
{
}

Test::~Test() {}

void Test::setSeed(int seed) {
  this->seed = seed;
  Generator::setSeed(seed);
}

void Test::addTest(int point_set, int dim, int no_of_points, int no_of_tests) {
  if(this->seed==0)
    this->seed = time(NULL);
  this->addPointSet(point_set, dim, no_of_points, no_of_tests);
}

void Test::addTest(std::string &filename, std::string &setname) {
  if(this->seed==0)
    this->seed = time(NULL);
  this->addPointSet(filename, setname);
}

void Test::addTest(std::string &filename) {
  if(this->seed==0)
    this->seed = time(NULL);
  this->addPointSet(filename);
}

void Test::doConcatSubgraphs(bool concat) {
  this->concat = concat;
}

void Test::doPostOptimise(bool po) {
  this->po = po;
}

void Test::doUseSpecialConcat(bool sct) {
  this->sct = sct;
}

void Test::doCollectStats(bool doCollect) {
  this->collect_stats = doCollect;
}

void Test::inclDelaunay(bool doDelaunay) {
  this->do_delaunay = doDelaunay;
}

void Test::setLoopTime(int sec) {
  this->loop_time = sec;
}

void Test::setSubgraphHeuristicOne(SubgraphHeuristic *sh, std::string &name) {
  this->sh_name = name;
  this->sh = sh;
}

void Test::clear() {
  this->results.clear();
  this->sets.clear();
}

std::vector<Result> *Test::getResults() {
  return &this->results;
}

void Test::doTestSubgraphAlgorithm(bool verbose) {
  assert(this->sh);

  std::string text = "Test of subgraph algorithm: " + this->sh_name;
  printHeader(text);

  std::cout << std::endl
	    << "  Algortihm: " << this->sh_name << std::endl
	    << "  No of tests: " << this->sets.size() << std::endl
	    << "  Seed: "        << this->seed << std::endl
	    << std::endl << " * Running tests *" << std::endl;;
  
  this->mean_ratio = 0;
  this->mean_time  = 0;
  this->mean_time_deviation = 0.0;
  this->no_of_bad_trees = 0;
  this->no_of_suboptimal_trees = 0;
  
  int no_of_tests = this->sets.size();
  for(int i = 0; i < no_of_tests; i++) {
    this->subgraphAlgTest(i, verbose);
  }
  
  this->mean_ratio /= no_of_tests;
  this->mean_time  /= no_of_tests;
  this->calcQuartilesRatio();

  std::cout << " * Tests done *" << std::endl
	    << std::endl
	    << "  " << this->sh_name << std::endl
	    << "    Avg time: " << this->mean_time << std::endl
	    << "    Time std. dev.: " << this->calcStdDevTime() << std::endl
	    << "    Ratio mean: " << this->mean_ratio << std::endl
	    << "    Ratio std. dev.: " << this->calcStdDevRatio() << std::endl
	    << "       min: " << this->quartiles[0] << std::endl
	    << "       1st quartile: " << this->quartiles[1] << std::endl
	    << "       median: " << this->quartiles[2] << std::endl
	    << "       3rd quartile: " << this->quartiles[3] << std::endl
	    << "       max: " << this->quartiles[4] << std::endl
    	    << "    No of bad trees: " << this->no_of_bad_trees << " ("
    	    << (100.0*((double)this->no_of_bad_trees)/((double)no_of_tests))
	    << "%)" << std::endl;
  
  std::cout << std::endl;
  text = "End of test";
  printHeader(text);  
}

void Test::doTestESMT(bool verbose) {
  assert(this->sh);

  std::string text = "Test of ESMT heuristic";
  printHeader(text);

  std::cout << std::endl
	    << "  Subgraph algortihm: " << this->sh_name << std::endl
	    << "  Concat: "             << this->concat << std::endl
	    << "  Post-Optimisation: "  << this->po << std::endl
	    << "  No of tests: " << this->sets.size() << std::endl
	    << "  Seed: "        << this->seed << std::endl
	    << std::endl << " * Running tests *" << std::endl;
  
  this->mean_ratio       = 0;
  this->mean_time        = 0;
  this->mean_time_deviation = 0;
  this->no_of_bad_trees     = 0;


  int no_of_tests = this->sets.size();
  for(int i = 0; i < no_of_tests; i++)
    this->ESMTTest(i, verbose, no_of_tests==1);

  this->mean_ratio       /= no_of_tests;
  this->mean_time        /= no_of_tests;
  this->mean_time_deviation /= no_of_tests;
  this->calcQuartilesRatio();
  
  std::cout << " * Tests done *" << std::endl
	    << std::endl
	    << "  " << this->sh_name << std::endl
	    << "    Avg time: " << this->mean_time << std::endl
	    << "    Time std. dev.: " << this->calcStdDevTime() << std::endl
	    << "    Ratio mean: " << this->mean_ratio << std::endl
	    << "    Ratio std. dev.: " << this->calcStdDevRatio() << std::endl
	    << "       min: " << this->quartiles[0] << std::endl
	    << "       1st quartile: " << this->quartiles[1] << std::endl
	    << "       median: " << this->quartiles[2] << std::endl
	    << "       3rd quartile: " << this->quartiles[3] << std::endl
	    << "       max: " << this->quartiles[4] << std::endl;
  if(this->no_of_bad_trees > 0) {
    std::cout << "    No of bad trees: " << this->no_of_bad_trees << " ("
	      << (100.0*((double)this->no_of_bad_trees)/((double)no_of_tests))
	      << "%)" << std::endl;
  }
  std::cout << std::endl;

  text = "End of test";
  printHeader(text);  
}

void Test::doTestESMTSpecial(int d, int n, int seed, bool verbose) {
  unsigned int i, j, sum;
  assert(this->sh);

  std::string text = "Special test of ESMT heuristic";
  printHeader(text);

  std::cout << std::endl << " * Running test *" << std::endl;

  Generator::setSeed(seed);
  IterativeConcat ic(d);
  std::vector<Point> points = Generator::randomFloatPoints(Point(d,100), Point(d,100), n);
  Utils::Delaunay del(points);
  ESMT esmt(del, &ic, true, true, false, verbose);

  std::vector<unsigned int> faces = del.getNumberOfFaces();

  std::cout << std::endl << " * Test done *" << std::endl << std::endl
	    << "  Delaunay faces: " << std::endl
	    << "K;No;Accum" << std::endl;
  sum = 0;
  for(i = 2; i < faces.size(); i++) {
    sum += faces[i];
    std::cout << i << ";" << faces[i] << ";" << sum << std::endl;
  }
  
  std::cout << std::endl << "  Covered faces (and sausages): " << std::endl
	    << "K;No;Accum;Frac" << std::endl;
  sum = 0;
  for(j = 0; j < esmt.getStats()->covered_faces.size(); j++) {
    unsigned int no = esmt.getStats()->covered_faces[j];
    if(no != 0) {
      sum += no;
      std::cout << j << ";" << no << ";" << sum << ";";
      if(j < faces.size())
	std::cout << (float)no / (float)faces[j];
      else 
	std::cout << "-";
      std::cout << std::endl;
    }
  }
  
  std::cout << std::endl << "  Added sub-trees:" << std::endl
	    << "K;No;Accum" << std::endl;
  sum = 0;
  for(j = 0; j < esmt.getStats()->added_sub_trees.size(); j++) {
    unsigned int no = esmt.getStats()->added_sub_trees[j];
    if(no != 0) {
      sum += no;
      std::cout << j << ";" << no << ";" << sum << std::endl;
    }
  }

  std::cout << std::endl << "  Steiner points:" << std::endl
	    << "Before;After;Added;Overlapping" << std::endl;
  std::cout << esmt.getStats()->no_of_sp << ";"
	    << esmt.getStats()->no_of_sp_post_optimisation << ";"
	    << (esmt.getStats()->no_of_sp_post_optimisation-esmt.getStats()->no_of_sp)
	    << ";" << esmt.getStats()->no_of_sp_overlapping << std::endl;
  
  std::cout << std::endl;

  text = "End of test";
  printHeader(text); 
}

void Test::subgraphAlgTest(int i, bool verbose) {
  double start_time, end_time;
  int iterations = 0;
  std::vector<Point> points = this->sets[i].points;

  // Make all possible edges. Needed for mst calculation
  std::vector<Edge> edges;
  for(int j = 0; j < (int)points.size(); j++)
    for(int k = 0; k < j; k++)
      edges.push_back(Edge(j,k,Utils::length(points[j],points[k])));
  Graph g(points, edges);
  g = Utils::MSTKruskal(g);
  SteinerTree *st = new SteinerTree();
  
  start_time = getTime();
  do {
    delete st;
    st = this->sh->findSteinerPoints(g);
    iterations++;
  }
  while((end_time = getTime()) < start_time + this->loop_time);
  
  Result res;
  res.time          = ((double)(end_time-start_time)) / iterations;
  res.ratio         = st->getSteinerRatio();
  this->mean_time  += res.time;
  this->mean_ratio += res.ratio;
  if(res.ratio*0.999999 > 1.0)
      this->no_of_bad_trees++;

  this->results.push_back(res);
  
  delete st;

  if(verbose)
    // Print result
    std::cout << " [ Test " << (i+1) << " | Time: " << res.time
	      << " | Ratio: " << res.ratio << " ]" << std::endl;
}

void Test::ESMTTest(int i, bool verbose, bool method_verbose) {
  double start_time, end_time;
  int iterations = 0;
  std::vector<Point> points = this->sets[i].points;
  ESMT *esmt = NULL;
  std::vector<unsigned int> faces;
  if(!this->do_delaunay) {
    Utils::Delaunay del(points);
    start_time = getTime();
    do {
      if(esmt)
	delete esmt;
      esmt = new ESMT(del, this->sh, this->concat, this->po, this->sct, method_verbose);
      iterations++;
    }
    while((end_time = getTime()) < start_time + this->loop_time);
  }
  else {
    start_time = getTime();
    do {
      if(esmt)
	delete esmt;
      esmt = new ESMT(points, this->sh, this->concat, this->po, this->sct, method_verbose);
      iterations++;
      method_verbose = false;
    }
    while((end_time = getTime()) < start_time + this->loop_time);
  }
  Result res;
  res.time          = ((double)(end_time-start_time)) / iterations;
  res.ratio         = esmt->getSteinerRatio();
  res.smt           = esmt->getSMTLength();
  res.mst           = esmt->getMSTLength();
  if(this->collect_stats)
    res.stat = *esmt->getStats();
  this->mean_time  += res.time;
  this->mean_ratio += res.ratio;
  if(res.ratio*0.99999 > 1.0)
      this->no_of_bad_trees++;

  this->results.push_back(res);
  delete esmt;

  // Print result
  if(verbose)
    std::cout << " [ Test " << (i+1) << " | Time: " << res.time
	      << " | Ratio: " << res.ratio << " ]" << std::endl;
}

/*
 * Helper function for generating point sets for tests.
 *
 * Generates no_of_tests sets, with no_of_points points in each set,
 * of type point_set.
 */
void Test::addPointSet(int point_set, int dim, int no_of_points, int no_of_tests = 1) {
  unsigned int i, j;
  Point max(dim, this->max);
  Point min(dim, this->min);

  if(point_set==TEST_POINT_SET_DELAUNAY_SIMPLICES) {
    //Do Delaunay
    std::vector<Point> points = Generator::randomFloatPoints(max, min, no_of_tests+dim);
    Utils::Delaunay del(points);
    for(i = 0; (int)i < no_of_tests && i < del.getSimplices()->size(); i++) {
      Set set;
      set.name = "Simplex from Delaunay";
      for(j = 0; j < (unsigned int)dim+1; j++) {
	Utils::Delaunay::Simplex s = (*del.getSimplices())[i];
	set.points.push_back(points[s.map[j]]);
      }
      this->sets.push_back(set);
    }
  }      
  else {
    for(i = 0; (int)i < no_of_tests; i++) {
      Set set;
      switch(point_set) {
      case TEST_POINT_SET_RANDOM_I:
	set.name   = "Random integer point set";
	set.points = Generator::randomIntPoints(max, min, no_of_points);
	break;
      case TEST_POINT_SET_RANDOM_D:
	set.name   = "Random double point set";
	set.points = Generator::randomFloatPoints(max, min, no_of_points);
	break;
      case TEST_POINT_SET_GRID:
	set.name = "Points in a grid";
	set.points = Generator::gridFromSide(dim, no_of_points);
	break;
      case TEST_POINT_SET_SAUSAGE:
	set.name   = "Smith/Du sausage";
	set.points = Generator::sausage(dim, no_of_points);
	break;
      case TEST_POINT_SET_SIMPLEX:
	set.name   = "Regular d-simplex";
	set.points = Generator::simplex(dim);
	break;
      default:
	std::cerr << "Unknown point set: " << point_set << std::endl;
	exit(1);
      }
      this->sets.push_back(set);
    }
  }
}

/*
 * Helper function for loading point sets for tests.
 */
void Test::addPointSet(std::string &filename, std::string &setname) {
  Set set;
  set.name   = setname;
  set.points = Generator::loadFromFile(filename, setname);
  this->sets.push_back(set);
}

/*
 * Helper for adding all sets from one file.
 */
void Test::addPointSet(std::string &filename) {
  std::vector < std::vector<Point> > loaded_sets;
  std::vector < std::string > names;
  unsigned i;
  loaded_sets = Generator::loadFromFile(filename, names);
  
  for(i = 0; i < loaded_sets.size(); i++) {
    Set set;
    set.name   = names[i];
    set.points = loaded_sets[i];
    this->sets.push_back(set);
  }
}

/* Helper function for printing a nice header/footer */
void Test::printHeader(const std::string &text) {
  int text_size = text.size();
  
  int width = text_size+4 > 80 ? text_size+4 : 80;
  
  for(int i = 0; i < width; i++)
    std::cout << "#";
  std::cout << std::endl << "# ";

  int indent = (width - text_size) / 2 - 2;
  for(int i = 0; i < indent; i++)
    std::cout << " ";

  std::cout << text;
  if(text_size % 2 == 1)
    std::cout << " ";

  for(int i = 0; i < indent; i++)
    std::cout << " ";

  std::cout << " #" << std::endl;

  for(int i = 0; i < width; i++)
    std::cout << "#";
  std::cout << std::endl;
}

/*
 * Gets the time since program start in seconds
 */
int Test::getTime() {
  return clock() / CLOCKS_PER_SEC;
}

/*
 * Prints a .dat file for use with pgfplots
 */
void Test::createDatFile(const std::string &fileName) {
  std::ofstream file;
  file.open (fileName.c_str());
  std::string sep = ";";
  int no_of_tests = this->results.size();
  file << "% Test results - ESMT heuristic" << std::endl
       << "% Seed: " << this->seed << std::endl
       << "% Subgraph alg: " << this->sh_name << std::endl
       << "% Loop_time = " << this->loop_time << std::endl
       << "% Concat sub-graphs: " << this->concat << std::endl
       << "% Post optimisation: " << this->po << std::endl
       << "% Special concat: " << this->sct << std::endl
       << "% No of tests: " << no_of_tests << std::endl
       << "% Test results: " << std::endl;
  
  file << "no" << sep << "name" << sep << "N"
       << sep << "d" << sep << "t" << sep << "l_mst" << sep
       << "l_smt" << sep << "rho";
  if(this->collect_stats)
    file << sep << "Del_Simplices" << sep << "Sub_Trees" << sep
	 << "Sub_Trees_Add" << sep << "No_Of_Sp" << sep << "No_Of_Sp_Po";
  file << std::endl;
  for(int i = 0; i < no_of_tests; i++) {
    file << (i+1) << sep << this->sets[i].name << sep
	 << this->sets[i].points.size() << sep
	 << this->sets[i].dim() << sep
	 << this->results[i].time << sep
	 << this->results[i].mst << sep
	 << this->results[i].smt << sep
	 << this->results[i].ratio;
  if(this->collect_stats)
    file << sep << this->results[i].stat.no_of_simplices
	 << sep << this->results[i].stat.sub_trees_in_queue
	 << sep << this->results[i].stat.add_sub_trees_total
	 << sep << this->results[i].stat.no_of_sp
	 << sep << this->results[i].stat.no_of_sp_post_optimisation;
    file << std::endl;
  }
  double stddevrat  = this->calcStdDevRatio();
  double stddevtime = this->calcStdDevTime();
  file << "% Mean ratio: " << this->mean_ratio << std::endl
       << "% Ratio std. dev.: " << stddevrat << " (" 
       << (this->mean_ratio-stddevrat) << ", "
       << (this->mean_ratio+stddevrat) << ")" << std::endl
       << "% Mean time: " << this->mean_time << std::endl
       << "% Time std. dev.: " << stddevtime << " (" 
       << (this->mean_time-stddevtime) << ", "
       << (this->mean_time+stddevtime) << ")" << std::endl
       << "% Quartiles: " << std::endl
       << "%    min: " << this->quartiles[0] << std::endl
       << "%    1st quartile: " << this->quartiles[1] << std::endl
       << "%    median: " << this->quartiles[2] << std::endl
       << "%    3rd quartile: " << this->quartiles[3] << std::endl
       << "%    max: " << this->quartiles[4] << std::endl;
  file.close();
  std::cout << " File " << fileName << " created!" << std::endl;
}

double Test::calcStdDevRatio() {
  double res = 0.0;
  unsigned int i;
  for(i = 0; i < this->results.size(); i++) {
    double tmp = this->mean_ratio - this->results[i].ratio;
    res += tmp*tmp;
  }
  res /= (double)(this->results.size()-1);
  res = std::sqrt(res);
  return res;
}

double Test::calcStdDevTime() {
  double res = 0.0;
  unsigned int i;
  for(i = 0; i < this->results.size(); i++) {
    double tmp = this->mean_time - this->results[i].time;
    res += tmp*tmp;
  }
  res /= (double)(this->results.size()-1);
  res = std::sqrt(res);
  return res;
}

void Test::calcQuartilesRatio() {
  std::vector<Result> s = this->results;
  std::sort(s.begin(), s.end(), Test::compareRatio);
  unsigned int m1, m2;
  this->quartiles[0] = s[0].ratio;
  this->quartiles[4] = s[s.size()-1].ratio;
  this->quartiles[2] = this->calcMedian(s, 0, s.size());
  m1 = s.size()/2;
  m2 = m1;
  if(s.size()%2==1)
    m2++;
  this->quartiles[1] = calcMedian(s, 0, m1);
  this->quartiles[3] = calcMedian(s, m2, s.size());
}

double Test::calcMedian(std::vector<Result> &s, int begin, int end) {
  int size = end-begin;
  if(size%2==1)
    return s[begin+size/2].ratio;
  else
    return (s[begin+size/2].ratio + s[begin+size/2-1].ratio) / 2.0;
}

bool Test::compareRatio(const Result &res1, const Result &res2) {
  return res1.ratio < res2.ratio;
}

bool Test::testTopology(SteinerTree &st1, SteinerTree &st2) {
  std::vector<Point> points1 = st1.getPoints();
  std::vector<Point> points2 = st2.getPoints();
  if(points1.size() != points2.size())
    return false;
  
  std::vector<Edge> edges1 = st1.getEdges();
  std::vector<Edge> edges2 = st2.getEdges();

  if(edges1.size() != edges2.size())
    return false;
  
  int n = 0, n1 = 0, n2 = 0, c1, c2;
  for(int i = 0; i < (int)points1.size(); i++) {
    c1 = 0;
    c2 = 0;
    for(int j = 0; j < (int)edges1.size(); j++) {
      if(edges1[j].i0==i || edges1[j].i1==i) {
	c1++;
	n1 = (edges1[j].i0==i?edges1[j].i1:edges1[j].i0);
      }
      if(edges2[j].i0==i || edges2[j].i1==i) {
	c2++;
	n2 = (edges2[j].i0==i?edges2[j].i1:edges2[j].i0);
      }
    }
    if(c1 != c2)
      return false;
    else if(c1 == 1) {
      n = i;
      break;
    }
  }
  return this->testTopoRec(st1, st2, n1, n2, n, n);
}

bool Test::testTopoRec(SteinerTree &st1, SteinerTree &st2,
		       int n1, int n2, int l1, int l2) {
  std::vector<Edge>  *e1 = st1.getEdgesPtr();
  std::vector<Edge>  *e2 = st2.getEdgesPtr();
  
  int e11, e12, e21, e22;
  e11 = -1;
  e12 = -1;
  e21 = -1;
  e22 = -1;

  int e = e1->size();

  for(int i = 0; i < e; i++) {
    if((*e1)[i].i0 == n1 && (*e1)[i].i1 != l1) {
      if(e11 < 0)
	e11 = (*e1)[i].i1;
      else
	e12 = (*e1)[i].i1;
    }
    else if((*e1)[i].i1 == n1 && (*e1)[i].i0 != l1) {
      if(e11 < 0)
	e11 = (*e1)[i].i0;
      else
	e12 = (*e1)[i].i0;
    }
    if((*e2)[i].i0 == n2 && (*e2)[i].i1 != l2) {
      if(e21 < 0)
	e21 = (*e2)[i].i1;
      else
	e22 = (*e2)[i].i1;
    }
    else if((*e2)[i].i1 == n2 && (*e2)[i].i0 != l2) {
      if(e21 < 0)
	e21 = (*e2)[i].i0;
      else
	e22 = (*e2)[i].i0;
    }
  }
  
  if(e11 < 0 && e21 < 0)
    return true;
  if((e11+1)*(e21+1)==0 || (e12+1)*(e22+1)==0)
    return false;
  
  if(e12 < 0) {
    // Only one exit.
    return this->testTopoRec(st1, st2, e11, e21, n1, n2);
  }
  else {
    // Two exits.
    bool b1, b2;
    b1 = this->testTopoRec(st1, st2, e11, e21, n1, n2);
    b2 = this->testTopoRec(st1, st2, e12, e22, n1, n2);
    if(b1 && b2)
      return true;
    // Permutation
    b1 = this->testTopoRec(st1, st2, e11, e22, n1, n2);
    b2 = this->testTopoRec(st1, st2, e12, e21, n1, n2);
    if(b1 && b2)
      return true;
  }
  return false;
}

int Test::Set::dim() {
  if(this->points.size() == 0)
    return -1;
  return this->points[0].dim();
}
