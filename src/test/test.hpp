/**
 * File with bulk-tests of the different heuristics.
 *
 * Options allow for different types of heuristics, point-sets, etc.
 *
 * All results are printed to std::cout.
 */
#ifndef STEINER_TEST_H
#define STEINER_TEST_H

#include <vector>

#include "steiner/steiner_tree.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/esmt.hpp"

#define TEST_POINT_SET_RANDOM_I           0
#define TEST_POINT_SET_RANDOM_D           1
#define TEST_POINT_SET_GRID               2
#define TEST_POINT_SET_SAUSAGE            3
#define TEST_POINT_SET_SIMPLEX            4
#define TEST_POINT_SET_DELAUNAY_SIMPLICES 5

typedef Utils::Point      Point;

class Test {
public:
  struct Result {
    double smt;
    double mst;
    double ratio;
    double time;
    ESMT::Stats stat;
  };

  Test();
  ~Test();
  
  void setSeed(int seed);
  void addTest(int point_set, int dim, int no_of_points, int no_of_tests = 1);
  void addTest(std::string &filename, std::string &setname);
  void addTest(std::string &filename);
  void doConcatSubgraphs(bool concat=true);
  void doPostOptimise(bool po=true);
  void doUseSpecialConcat(bool sct=true);
  void doCollectStats(bool doCollect=true);
  void inclDelaunay(bool doDelaunay);
  void setLoopTime(int sec);
  void setSubgraphHeuristicOne(SubgraphHeuristic *sh, std::string &name);

  void clear();

  std::vector<Result> *getResults();
  
  void doTestSubgraphAlgorithm(bool verbose = true);
  void doTestESMT(bool verbose = true);
  void doTestESMTSpecial(bool verbose = true);
  
  void createDatFile(const std::string &fileName);  

  bool testTopology(SteinerTree &st1, SteinerTree &st2);
protected:
private:
  void subgraphAlgTest(int i, bool verbose);
  void ESMTTest(int i, bool verbose, bool method_verbose=false);
  
  void addPointSet(int point_set, int dim, int no_of_points, int no_of_tests);
  void addPointSet(std::string &filename, std::string &setname);
  void addPointSet(std::string &filename);

  void printHeader(const std::string &text);
  int getTime();

  double calcStdDevTime();
  double calcStdDevRatio();
  void calcQuartilesRatio();
  double calcMedian(std::vector<Result> &s, int begin, int end);
  static bool compareRatio(const Result &res1, const Result &res2);

  bool testTopoRec(SteinerTree &st1, SteinerTree &st2,
		   int n1, int n2, int l1, int l2);

  /* Config values */
  int seed;
  double max;
  double min;
  int loop_time;
  bool concat;
  bool po;
  bool sct;
  bool collect_stats;
  bool do_delaunay;
  double quartiles[5];

  /* Subgraph heuristic*/
  std::string sh_name;
  SubgraphHeuristic *sh;
  
  /* Sets and results */
  struct Set {
    std::string name;
    std::vector<Point> points;
    int dim();
  };

  std::vector<Test::Set> sets;
  std::vector<Test::Result> results;

  /* Stats */
  double mean_ratio;
  double mean_time;

  double mean_time_deviation;

  int no_of_bad_trees;
  int no_of_suboptimal_trees;
};

#endif /* STEINER_TEST_H */
