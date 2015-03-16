/**
 * Test of the gdl_window
 */
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "test/test.hpp"
#include "steiner/utils/delaunay.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/point_set_generator.hpp"
#include "steiner/graph.hpp"
#include "steiner/esmt.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/heuristics/steiner_finder.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/heuristics/concat.hpp"

#define MAX_DIM     12

typedef Utils::Point      Point;
typedef Utils::Edge       Edge;
typedef Utils::Generator  Generator;

void usage();
std::string getString(std::string tag, int argc, char *argv[], bool mandatory=false);
bool getBool(std::string tag, int argc, char *argv[]);
int getInt(std::string tag, int argc, char *argv[], int def=-1, bool mandatory=false);
std::vector<Point> getSet(int argc, char *argv[]);


int main(int argc, char *argv[]) {
  if(argc < 2)
    usage();
  if (std::string("testtopo").compare(0, 20, argv[1]) == 0) {
    std::cout << "Running topo test!" << std::endl;
    std::string name("Smith iterative");
    std::string resFile("/home/stephan/Desktop/result.dat");
    std::string path("../data/protein3D/W1.stp");
    std::string sname("W1");
    std::vector<Point> points1 = Generator::loadFromFile(path, sname);
    path = "../data/protein3D/W3.stp";
    sname = "W3";
    std::vector<Point> points2 = Generator::loadFromFile(path, sname);

    ESMT e1(points1);
    ESMT e2(points2);

    Test test;
    std::cout << test.testTopology(e1,e2) << std::endl;

    return 0;
  }
  else if(std::string("fulltest").compare(0, 20, argv[1]) == 0) {
    Test test;
    test.doTestESMTSpecial();
    return 0;
  }
  else {
    /*
     * See usage() below for options
     */
    std::string v("-v");
    std::string pt("-pt");
    std::string npo("-npo");
    std::string nsc("-nsc");
    std::string nd("-nd");
    std::string cs("-cs");
    std::string out("-out");
    std::string lt("-lt");
    std::string alg("-alg");
    std::string s("-s");
    std::string sct("-sct");
    std::string val("-val");
    bool verbose = getBool(v, argc, argv);
    bool print   = getBool(pt, argc, argv);
    int n, t, tests, dim;
    std::string algname, outfile, filename, setname;
    SubgraphHeuristic *sh;
    SubgraphHeuristic *shsp = NULL;
    algname = getString(alg, argc, argv);
    if(algname.size()==0 || algname.compare(0, 10, "NO") == 0)
      sh = new IterativeSmith(MAX_DIM);
    else if(algname.compare(0, 10, "RNO") == 0)
      sh = new IterativeConcat(MAX_DIM);
    else if(algname.compare(0, 10, "SP") == 0) {
      shsp = new IterativeConcat(MAX_DIM);
      sh = new SteinerFinder(shsp);
    }
    else {
      std::cerr << "Unknown sub-graph algorithm: " << algname << std::endl;
      usage();
    }
    if(std::string("test").compare(0, 20, argv[1]) == 0) {
      if(argc < 3)
	usage();
      Test test;
      test.setSubgraphHeuristicOne(sh, algname);
      test.doConcatSubgraphs(!getBool(nsc, argc, argv));
      test.doPostOptimise(!getBool(npo, argc, argv));
      test.doUseSpecialConcat(getBool(sct, argc, argv));
      test.inclDelaunay(!getBool(nd, argc, argv));
      test.setLoopTime(getInt(lt, argc, argv, 0));
      test.setSeed(getInt(s, argc, argv, time(NULL)));
      test.doCollectStats(getBool(cs, argc, argv));
      outfile = getString(out, argc, argv);
      // Add tests
      bool tests_ok = false;
      for(int i = 2; i < argc; i++) {
	if(std::string("-in").compare(0, 20, argv[i]) == 0) {
	  if(i >= argc-2) {
	    std::cerr << "Missing arguments for -in" << std::endl;
	    usage();
	  }
	  filename = argv[i+1];
	  setname = argv[i+2];
	  test.addTest(filename,setname);
	  tests_ok = true;
	}
	else if(std::string("-ina").compare(0, 20, argv[i]) == 0) {
	  if(i >= argc-1) {
	    std::cerr << "Missing arguments for -ina" << std::endl;
	    usage();
	  }
	  filename = argv[i+1];
	  test.addTest(filename);
	  tests_ok = true;
	}
	else if(std::string("-g").compare(0, 20, argv[i]) == 0
		|| std::string("-gn").compare(0, 20, argv[i]) == 0) {
	  if(i >= argc-3) {
	    std::cerr << "Missing arguments for -g" << std::endl;
	    usage();
	  }
	  setname = argv[i+1];
	  dim     = atoi(argv[i+2]);
	  n       = atoi(argv[i+3]);
	  if(setname.compare(0, 20, "random") == 0)
	    t = TEST_POINT_SET_RANDOM_D;
	  else if(setname.compare(0, 20, "simplex") == 0)
	    t = TEST_POINT_SET_SIMPLEX;
	  else if(setname.compare(0, 20, "sausage") == 0)
	    t = TEST_POINT_SET_SAUSAGE;
	  else if(setname.compare(0, 20, "grid") == 0)
	    t = TEST_POINT_SET_GRID;
	  else if(setname.compare(0, 20, "delaunay") == 0)
	    t = TEST_POINT_SET_DELAUNAY_SIMPLICES;
	  else {
	    std::cerr << "Unknown point set used with "
		      << argv[i] << ": " << setname << std::endl;
	    usage();
	  }
	  tests = 1;
	  if(std::string("-gn").compare(0, 20, argv[i]) == 0) {
	    if(i >= argc-4) {
	      std::cerr << "Missing arguments for -gn" << std::endl;
	      usage();
	    }
	    tests = atoi(argv[i+4]);
	  }
	  test.addTest(t, dim, n, tests);
	  tests_ok = true;
	}
      }
      if(!tests_ok) {
	std::cerr << "You must include at least one point set for testing (use -in, -ina, -g or -gn)." << std::endl;
	usage();
      }
      if(std::string("esmt").compare(0, 20, argv[2]) == 0)
	test.doTestESMT(verbose);
      else if(std::string("sgh").compare(0, 20, argv[2]) == 0)
	test.doTestSubgraphAlgorithm(verbose);
      else
	usage();
      if(outfile.size() > 0)
	test.createDatFile(outfile);
    }
    else if(std::string("esmt").compare(0, 20, argv[1]) == 0) {
      Generator::setSeed(getInt(s, argc, argv, time(NULL)));
      std::vector<Point> points = getSet(argc, argv);
      if(points.size() == 0) {
	std::cerr << "Empty or no point set given!" << std::endl;
	usage();
      }
      ESMT esmt(points, sh, !getBool(nsc, argc, argv),
		!getBool(npo, argc, argv), getBool(sct, argc, argv), verbose);
      
      if(getBool(val, argc, argv)) {
	if(Utils::validate(esmt))
	  std::cout << "Validate OK!" << std::endl;
	else
	  std::cout << "Validate ERROR!" << std::endl;
      }

      if(print)
	std::cout << "## RESULT ##" << std::endl << esmt << std::endl;
      else {
	std::cout << "Done!" << std::endl
		  << "  |MST| = " << esmt.getMSTLength() << std::endl
		  << "  |SMT| = " << esmt.getSMTLength() << std::endl
		  << "  Ratio = " << esmt.getSteinerRatio() << std::endl;
      }
    }
    else if(std::string("sgh").compare(0, 20, argv[1]) == 0) {
      std::cerr << "Not implemented yet." << std::endl;
      exit(1);
    }
    delete sh;
    delete shsp;
    return 0;
  }
  std::cout << "argument fell through, functionality not yet implemented!" << std::endl;
  return 1;
}

void usage() {
  std::cout << "Usage:" << std::endl
	    << " esmt [options] <points>" << std::endl
	    << " test esmt [options] <points>" << std::endl
	    << " test sgh [options] <points>" << std::endl
	    << std::endl
	    << " options:" << std::endl
	    << "   -v           Verbose" << std::endl
	    << "   -npo         Disable post optimisation" << std::endl
	    << "   -nsc         Disable sub-graph concatenation" << std::endl
	    << "   -sct         Redo concatenation, adding other, non-covered FSTs" << std::endl
	    << "   -alg <name>  Sub-graph algorithm:" << std::endl
	    << "                  NO  = Numerical optimisation" << std::endl
	    << "                  RNO = Restricted numerical optimisation" << std::endl
	    << "                  SP  = Simplex partitioning" << std::endl
	    << "                  Default NO" << std::endl
	    << "   -s           Seed for random generation of point sets." << std::endl
	    << "  # esmt only" << std::endl
	    << "   -pt          Print resulting tree" << std::endl
	    << "   -val         Validate the resulting tree" << std::endl
	    << "  # test only" << std::endl
	    << "   -nd          Do not include Delaunay tesselation in time measurement." << std::endl
	    << "   -cs          Collect extra statistics (no of simplices, etc.)." << std::endl
	    << "   -out <path>  Output file for test results (CSV format)." << std::endl
	    << "   -lt <sec>    Loop time" << std::endl
	    << std::endl
	    << " points:" << std::endl
	    << "   -in <file name> <set name>          Read set <set name>" << std::endl
	    << "                                       from file <file name>" << std::endl
	    << "   -g <set name> <dim> <no of points>  Generate a point set with"<<std::endl
	    << "                                       the given number of points in"<<std::endl
	    << "                                       the given dimension."<<std::endl
	    << "                                       Possible set names:" << std::endl
	    << "                                           random, sausage, simplex, grid" << std::endl
	    << " points (test only):" << std::endl
	    << "   -ina <file name>              Read all sets from file <file name>" << std::endl
	    << "   -gn <set name> <dim> <no of points> <no of tests>" << std::endl
	    << "                                 Generate the given number of" << std::endl
	    << "                                 point set with the given"<<std::endl
	    << "                                 number of points in the given"<<std::endl
	    << "                                 dimension."<<std::endl;
  exit(1);
}

std::string getString(std::string tag, int argc, char *argv[], bool mandatory) {
  for(int i = 2; i < argc; i++) {
    if(tag.compare(0, 20, argv[i]) == 0) {
      if(i == argc-1) {
	std::cerr << "Missing argument for " << tag << std::endl;
	usage();
      }
      return std::string(argv[i+1]);
    }
  }
  if(mandatory)  {
    std::cerr << "Missing argument " << tag << std::endl;
    usage();
  }
  return "";
}

bool getBool(std::string tag, int argc, char *argv[]) {
  for(int i = 2; i < argc; i++) {
    if(tag.compare(0, 20, argv[i]) == 0)
      return true;
  }
  return false;
}

int getInt(std::string tag, int argc, char *argv[], int def, bool mandatory) {
  for(int i = 2; i < argc; i++) {
    if(tag.compare(0, 20, argv[i]) == 0) {
      if(i == argc-1) {
	std::cerr << "Missing argument for " << tag << std::endl;
	usage();
      }
      return atoi(argv[i+1]);
    }
  }
  if(mandatory)  {
    std::cerr << "Missing argument " << tag << std::endl;
    usage();
  }
  return def;
}

std::vector<Point> getSet(int argc, char *argv[]) {
  std::vector<Point> result;
  std::string file, name;
  int dim, n;
  for(int i = 2; i < argc; i++) {
    if(std::string("-in").compare(0, 20, argv[i]) == 0) {
      if(i >= argc-2) {
	std::cerr << "Missing arguments for -in" << std::endl;
	usage();
      }
      file = argv[i+1];
      name = argv[i+2];
      result = Generator::loadFromFile(file,name);
    }
    else if(std::string("-g").compare(0, 20, argv[i]) == 0) {
      if(i >= argc-3) {
	std::cerr << "Missing arguments for -g" << std::endl;
	usage();
      }
      name = argv[i+1];
      dim  = atoi(argv[i+2]);
      n    = atoi(argv[i+3]);
      if(name.compare(0, 20, "random") == 0)
	result  = Generator::randomFloatPoints(Point(dim, 100), Point(dim, -100), n);
      else if(name.compare(0, 20, "simplex") == 0)
	result  = Generator::simplex(dim);
      else if(name.compare(0, 20, "sausage") == 0)
	result  = Generator::sausage(dim, n);
      else if(name.compare(0, 20, "grid") == 0)
	result  = Generator::gridFromSide(dim, n);
    }
  }
  return result;
}
