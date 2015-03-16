#include <vector>
#include <cstdlib>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "steiner/utils/point.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point_set_generator.hpp"

// States for STP parse operation
#define PSG_PARSE_STP_INITIAL            0
#define PSG_PARSE_STP_SEC_COMMENTS_WAIT  1
#define PSG_PARSE_STP_SEC_COMMENTS       2
#define PSG_PARSE_STP_SEC_GRAPH_WAIT     3
#define PSG_PARSE_STP_SEC_GRAPH          4
#define PSG_PARSE_STP_SEC_COORDS_WAIT    5
#define PSG_PARSE_STP_SEC_COORDS         6
#define PSG_PARSE_STP_EOF_WAIT           7

typedef Utils::Point     Point;

/* Seed holder */
unsigned int Utils::Generator::seed = 0;

/* Implementation of setSeed(...) */
void Utils::Generator::setSeed(int seed) {
  Utils::Generator::seed = seed;
  srand(seed);
}

/* Implementation of getSeed() */
int Utils::Generator::getSeed() {
  return Utils::Generator::seed;
}

/* Implementation of randomIntPoints(...) */
std::vector<Point> Utils::Generator::randomIntPoints(Point max, Point min,
						     int no_of_points) {
  int dimension = (int)max.dim();
  // Validate
  assert(dimension == (int)min.dim());
  for (int i = 0; i < dimension; i++)
    assert(max[i] > min[i]);
  
  // Generate
  std::vector<Point> result;
  for(int i = 0; i < no_of_points; i++) {
    std::vector<double> list;
    for (int j = 0; j < dimension; j++) {
      list.push_back((int)(rand() % (int)(max[j]-min[j]+1)) + min[j]);
    }
    result.push_back(Point(list));
  }
  return result;
}
    
/* Implementation of randomFloatPoints(...) */
std::vector<Point> Utils::Generator::randomFloatPoints(Point max, Point min,
						       int no_of_points) {
  int dimension = (int)max.dim();
  // Validate
  assert(dimension == (int)min.dim());
  for (int i = 0; i < dimension; i++)
    assert(max[i] >= min[i]);
  
  // Generate
  std::vector<Point> result;
  for(int i = 0; i < no_of_points; i++) {
    // Is this uniformly distributed? Probably not...
    std::vector<double> list;
    for (int j = 0; j < dimension; j++) {
      list.push_back((double)(rand() % (int)(max[j]-min[j]+1) + min[j]) + Utils::frand());
    }
    result.push_back(Point(list));
  }

  return result;
}

/* Implementation of grid(...) (no_of_points version) */
std::vector<Point> Utils::Generator::grid(unsigned int dim,
					  unsigned int no_of_points) {
  double side = pow(no_of_points, 1.0 / (double)dim);
  return Utils::Generator::gridFromSide(dim, (int)side);
}
  
/* Implementation of grid(...) */
std::vector<Point> Utils::Generator::gridFromSide(unsigned int dim,
						  unsigned int side) {
  unsigned int i;
  // Validate
  int no_of_points = 1;
  for (i = 0; i < dim; i++)
    no_of_points *= side;
  assert(no_of_points < 20000);

  std::vector<Point> result;
  Point p(dim);
  Utils::Generator::gridRec(p, result, dim, side);
  return result;
}

/* Implementation of gridRec(...) */
void Utils::Generator::gridRec(Point &cur, std::vector<Point> &result,
			       int dim, int side) {
  if(dim==0) {
    result.push_back(cur);
    return;
  }
  dim--;
  double k;
  for(k = 0; k < side; k++) {
    cur[dim] = k;
    Utils::Generator::gridRec(cur, result, dim, side);
  }
}
    
/* Implementation of sausage(...) */
std::vector<Point> Utils::Generator::sausage(int d, int N) {
  int i;
  std::vector<Point> result;
  if(N <= d) {
    std::cerr << "Error when generating sausage. Bad argument, N = "
	      << N << " <= " << d << " = d" << std::endl;
    return result;
  }
  result = Utils::Generator::simplex(d);
  N -= (d+1);
  while(N > 0) {
    // Calculate centroid of last d points
    Point c(d,0);
    for(i = 0; i < d; i++) {
      c += result[result.size()-1-i];
    }
    c *= 1.0/(double)d;
    // Move the d+1 last point through the centroid.
    Point p = result[result.size()-1-d];
    Point v = c-p;
    p = c + v;
    result.push_back(p);
    N--;
  }
  return result;
}

/* Implementation of simplex(d) */
std::vector<Point> Utils::Generator::simplex(int d) {
  int i, j;
  double val, dot;
  std::vector<Point> result;
  dot = -(1.0/(double)d);
  for(i = 0; i < d+1; i++)
    result.push_back(Point(d));
  val = 1.0;
  for(i = 0; i < d-1; i++) {
    // i'th is set to val, rest to zero
    result[i][i] = val;
    for(j = i+1; j < d; j++)
      result[i][j] = 0.0;
    // Calculate next entries from dot.
    // dot == r[i][0]*r[i+1][0]...r[i][i]*x == -1/d
    // <=> x == (-1/d + r[i][0]*r[i+1][0]...r[i][i-1]*r[i+1][i-1]) / r[i][i]
    val = dot;
    for(j = 0; j < i; j++)
      val -= result[i][j]*result[i+1][j];
    val /= result[i][i];
    for(j = i+1; j < d+1; j++)
      result[j][i] = val;
    // Calculate val == r[i+1][i+1]
    // 1 == r[i+1][0]^2+...+r[i+1][i]^2+x^2
    // <=> x = sqrt(1-(r[i+1][0]^2+...+r[i+1][i]^2))
    val = 1.0;
    for(j = 0; j <= i; j++)
      val -= result[i+1][j]*result[i+1][j];
    val = std::sqrt(val);
  }
  result[d-1][d-1] = val;
  result[d][d-1]   = -val;
  return result;
}

/* Implementation of loadFromFile */
std::vector<Point> Utils::Generator::loadFromFile(std::string &path,
						  std::string &set_name,
						  bool verbose) {
  std::ifstream file(path.c_str());
  std::string buffer, comments, set_name_quotes;
  std::vector<Point> result;
  int state       = PSG_PARSE_STP_INITIAL;
  bool done       = false;
  bool found      = false;
  int no_of_nodes = 0;
  set_name_quotes = "\"" + set_name + "\"";
  if(file.is_open()) {
    while(getline(file, buffer) && !done) {
      std::istringstream iss(buffer);
      std::vector<std::string> tokens;
      std::copy(std::istream_iterator<std::string>(iss),
		std::istream_iterator<std::string>(),
		std::back_inserter< std::vector<std::string> >(tokens));
      if(tokens.size() == 0)
	// Empty line
	continue;
      switch(state) {
      case PSG_PARSE_STP_INITIAL:
	// Check file magic number
	if(tokens[0] != "33D32945") {
	  std::cerr << "Utils::Generator::loadFromFile(...) - Not a STP file. "
		    << "Magic number is wrong: " << tokens[0]
		    << ". Exiting!" << std::endl;
	  exit(1);
	}
	else
	  state = PSG_PARSE_STP_SEC_COMMENTS_WAIT;
	break;
      case PSG_PARSE_STP_SEC_COMMENTS_WAIT:
	if(tokens.size() != 2
	   || tokens[0] != "SECTION"
	   || tokens[1] != "Comments") {
	  std::cerr << "Utils::Generator::loadFromFile(...) - Format error. "
		    << "Expected comments section, but read: " << std::endl
		    << buffer << std::endl
		    << "Exiting!" << std::endl;
	  exit(1);
	}
	else
	  state = PSG_PARSE_STP_SEC_COMMENTS;
	break;
      case PSG_PARSE_STP_SEC_COMMENTS:
	// First entry in comments should be name
	if(tokens[0] == "END") {
	  state = PSG_PARSE_STP_SEC_GRAPH_WAIT;
	}
	else if(tokens.size() < 2) {
	  std::cerr << "Utils::Generator::loadFromFile(...) - Format error. "
		    << "Expected comment tag, but read: " << std::endl
		    << buffer << std::endl
		    << "Exiting!" << std::endl;
	  exit(1);
	}
	else if(tokens[0] == "Name" && tokens[1] == set_name_quotes) {
	  // Found the set
	  found = true;
	  comments += buffer + "\n";
	}
	else if(found) {
	  // A comment line
	  comments += buffer + "\n";
	}
	// This is not the correct set, simply skip
	break;
      case PSG_PARSE_STP_SEC_GRAPH_WAIT:
	if(tokens.size() != 2
	   || tokens[0] != "SECTION"
	   || tokens[1] != "Graph") {
	  std::cerr << "Utils::Generator::loadFromFile(...) - Format error. "
		    << "Expected graph section, but read: " << std::endl
		    << buffer << std::endl
		    << "Exiting!" << std::endl;
	  exit(1);
	}
	else
	  state = PSG_PARSE_STP_SEC_GRAPH;
	break;
      case PSG_PARSE_STP_SEC_GRAPH:
	if(tokens[0] == "END") {
	  state = PSG_PARSE_STP_SEC_COORDS_WAIT;
	}
	else if(found) {
	  if(tokens.size() != 2 || tokens[0] != "Nodes") {
	    std::cerr << "Utils::Generator::loadFromFile(...) - Format error. "
		      << "Expected 'Nodes' tag, but read: " << std::endl
		      << buffer << std::endl
		      << "Exiting!" << std::endl;
	    exit(1);
	  }
	  else {
	    std::istringstream(tokens[1]) >> no_of_nodes;
	  }
	}
	// Not the right set - skip
	break;
      case PSG_PARSE_STP_SEC_COORDS_WAIT:
	if(tokens.size() != 2
	   || tokens[0] != "SECTION"
	   || tokens[1] != "Coordinates") {
	  std::cerr << "Utils::Generator::loadFromFile(...) - Format error. "
		    << "Expected coordinates section, but read: " << std::endl
		    << buffer << std::endl
		    << "Exiting!" << std::endl;
	  exit(1);
	}
	else
	  state = PSG_PARSE_STP_SEC_COORDS;
	break;
      case PSG_PARSE_STP_SEC_COORDS:
	if(tokens[0] == "END") {
	  state = PSG_PARSE_STP_EOF_WAIT;
	}
	else if(found) {
	  // Syntax: DD   index x y
	  //         DDD  index x y z
	  //         DDDD ...
	  unsigned int dim = tokens[0].size();
	  if(dim != tokens.size()-2) {
	    // Error
	    std::cerr << "Utils::Generator::loadFromFile(...) - Format error. "
		      << "Wrong number of tokens for point, read: " << std::endl
		      << buffer << std::endl
		      << "Exiting!" << std::endl;
	    exit(1);
	  }
	  Point p(dim);
	  for(unsigned int i=0; i<dim; i++) {
	    std::istringstream(tokens[i+2]) >> p[i]; 
	  }
	  result.push_back(p);
	}
	// Not the right set, skip
	break;
      case PSG_PARSE_STP_EOF_WAIT:
	if(tokens.size() != 1 || tokens[0] != "EOF") {
	  std::cerr << "Utils::Generator::loadFromFile(...) - Format error. "
		    << "Expected EOF token, but read: " << std::endl
		    << buffer << std::endl
		    << "Exiting!" << std::endl;
	  exit(1);
	}
	else if(found)
	  done = true;
	else
	  state = PSG_PARSE_STP_INITIAL;
	break;
      default:
	std::cerr << "Utils::Generator::loadFromFile(...) - unexpected state. "
		  << "Exiting!" << std::endl;
	exit(1);
      }
    }
    if(!found || !done) {
      std::cerr << "Utils::Generator::loadFromFile(...) - Set " << set_name
		<<" not found in file, or file incomplete... Exiting!"
		<< std::endl;
      exit(1);
    }
    else {
      // Print comments section
      if(verbose) {
	std::cout << "*** Read set: " << set_name_quotes << " ***"
		  << std::endl << comments;
	if(no_of_nodes != (int)result.size())
	  std::cout << "Warning: no_of_nodes != result.size() !!!" << std::endl;
	for(int i = 0; i < ((int)set_name_quotes.size() + 18); i++) {
	  std::cout << "*";
	}
	std::cout << std::endl;
      }
    }
  }
  else {
    std::cerr << "Utils::Generator::loadFromFile(...) - could not open file: "
	      << path << ", exiting!" << std::endl;
    exit(1);
  }
  // Close file
  file.close();
  return result;
}

std::vector< std::vector<Point> > Utils::Generator::loadFromFile(std::string &path, std::vector<std::string> &names, bool verbose) {
  std::ifstream file(path.c_str());
  std::string buffer;
  unsigned int i;

  names.clear();
  
  if(file.is_open()) {
    while(getline(file, buffer)) {
      std::istringstream iss(buffer);
      std::vector<std::string> tokens;
      std::copy(std::istream_iterator<std::string>(iss),
		std::istream_iterator<std::string>(),
		std::back_inserter< std::vector<std::string> >(tokens));
      if(tokens.size() != 2)
	// Wrong line
	continue;
      if(tokens[0]=="Name")
	names.push_back(tokens[1].substr(1,tokens[1].size()-2));
    }
  }
  file.close();
  

  std::vector< std::vector<Point> > result;
  for(i = 0; i < names.size(); i++)
    result.push_back(Utils::Generator::loadFromFile(path, names[i], verbose));

  return result;
}
  
