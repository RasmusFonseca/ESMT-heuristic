#ifndef UTILS_DELAUNAY_H
#define UTILS_DELAUNAY_H

#include <vector>

#include "steiner/graph.hpp"
#include "steiner/utils/point.hpp"

/*
 * Namespace Utils.
 * Contains utility functions for the heuristic ESMT algorithm
 */
namespace Utils {

  /**
   * @class Delaunay
   *
   * This is a simple interface for performing the Delaunay triangulation in
   * d dimensions.
   */
  class Delaunay : public Graph {
  public:
    /**
     * Creates the Delaunay tesselation of the given point set in d dimensions.
     * Using the QHULL package.
     *
     * @param points  The point set.
     */
    Delaunay(std::vector<Utils::Point> &points);
    
    /**
     * Destructor.
     */
    ~Delaunay();
    
    /**
     * @struct PointHandle
     *
     * For each point, stores a list of simplex indices.
     */
    struct PointHandle {
      std::vector<int> simplices;
    };

    /**
     * @struct Simplex
     *
     * Represents a simplex in the tesselation, containing d+1 points.
     *
     * The points a represented as a list of indices. Furthermore, indices
     * of neighbouring simplices and vertices.
     */
    struct Simplex {
      /**
       * Constructor.
       *
       * Allocates the structure.
       *
       * @param d   The dimension.
       */
      Simplex(int d);
      
      /**
       * Copy constructor.
       *
       * @param s   The simplex to copy.
       */
      Simplex(const Simplex &s);
      
      /**
       * Destructor.
       */
      ~Simplex();
      
      unsigned int n;
      int *map;
      int *nextSimplex;
      int *nextVertex;
    };

    /**
     * Getter for the list of simplices.
     */
    std::vector<Simplex>     *getSimplices();
    
    /**
     * Getter for the list of point handles.
     */
    std::vector<PointHandle> *getPointHandles();
    
    /**
     * Stats: Count all faces
     */
    std::vector<unsigned int> getNumberOfFaces();

  protected:
  private:
    /**
     * Performs the Delaunay tesselation using qhull.
     *
     * @param d   The dimension.
     */
    void doDelaunayQHull(int d);

    /**
     * Sub procedure for doDelaunayQHull.
     *
     * Creates an input file for the qhull program.
     *
     * @param d   The dimension.
     */
    void doPrepInputQHull(int d);
    
    /**
     * Sub procedure for doDelaunayQHull.
     *
     * Parses the output file from the qhull program.
     *
     * @param d   The dimension.
     */
    void doParseOutputQHull(int d);

    /**
     * Recursive procedure for finding all faces
     */
    void findFaces(Simplex &simplex, std::vector<unsigned int> &cur_set,
		   std::unordered_map<unsigned long, bool> &flag,
		   std::vector<unsigned int> &result);

    /** The list of simplices */
    std::vector<Simplex>      simplices;
    /** The list of point handles */
    std::vector<PointHandle>  point_handles;

	/** qdelaunay input file name */
	std::string qdelaunay_input;

	/** qdelaunay output file name */
	std::string qdelaunay_output;
  };

}

#endif // UTILS_DELAUNAY_H
