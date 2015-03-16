#ifndef POINT_SET_GENERATOR_H
#define POINT_SET_GENERATOR_H

#include <vector>

#include "steiner/utils/point.hpp"

typedef Utils::Point     Point;

/**
 * Namespace Utils.
 * Contains utility functions for the heuristic ESMT algorithm
 */
namespace Utils {
  
  /**
   * @class Generator
   * @authors Stephan Lorenzen and Andreas Olsen
   * 
   * Utility for generating different point sets, e.g. random sets and
   * various simplexes.
   */
  class Generator {
  public:
    
    /**
     * Set the seed for the generator
     *
     * @param seed   Seed for the PRNG
     */
    static void setSeed(int seed);
    
    /**
     * Get the seed for the generator
     *
     * @return       The seed used for the PRNG
     */
    static int getSeed();

    /**
     * Generates a set of randomly distributed points with integer coordinates
     *
     * The points are limited in the cube defined by the max/min values.
     *
     * @param max           Maximum value. max must be greater than min.
     * @param min           Minimum value.
     * @param no_of_points  Number of points to be generated.
     *
     * @return              A set of no_of_points random points (integer valued).
     */
    static std::vector<Point> randomIntPoints(Point max, Point min,
					      int no_of_points);
    
    
    /**
     * Generates a set of randomly distributed points with double coordinates
     *
     * The points are limited in the cube defined by the max/min values.
     *
     * @param max           Maximum value. max must be greater than min.
     * @param min           Minimum value.
     * @param no_of_points  Number of points to be generated.
     *
     * @return              A set of no_of_points random points (double valued).
     */
    static std::vector<Point> randomFloatPoints(Point max, Point min,
						int no_of_points);

    /**
     * Generates a grid of points inside the given limits
     *
     * @param dim          Dimension
     * @param no_of_points The number of points to generate. If the number
     *                     is not of the form a^dim, the number will be
     *                     rounded down to such a number 
     *
     * @return             A set of points arranged in a grid.
     */
    static std::vector<Point> grid(unsigned int dim, unsigned int no_of_points);
  
    /**
     * Generates a grid of points inside the given limits with the given
     * distance between each point.
     *
     * @param dim         Dimension
     * @param side        The side-length of the generated grid.
     *                    Total number of points generated is side^dim
     *
     * @return            A set of points arranged in a grid.
     */
    static std::vector<Point> gridFromSide(unsigned int dim, unsigned int side);
    
    /**
     * Generates the sausage point set described by Du and Smith in
     * 'Three disproofs of the Gilbert-Pollak conjecture on Steiner
     * ratio in three or more dimensions'.
     * 
     * The set is created as follows:
     *  1. Create a d-simplex around the origin.
     *  2. For the face with max indices, create a new simplex on that face.
     *  3. Repeat two in the same direction.
     *
     * @param d            The dimension
     * @param N            The number of points in the sausage
     *
     * @return             A sausage point set as described by Smith
     */
    static std::vector<Point> sausage(int d, int N);
    
    /**
     * Generates a d-dimensional regular simplex.
     * The simplex is generated with center in the origin with distance == 1
     * to each vertex.
     *
     * @param d    The dimension.
     *
     * @return     A point set of size d+1, describing a regular d-simplex.
     */
    static std::vector<Point> simplex(int d);

    /**
     * Load a point set from a STP formatted file.
     *
     * @param file     The path to the STP formatted file.
     * @param set_name The name of the point set in the file.
     * @param verbose  Verbose
     *
     * @return         The set of points read from the file.
     */
    static std::vector<Point> loadFromFile(std::string &file, std::string &set_name, bool verbose = false);

    /**
     * Load all point sets from a STP formatted file.
     *
     * @param file     The path to the STP formatted file.
     * @param names    A vector to be filled with names of loaded files.
     * @param verbose  Verbose
     *
     * @return         A set of point sets read from the file.
     */
    static std::vector< std::vector<Point> > loadFromFile(std::string &file, std::vector<std::string> &names, bool verbose = false);
  
  private:
    /** The seed for this Generator. */
    static unsigned int seed;

    /** Recursive function for generating d-dimensional grids */
    static void gridRec(Point &cur, std::vector<Point> &result, int dim, int side);
  };
}

#endif // POINT_SET_GENERATOR_H
