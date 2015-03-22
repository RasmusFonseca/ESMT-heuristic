#ifndef POINT_H
#define POINT_H

#include <vector>
#include <iostream>

#define COMP_ERR 0.01

/*
 * Namespace Utils.
 */
namespace Utils {
  /**
   * @class Point
   *
   * Represents a three dimensional point.
   * Contains a flag, that indicates if the point is a Steiner point.
   */
  class Point {
  public:
    /**
     * Default constructor
     *
     * Creates a point in the origin of a three dimensional space
     */
    Point();

    /**
     * Constructor
     *
     * @param dimension      The dimension of the point
     * @param value          The value of the point
     *
     * Creates a point of the specified dimension with all values as the specified value.
     */
    Point(int dimension, double value = 0);

    /**
     * Constructor.
     *
     * Creates a point (x,y,z), and sets the isSteiner
     * flag to isSteiner.
     * 
     * @param x           The x-coordinate
     * @param y           The y-coordinate
     * @param z           The z-coordinate
     * @param isSteiner   The initial value of the isSteiner flag.
     */
    Point(double x, double y, double z, bool isSteiner=false);
    
    /**
     * Constructor
     * 
     * Creates a point of dimension list.size
     * and sets the isStiner flag to isSteiner
     *
     * @param list       The list of coordinates
     * @param isSteiner  The initial value of the isSteiner flag
     */
    Point(std::vector<double> &list, bool isSteiner = false);

    /**
     * Copy constructor (=)
     *
     * @param other  The point to copy.
     */
    Point(Point const& other);
    
    /**
     * Destructor
     */
    ~Point();
    
    /**
     * Setter for the isSteiner flag.
     */
    void setSteiner(bool isSteiner = true);
    
    /**
     * Getter for the coordinate list.
     * 
     * @return The coordinate list.
     */
    double *getList() const;

    /** Getter for a pointer to the coordinate list.
     *
     * @return pointer to the coordinate list.
     */
    //double **getListPtr() const;

    /**
     * Getter for the isSteiner flag
     *
     * @return true if point is a Steiner point,
     *         otherwise false.
     */
    bool isSteiner() const;

    /**
     * Getter for the dimension of this point
     *
     * @return the dimension of this point
     */
    unsigned int dim() const;

    /**
     * Assignment operator. Copies all attributes.
     *
     * @param other   The point to assign.
     *
     * @return        A copy of this point.
     */
    Point& operator=(Point const& other);
    
    /**
     * Adds one point to this point (vector additon)
     *
     * @param other   The point to add.
     *
     * @return        This point plus other.
     */
    Point& operator+=(Point const& other);
    
    /**
     * Subtracts one point from this point (vector subtraction)
     *
     * @param other   The point to subtract by.
     *
     * @return        This point minus other.
     */
    Point& operator-=(Point const& other);
    
    /**
     * Multiplies this point by a scalar
     *
     * @param a    The scalar to multiply by.
     * 
     * @return     This point, multiplied by a
     */
    Point& operator*=(double const& a);
    
    /**
     * Compares this point to another point.
     * Returns true if two points are within
     * COMP_ERR distance of eachother.
     *
     * @param other  The point to compare against.
     *
     * @return       True if this point is equal
     *               to other
     */
    bool operator==(Point const& other);
    
    /**
     * Compares this point to another point.
     * Returns false if two points are within
     * COMP_ERR distance of eachother.
     *
     * @param other  The point to compare against.
     *
     * @return       True if this point is different
     *               from other
     */
    bool operator!=(Point const& other);
    
    /**
     * Allows for accessing the x, y and z
     * coordinates with the [] operator.
     *
     * @param i   The index, must be less than the size of the coordinate list
     *
     * @return    The value of the (i+1)th coordinate.
     */
    double const& operator[](const int i) const;
    
    /**
     * Allows for setting coordinates with the [] operator.
     *
     * @param i   The index, must less than the size of the coordinate list.
     *
     * @return    A reference to the (i+1)th coordinate.
     */
    double& operator[](const int i);
      
  protected:
    unsigned int _dim;
    double* _list;
    bool _isSteiner;
  private:
    /** the list of coordinates. */
    //std::vector<double> _list;
    /** isSteiner flag. True if this point is a Steiner point. */
    //bool _isSteiner;
  };

  /**
   * Returns the distance between point a and b
   *
   * @param a   The first point
   * @param b   The second point
   *
   * @return    The distance between a and b
   */
  double length(const Point &a, const Point &b);

  /**
   * Returns the distance between point a and the origin
   * (= the length of the vector a).
   * 
   * @param a   The point
   *
   * @return    The length of a / Distance from a to (0,0,0).
   */
  double length(const Point &a);

  /**
   * Returns the dot product of the points
   * a and b (vector dot product).
   *
   * @param a   The first point
   * @param b   The second point
   *
   * @return    The dot product of a and b
   */
  double dot(const Point &a, const Point &b);

  /**
   * Returns the cross product of the points
   * a and b (vector cross product).
   *
   * @param a   The first point
   * @param b   The second point
   *
   * @return    The cross product of a and b
   */
  Point cross(const Point &a, const Point &b);

  /**
   * Returns the normalised vector along the vector a.
   *
   * @param a    The vector to be normalised.
   *
   * @return     A vector with length == 1.
   */
  Point normalise(const Point &a);

  /**
   * Rotates a point t radians around a given axis.
   *
   * @param p   The point to be rotated.
   * @param v   A point on the axis of rotation.
   * @param u   A point on the axis of rotation different from v.
   * @param t   The number of radians to rotate.
   *
   * @return    p rotated t radians around the line segment vu.
   *
  Point rotate(Point p, Point v, Point u, double t);*/

  /**
   * @struct Edge
   *
   * Represents an edge in a graph.
   */
  struct Edge {
    /**
     * Constructor.
     *
     * @param i0   Index of first end point.
     * @param i1   Index of second end point.
     */
    Edge(int i0, int i1);
    
    /**
     * Constructor.
     *
     * @param i0   Index of first end point.
     * @param i1   Index of second end point.
     * @param len  Length of the edge. Needed for
     *             MSTKruskal and other.
     */
    Edge(int i0, int i1, double len);

    /**
     * Return a hashmap-key for this edge.
     */
    unsigned long key();

    /**
     * Generate a hashmap-key for this edge.
     */
    static unsigned long key(int i0, int i1);
    
    /** Index of first end point. */
    int i0;
    /** Index of second end point. */
    int i1;
    /** Length of this edge. */
    double length;
  };

  /**
   * Calculates arctan(abc)
   *
   * @param a    End of one line
   * @param b    Center point of angle
   * @param c    End of other line
   *
   * @return arctan(abc)
   */
  double angle(Utils::Point const &a, Utils::Point const &b, Utils::Point const &c);
};

/**
 * Adds two points.
 *
 * @param a   The first point.
 * @param b   The second point.
 * 
 * @return    a + b (vector addition).
 */ 
Utils::Point operator+(Utils::Point const& a, Utils::Point const& b);

/**
 * Subtracts two points.
 *
 * @param a   The first point.
 * @param b   The second point.
 * 
 * @return    a - b (vector subtraction).
 */
Utils::Point operator-(Utils::Point const& a, Utils::Point const& b);

/**
 * Right multiplies a point by a scalar.
 *
 * @param a   The scalar.
 * @param b   The point.
 * 
 * @return    a * b (scalar-vector multiplication).
 */
Utils::Point operator*(double const& a, Utils::Point const& b);

/**
 * Left multiplies a point by a scalar.
 *
 * @param a   The scalar.
 * @param b   The point.
 * 
 * @return    b * a (scalar-vector multiplication).
 */
Utils::Point operator*(Utils::Point const& b, double const& a);

/**
 * Prints the point to os in nice formatting.
 *
 * @param os     The out-stream
 * @param point  The point to print.
 *
 * @return       A reference to the stream.
 */
std::ostream& operator<<(std::ostream& os, Utils::Point const& point);

#endif // POINT_H
