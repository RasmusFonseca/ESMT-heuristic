#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

#include <vector>

#include "steiner/utils/point.hpp"

namespace Utils {
  /**
   * This class implements a disjoint set data structure, as described
   * in T. H. Cormens et al. (Introduction to Algorithms).
   *
   */
  template <class T>
  class DisjointSet {
  public:
  
    /**
     * Constructor
     * Equivalent to T.H. Cormens et al. Make-Set
     * Initialises the set with one point.
     *
     * @param value   The value of this Disjoint set.
     */
    DisjointSet(T value);
  
    /**
     * Copy constructor
     *
     * @param set  The DisjointSet to copy.
     */
    DisjointSet(const DisjointSet<T> &set);
  
    /**
     * Destructor
     */
    ~DisjointSet();

    /**
     * The Find-Set procedure
     *
     * @return  The representant for this set.
     */
    DisjointSet<T> *findSet();

    /**
     * Returns the value of the given set
     *
     * @return The value for this set.
     */
    T getValue();

    /**
     * The Union procedure
     *
     * @param otherSet  The DisjointSet to union with.
     */
    void setUnion(DisjointSet<T> *otherSet);

    /**
     * Equallity operator.
     *
     * @param otherSet  The set to compare with.
     * @return          True if otherSet.value == this.value.
     */
    bool operator==(DisjointSet<T> const &otherSet);
    
    /**
     * In-equallity operator.
     *
     * @param otherSet  The set to compare with.
     * @return          True if otherSet.value != this.value.
     */
    bool operator!=(DisjointSet<T> const &otherSet);

  protected:
  private:
    /**
     * The linking procedure.
     *
     * @param otherSet  The set to link with.
     */
    void link(DisjointSet<T> *otherSet);

    /** The value of this set. */
    T value;
    /** The rank of this set. */
    int rank;
    /** The parent_set of this set. */
    DisjointSet *parent_set;
  };
}
#endif // DISJOINT_SET_H
