#include <vector>
#include <iostream>

#include "steiner/utils/disjoint_set.hpp"
#include "steiner/utils/point.hpp"

/*
 * Constructor
 * Equivalent to T.H. Cormens et al. Make-Set
 * Initialises the set with one point.
 */
template <class T>
Utils::DisjointSet<T>::DisjointSet(T value) {
  this->rank       = 0;
  this->value      = value;
  this->parent_set = this;
}

/*
 * Copy constructor
 */
template <class T>
Utils::DisjointSet<T>::DisjointSet(const DisjointSet<T> &set) {
  this->rank       = set.rank;
  this->value      = set.value;
  this->parent_set = this;
}

/*
 * Destructor
 */
template <class T>
Utils::DisjointSet<T>::~DisjointSet() {
}

/*
 * The Find-Set procedure
 */ 
template <class T>
Utils::DisjointSet<T> *Utils::DisjointSet<T>::findSet() {
  if(*this != *(this->parent_set))
    this->parent_set = this->parent_set->findSet();
  return this->parent_set;
}

/*
 * Returns the value of the given set
 */
template <class T>
T Utils::DisjointSet<T>::getValue() {
  return this->value;
}

/*
 * The Union procedure
 */
template <class T>
void Utils::DisjointSet<T>::setUnion(DisjointSet<T> *otherSet) {
  Utils::DisjointSet<T> *a = this->findSet();
  Utils::DisjointSet<T> *b = otherSet->findSet();
  a->link(b);
}

/*
 * The linking procedure
 */ 
template <class T>
void Utils::DisjointSet<T>::link(Utils::DisjointSet<T> *otherSet) {
  if(this->rank > otherSet->rank) {
    otherSet->parent_set = this;
  }
  else {
    this->parent_set = otherSet;
    if(this->rank == otherSet->rank) {
      otherSet->rank++;
    }
  }
}

/*
 * Comparison operators
 */ 
template <class T>
bool Utils::DisjointSet<T>::operator==(Utils::DisjointSet<T> const &otherSet) {
  return this->value == otherSet.value;
}
template <class T>
bool Utils::DisjointSet<T>::operator!=(Utils::DisjointSet<T> const &otherSet) {
  return !(this->value == otherSet.value);
}

// Needed for each template use (or linker will go bananas)
template class Utils::DisjointSet<Utils::Point*>;
template class Utils::DisjointSet<int>;
