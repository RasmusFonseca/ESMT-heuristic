#include <vector>
#include <iostream>
#include <cmath>

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "steiner/utils/point.hpp"

Utils::Point::Point() {
  this->_dim  = 3;
  this->_list = new double[3];
  for (int i = 0; i < 3; i++)
    this->_list[i] = 0.0;
  this->_isSteiner = false;
}

Utils::Point::Point(int dimension, double value) {
  this->_dim  = dimension;
  this->_list = new double[dimension];
  for (int i = 0; i < dimension; i++)
    this->_list[i] = value;
  this->_isSteiner = false;
}


Utils::Point::Point(double x, double y, double z, bool isSteiner) {
  this->_dim  = 3;
  this->_list = new double[3];
  this->_list[0] = x;
  this->_list[1] = y;
  this->_list[2] = z;
  this->_isSteiner = isSteiner;
}

Utils::Point::Point(std::vector<double> &list, bool isSteiner) {
  this->_dim  = list.size();
  this->_list = new double[list.size()];
  for(unsigned int i = 0; i < list.size(); i++)
    this->_list[i] = list[i];
  this->_isSteiner = isSteiner;
}

Utils::Point::Point(Utils::Point const& other) {
  this->_dim  = other._dim;
  this->_list = new double[other._dim];
  for(unsigned int i = 0; i < other._dim; i++)
    this->_list[i] = other._list[i];
  this->_isSteiner = other._isSteiner;
}

Utils::Point::~Point() {
  delete this->_list;
}

void Utils::Point::setSteiner(bool isSteiner) {
  this->_isSteiner = isSteiner;
}

bool Utils::Point::isSteiner() const {
  return this->_isSteiner;
}

unsigned int Utils::Point::dim() const {
  return this->_dim;
}

double *Utils::Point::getList() const {
  return this->_list;
}

Utils::Point& Utils::Point::operator=(Utils::Point const& other) {
  if(&other == this)
    return *this;
  if(other._dim != this->_dim) {
    delete this->_list;
    this->_list = new double[other._dim];
    this->_dim  = other._dim;
  }
  for(unsigned int i = 0; i < this->_dim; i++)
    this->_list[i] = other._list[i];
  this->_isSteiner = other._isSteiner;

  return *this;
}

bool Utils::Point::operator!=(Utils::Point const& other) {
  return !(*this==other);
}

Utils::Point& Utils::Point::operator+=(Utils::Point const& other) {
  for (unsigned int i = 0; i < this->_dim; i++)
    this->_list[i]+=other._list[i];
  return *this;
}

Utils::Point& Utils::Point::operator-=(Utils::Point const& other) {
  for (unsigned int i = 0; i < this->_dim; i++)
    this->_list[i]-=other._list[i];
  return *this;
}

Utils::Point& Utils::Point::operator*=(double const& a) {
  for (unsigned int i = 0; i < this->_dim; i++)
    this->_list[i] *= a;
  return *this;
}

bool Utils::Point::operator==(Utils::Point const& other) {
  if (this->_dim != other._dim)
    return false;
  for (unsigned int i = 0; i < this->_dim; i++)
    if (std::abs(this->_list[i]-other._list[i]) > COMP_ERR)
	return false;
  return true;
}

Utils::Point operator+(Utils::Point const& a, Utils::Point const& b) {
  Utils::Point result(a.dim());
  for (unsigned int i = 0; i < a.dim(); i++)
    result[i] = a[i]+b[i];
  return result;
}

Utils::Point operator-(Utils::Point const& a, Utils::Point const& b) {
  Utils::Point result = Utils::Point(a.dim());
  for (unsigned int i = 0; i < a.dim(); i++)
    result[i] = a[i]-b[i];
  return result;
}

Utils::Point operator*(double const& a, Utils::Point const& b) {
  Utils::Point result = Utils::Point(b.dim());
  for (unsigned int i = 0; i < b.dim(); i++)
    result[i] = a*b[i];
  return result;
}

Utils::Point operator*(Utils::Point const& b, double const& a) {
  Utils::Point result = Utils::Point(b.dim());
  for (unsigned int i = 0; i < b.dim(); i++)
    result[i] = a*b[i];
  return result;
}

double const& Utils::Point::operator[](const int i) const {
  assert(i >= 0 && (unsigned int)i < this->_dim);
  return this->_list[i];
}

double& Utils::Point::operator[](const int i) {
  assert(i >= 0 && (unsigned int)i < this->_dim);
  return this->_list[i];
}

std::ostream& operator<<(std::ostream& os, Utils::Point const& point) {
  os << "(";
  for (unsigned int i = 0; i < point.dim(); i++) {
    os << point[i];
    if (i != point.dim() - 1)
      os << ", ";
  }
  os << ")";
  return os;
}

double Utils::length(const Utils::Point &a, const Utils::Point &b) {
  double s = 0, diff = 0;
  for (unsigned int i = 0; i < a.dim(); i++) {
    diff = a[i]-b[i];
    s += diff*diff;
  }
  return sqrt(s);
}

double Utils::length(const Utils::Point &a) {
  double s = 0;
  for (unsigned int i = 0; i < a.dim(); i++)
    s += a[i]*a[i];
  return sqrt(s);
}

double Utils::dot(const Utils::Point &a, const Utils::Point &b) {
  if (a.dim() != b.dim())
    return 0;
  double result = 0;
  for (unsigned int i = 0; i < a.dim(); i++)
    result += a[i]*b[i];
  return result;
}

Utils::Point Utils::cross(const Utils::Point &a, const Utils::Point &b) {
  assert(a.dim() == 3 && b.dim() == 3);
  return Utils::Point(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

Utils::Point Utils::normalise(const Point &a) {
  return (1/Utils::length(a, Utils::Point()))*a;
}

/* Implementation of angle */
double Utils::angle(Utils::Point const &a, Utils::Point const &b, Utils::Point const &c) {
  return Utils::dot(b-a, b-c) / (Utils::length(b, a) * Utils::length(b, c));
}


////////////////////////////////////////////////////
// Edge functions

/* Constructors */
Utils::Edge::Edge(int i0, int i1) : length(-1) {
  if(i0 > i1) {
    this->i0 = i1;
    this->i1 = i0;
  }
  else {
    this->i0 = i0;
    this->i1 = i1;
  }
}
Utils::Edge::Edge(int i0, int i1, double len) : length(len) {
  if(i0 > i1) {
    this->i0 = i1;
    this->i1 = i0;
  }
  else {
    this->i0 = i0;
    this->i1 = i1;
  }
}
    
/* Hashkey */
unsigned long Utils::Edge::key() {
  unsigned long res = this->i1;
  res = (res << 32) | this->i0;
  return res;
}
unsigned long Utils::Edge::key(int i0, int i1) {
  unsigned long a, b;
  if(i0 < i1) {
    a = i0;
    b = i1;
  }
  else {
    a = i1;
    b = i0;
  }
  return (b << 32) | a;
}
