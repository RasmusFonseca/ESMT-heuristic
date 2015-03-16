#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/heuristics/steiner_finder.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/fermat.hpp"

typedef Utils::Point Point;
typedef Utils::Edge Edge;

SteinerFinder::SteinerFinder(SubgraphHeuristic* sh) {
  this->sh = sh;
}

SteinerFinder::~SteinerFinder() {}

SteinerTree* SteinerFinder::findSteinerPoints(Graph &subgraph) {
  std::vector<Point> *points = subgraph.getPointsPtr();
  this->mst_length = subgraph.getLength();
  this->dim        = subgraph.dimension();
  
  assert(points->size() > 2);

  if(points->size() == 3)
    return Utils::getFermatSMT(subgraph);
  
  this->best_ratio  = 2.0;
  this->best_tree   = NULL;
  this->best_length = this->mst_length;
  unsigned int i;
  for(i = 0; i < points->size(); i++) {
    std::vector<int> part;
    part.push_back(i);
    this->simplexPartitionRec(*points, part, i);
  }
  return best_tree;
}

void SteinerFinder::simplexPartitionRec(std::vector<Point> &points, std::vector<int> &curPart, unsigned int index) {
  unsigned int i;
  SteinerTree *st;
  if(points.size()-curPart.size() < 3)
    return;
  for(i = index+1; i < points.size(); i++) {
    curPart.push_back(i);
    st = this->merge(points, curPart);
    if(st != NULL && st->getSteinerRatio() < this->best_ratio) {
      if(this->best_tree)
	delete this->best_tree;
      this->best_tree   = st;
      this->best_length = st->getSMTLength();
      this->best_ratio  = st->getSteinerRatio();
    }
    else
      delete st;
    this->simplexPartitionRec(points, curPart, i);
    curPart.pop_back();
  }
}

/* partOne may be assumed to be sorted! */
SteinerTree *SteinerFinder::merge(std::vector<Point> &points,
				  std::vector<int> &partOne) {
  Graph p[2];
  std::vector<int> map[2];
  unsigned int i, j;
  SteinerTree *st[2];
  for(i = 0, j = 0; i < points.size(); i++) {
    if(j < partOne.size() && partOne[j] == (int)i) {
      p[0].getPointsPtr()->push_back(points[i]);
      map[0].push_back(i);
      j++;
    }
    else {
      p[1].getPointsPtr()->push_back(points[i]);
      map[1].push_back(i);
    }
  }
  /*
  std::cout << "P1: ";
  for(i = 0; i < map[0].size(); i++)
    std::cout << map[0][i] << " ";
  std::cout << std::endl;

  std::cout << "P2: ";
  for(i = 0; i < map[1].size(); i++)
    std::cout << map[1][i] << " ";
  std::cout << std::endl;
  */
  p[0].setMSTLength(this->mst_length);
  p[1].setMSTLength(this->mst_length);

  Point c = centroid(p[1]);
  //std::cout << "C: " << c << std::endl;
  p[0].getPointsPtr()->push_back(c);
  map[0].push_back(-1);
  if(map[0].size() > 3)
    st[0] = this->sh->findSteinerPoints(p[0]);
  else
    st[0] = Utils::getFermatSMT(p[0], true);
  
  if(st[0]->getSMTLength() > this->best_length) {
    delete st[0];
    return NULL;
  }
  
  // We know, that val(c) == 1
  int c_pos = map[0].size()-1;
  unsigned int s = 0, e = 0;
  for (std::vector<Edge>::iterator it = st[0]->getEdgesPtr()->begin();
       it != st[0]->getEdgesPtr()->end(); ++it) {
    if(it->i0 == c_pos || it->i1 == c_pos) {
      s = it->i0 == c_pos ? it->i1 : it->i0;
      break;
    }
    ++e;
  }
  //std::cout << "s= " << s << std::endl;

  p[1].getPointsPtr()->push_back((*st[0]->getPointsPtr())[s]);
  (*p[1].getPointsPtr())[p[1].getPointsPtr()->size()-1].setSteiner(false);
  map[1].push_back(-2);
  if(map[1].size() > 3)
    st[1] = this->sh->findSteinerPoints(p[1]);
  else
    st[1] = Utils::getFermatSMT(p[1], true);
  
  if(st[1]->getSMTLength() > this->best_length) {
    delete st[0];
    delete st[1];
    return NULL;
  }

  //std::cout << std::endl << "___________________________\n" << *st[1] << std::endl << "___________________________\n";

  
  SteinerTree *result = new SteinerTree();
  for(std::vector<Point>::iterator it = points.begin();
      it != points.end(); ++it) {
    result->getPointsPtr()->push_back(*it);
  }
  
  int offsetp1 = p[0].getPointsPtr()->size();
  int offsetp1sp = points.size()-offsetp1;
  
  for(i = offsetp1; i < st[0]->getPointsPtr()->size(); i++)
      result->getPointsPtr()->push_back((*st[0]->getPointsPtr())[i]);

  int offsetp2 = p[1].getPointsPtr()->size();
  int offsetp2sp = result->getPointsPtr()->size()-offsetp2;
  for(i = offsetp2; i < st[1]->getPointsPtr()->size(); i++)
      result->getPointsPtr()->push_back((*st[1]->getPointsPtr())[i]);
  
  for(i = 0; i < st[0]->getEdgesPtr()->size(); i++) {
    if(i == e)
      continue;
    Edge edge = (*st[0]->getEdgesPtr())[i];
    //std::cout << "l1: " << edge.i0 << " - " << edge.i1 << std::endl;
    edge.i0 = edge.i0 >= offsetp1 ? edge.i0+offsetp1sp : map[0][edge.i0];
    edge.i1 = edge.i1 >= offsetp1 ? edge.i1+offsetp1sp : map[0][edge.i1];
    result->getEdgesPtr()->push_back(edge);
  }
  
  for(i = 0; i < st[1]->getEdgesPtr()->size(); i++) {
    Edge edge = (*st[1]->getEdgesPtr())[i];
    //std::cout << "l2: " << edge.i0 << " - " << edge.i1 << std::endl;
    edge.i0 = edge.i0 >= offsetp2 ? edge.i0+offsetp2sp : map[1][edge.i0];
    edge.i1 = edge.i1 >= offsetp2 ? edge.i1+offsetp2sp : map[1][edge.i1];
    edge.i0 = edge.i0 == -2 ? s+offsetp1sp : edge.i0;
    edge.i1 = edge.i1 == -2 ? s+offsetp1sp : edge.i1;
    result->getEdgesPtr()->push_back(edge);
  }
  
  delete st[0];
  delete st[1];
  
  result->setMSTLength(this->mst_length);
  return result;
}

Point SteinerFinder::centroid(Graph& subgraph) {
  Point result(this->dim);
  int size = subgraph.getPointsPtr()->size();
  for (int i = 0; i < size; i++)
    for (int j = 0; j < this->dim; j++)
      result[j] += (*subgraph.getPointsPtr())[i][j] / size;

  return result;
}

