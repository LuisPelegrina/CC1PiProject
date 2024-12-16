//
// Created by Luis Pelegrina Gutiérrez on 30/10/24.
//

#include "DBSCAN.h"
DBSCAN::DBSCAN(double epsilon, int min_pts)
  :epsilon(epsilon),
   min_pts(min_pts)
{}

void DBSCAN::setDistanceFunction(double (*distFunc)( SPoint&,  SPoint&)) {
  distanceFunction = distFunc;
}


std::vector<int>& DBSCAN::getClusterAssignment()  {
  return clusterAssignment;
}

void DBSCAN::fit( std::vector<SPoint>& points) {
  clusterAssignment.assign(points.size(), 0);  // 0 indicates unassigned

  int clusterIdx = 0;
  for (size_t i = 0; i < points.size(); ++i) {
    if (clusterAssignment[i] == 0) {  // Not assigned to any cluster
      std::vector<int> neighbors = regionQuery(points, points[i]);
      if ((int)neighbors.size() < min_pts) {
        clusterAssignment[i] = -1;  // Mark as noise
      } else {
        ++clusterIdx;
        expandCluster(points, i, clusterIdx);
      }
    }
  }
}

std::vector<int> DBSCAN::regionQuery( std::vector<SPoint>& points,  SPoint& p) {
  std::vector<int> neighbors;
  for (size_t i = 0; i < points.size(); ++i) {
    if (distanceFunction(p, points[i]) <= epsilon) {
      neighbors.push_back(i);
    }
  }
  return neighbors;
}


void DBSCAN::expandCluster( std::vector<SPoint>& points, int pointIdx, int clusterIdx) {
  //Gets the neighbours
  std::vector<int> seeds = regionQuery(points, points[pointIdx]);

  //Assings the point that started the cluster the corresponding ID
  clusterAssignment[pointIdx] = clusterIdx;

  //Do a loop through all the neighbours
  for (size_t i = 0; i < seeds.size(); ++i) {
    //Get the point ID of the seed
    int currentIdx = seeds[i];

    //If seed is not asigned...
    if (clusterAssignment[currentIdx] == 0) {
      //Assing it to the cluster we are expanding
      clusterAssignment[currentIdx] = clusterIdx;

      //If the neighbours are enought expand the current cluster
      std::vector<int> currentSHitNeighbors = regionQuery(points, points[currentIdx]);
      if ((int)currentSHitNeighbors.size() >= min_pts) {
        expandCluster(points, currentIdx, clusterIdx);
      }
    }

    // If the point is considered noise
    if (clusterAssignment[currentIdx] == -1) {
      clusterAssignment[currentIdx] = clusterIdx;
    }
  }
}

// --------------- Distance functions definitions for DBSCAN

double DBSCANHitEuclidianDistance( SPoint& p1,  SPoint& p2) {
  return std::sqrt((p1.get_x() - p2.get_x()) * (p1.get_x() - p2.get_x()) + (p1.get_y() - p2.get_y()) * (p1.get_y() - p2.get_y()));
}

double DBSCANHitWidthDistance(SPoint& p1,  SPoint& p2) {
  //std::cout<<" In width distance\n";
  double dX = std::pow(p1.get_x() - p2.get_x(), 2);
  double d0 = std::sqrt( dX + std::pow(p1.get_y() - p2.get_y(), 2));

  double dY1 =  std::min( std::abs(p1.get_y() + p1.get_err_y() - p2.get_y()),  std::abs(p1.get_y() - p1.get_err_y() - p2.get_y()));
  double dY2 =  std::min( std::abs( p1.get_y() - p2.get_y() + p2.get_err_y()),  std::abs(p1.get_y() - p2.get_y() - p2.get_err_y()) );

  double d1 = std::sqrt( dX + std::pow(dY1, 2));
  double d2 = std::sqrt( dX + std::pow(dY2, 2));

  //return std::min(d0, std::min( d1, d2 ));
  return std::min(d0, 0.5*(d1+d2));
}