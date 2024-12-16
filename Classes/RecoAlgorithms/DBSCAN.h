//
// Created by Luis Pelegrina Gutiérrez on 30/10/24.
//

#ifndef CC1PIPROJECT_DBSCAN_H
#define CC1PIPROJECT_DBSCAN_H
#include "SPoint.h"

class DBSCAN {
public:
  DBSCAN(double epsilon, int min_pts);
  void fit( std::vector<SPoint>& points);
  void setDistanceFunction(double (*distFunc)(SPoint&,  SPoint&));

  std::vector<int>& getClusterAssignment();
  // --------------- Distance functions for DBSCAN
  double DBSCANHitEuclidianDistance( SPoint& p1,  SPoint& p2);
  double DBSCANHitWidthDistance( SPoint& p1,  SPoint& p2);

private:
  double epsilon;
  int min_pts;
  std::vector<int> clusterAssignment;
  std::vector<int> regionQuery( std::vector<SPoint>& points,  SPoint& p);
  double (*distanceFunction)( SPoint&,  SPoint&);
  void expandCluster( std::vector<SPoint>& points, int point_Idx, int cluster_Idx);


};


#endif //CC1PIPROJECT_DBSCAN_H
