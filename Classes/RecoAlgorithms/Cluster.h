//
// Created by Luis Pelegrina Gutiérrez on 9/10/24.
//
#include "../Hit.h"
#include "TMath.h"
#include "reco_config.h"
#include "SPoint.h"
#include "LineEquation.h"
#include "TVector2.h"

#include <map>


#ifndef CC1PIPROJECT_CLUSTER_H
#define CC1PIPROJECT_CLUSTER_H

class Cluster {
public:
  Cluster(int ID = -1, int TPC_ID = -1, std::vector<Hit> input_hit_vec = {}, LineEquation line_eq = LineEquation(-1,-1,-1), SPoint pca_start_point = SPoint(0,0,0,0), SPoint pca_end_point = SPoint(0,0,0,0)):
    ID(ID),
    TPC_ID(TPC_ID),
    hit_vec(input_hit_vec),
    PCA_line_equation(line_eq),
    pca_start_point(pca_start_point),
    pca_end_point(pca_end_point)
  {}

  int ID;
  int TPC_ID;
  std::vector<Hit> hit_vec;
  SPoint pca_start_point;
  SPoint pca_end_point;
  LineEquation PCA_line_equation;

  double get_min_distance_to(int x, double y, int TPC_ID);
  double get_min_distance_to_cluster(Cluster cluster);
  double get_min_distance_to_cluster_w_err(Cluster cluster);
  double get_hit_mean_ocuppation();
  double get_mean_hit_width();
  LineEquation calculate_PCA_line_equation();
  std::map<int, double> calculate_connectedness_map();
};


#endif //CC1PIPROJECT_CLUSTER_H
