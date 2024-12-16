//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_HIT_H
#define MUPIPROJECT_HIT_H

#include "RecoAlgorithms/SPoint.h"

class Hit {
public:
  Hit()
    :ID(-1)
    ,integral(-1)
    ,integral_sigma(-1)
    ,drift_t(-1)
    ,drift_t_sigma(-1)
    ,TPC_ID(-1)
    ,Plane_ID(-1)
    ,Wire_ID(-1)
    ,channel_ID(-1)
    ,associated_pfp_ID(-1)
    ,plane_point(SPoint(-1.,-1.,0.,-1.))
  {}

  int ID;
  double integral;
  double integral_sigma;
  double drift_t;
  double drift_t_sigma;
  int TPC_ID;
  int Plane_ID;
  int Wire_ID;
  int channel_ID;
  int associated_pfp_ID;
  SPoint plane_point;

  void set_hit(double integral, double sigma_integral, double drift_t, double drift_t_sigma, int TPC_ID, int Plane_ID, int Wire_ID, int Channel_ID, int associated_pfp_ID);
};




#endif //MUPIPROJECT_HIT_H
