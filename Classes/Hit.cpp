//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#include "Hit.h"


void Hit::set_hit(double integral, double sigma_integral, double drift_t, double drift_t_sigma, int TPC_ID, int Plane_ID, int Wire_ID, int Channel_ID, int associated_pfp_ID) {
  this->integral = integral;
  this->integral_sigma = sigma_integral;
  this->drift_t = drift_t;
  this->drift_t_sigma = drift_t_sigma;
  this->TPC_ID = TPC_ID;
  this->Wire_ID = Wire_ID;
  this->Plane_ID = Plane_ID;
  this->channel_ID = Channel_ID;
  this->associated_pfp_ID = associated_pfp_ID;
  SPoint p(Wire_ID, drift_t/4, 0., drift_t_sigma/4) ;
  this->plane_point = p;
}