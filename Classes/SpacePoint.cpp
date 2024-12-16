//
// Created by Luis Pelegrina Gutiérrez on 17/9/24.
//

#include "SpacePoint.h"

void SpacePoint::set_space_point(double x, double y, double z, double integral, double sigma_integral, double associated_pfp_ID) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->integral = integral;
  this->sigma_integral = sigma_integral;
  this->associated_pfp_ID = associated_pfp_ID;
}