//
// Created by Luis Pelegrina Gutiérrez on 10/10/24.
//

#include "SPoint.h"

double SPoint::get_distance_to_point_w_erry(SPoint p) {
  // distance in X
  double dx = std::pow(p.get_x() - this->get_x(), 2);

  // distance between centers
  double d0 = std::sqrt(dx + std::pow(p.get_y() - this->get_y(), 2));

  // dist of hit 1 to width of hit 2
  double dp1 = std::sqrt(dx + std::pow(p.get_y() - this->get_y() + this->get_err_y(), 2));
  double dm1 = std::sqrt(dx + std::pow(p.get_y() - this->get_y() - this->get_err_y(), 2));
  double d1 = std::min(dp1, dm1);

  // dist of hit 2  to width of hit 1
  double dp2 = std::sqrt(dx + std::pow(p.get_y() - this->get_y() + p.get_err_y(), 2));
  double dm2 = std::sqrt(dx + std::pow(p.get_y() - this->get_y() - p.get_err_y(), 2));
  double d2 = std::min(dp2, dm2);

  // return min
  return std::min(d0, (d1 + d2) / 2.);
}

double SPoint::get_distance_to_point_center(SPoint p) {
  return sqrt(pow(fX - p.get_x(),2) + pow(fY - p.get_y(),2));
}