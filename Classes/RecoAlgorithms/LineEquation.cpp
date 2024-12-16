//
// Created by Luis Pelegrina Gutiérrez on 10/10/24.
//

#include "LineEquation.h"



SPoint LineEquation::get_line_closest_point(double a, double b, double c, SPoint p){
  // line in form ax+by+c=0
  float x = (b * (b * p.get_x() - a * p.get_y()) - a * c) / (pow(a,2) + pow(b,2));
  float y = (a * (-b * p.get_x() + a * p.get_y()) - b * c) / (pow(a,2) + pow(b,2));

  return SPoint(x, y, 0, 0);
}

SPoint LineEquation::get_line_closest_point(SPoint P) {
  return get_line_closest_point(-fM, 1, -fN, P);
}


float LineEquation::get_distance(SPoint p) {
  SPoint p_proj = get_line_closest_point(-fM, 1, -fN, p);
  return sqrt(pow(p.get_x() - p_proj.get_x(),2) + pow(p.get_y() - p_proj.get_y(),2));
}

float LineEquation::evaluate(SPoint p) {
  return fM * p.get_x() + fN;
}


float LineEquation::evaluate_x(double x) {
  return fM * x + fN;
}