//
// Created by Luis Pelegrina Gutiérrez on 10/10/24.
//

#include "TMath.h"

#ifndef CC1PIPROJECT_SPOINT_H
#define CC1PIPROJECT_SPOINT_H


class SPoint {
private:
  float fX;
  float fY;
  float fdx;
  float fdy;

public:


  SPoint(int x, double y, double err_x, double err_y)
    : fX((float)x),
      fY((float)y),
      fdx((float)err_x),
      fdy((float)err_y)
  {}

  SPoint(double x, double y, double err_x, double err_y)
    : fX((float)x),
      fY((float)y),
      fdx((float)err_x),
      fdy((float)err_y)
  {}

  SPoint(float x, float y, float err_x, float err_y)
    : fX(x),
      fY(y),
      fdx(err_x),
      fdy(err_y)
  {}

  float get_x() const {return fX;}
  float get_y() const {return fY;}
  float get_err_x() const {return fdx;}
  float get_err_y() const {return fdy;}
  void set_x(float x) { fX = x;}
  void set_y(float y) { fY = y;}
  void set_err_x(float err_x) { fdx = err_x;}
  void set_err_y(float err_y) { fdy = err_y;}

  double get_distance_to_point_center(SPoint p);
  double get_distance_to_point_w_erry(SPoint p);
};



#endif //CC1PIPROJECT_SPOINT_H
