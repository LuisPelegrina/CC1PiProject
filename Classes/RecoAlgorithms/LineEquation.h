//
// Created by Luis Pelegrina Gutiérrez on 10/10/24.
//

#ifndef CC1PIPROJECT_LINEEQUATION_H
#define CC1PIPROJECT_LINEEQUATION_H

#include "SPoint.h"
#include "TMath.h"

class LineEquation {
private:
  float fM;
  float fN;
  float fGoodness;

public:
  LineEquation(float slope=0, float intercept=0, float goodness=0):
    fM(slope),
    fN(intercept),
    fGoodness(goodness)
  {}
  //void SetGoodness(float goodness);

  float get_slope() {return fM;}
  float get_intercept() {return fN;}
  float get_goodness() {return fGoodness;}

  float get_distance(SPoint p);
  SPoint get_line_closest_point(double a, double b, double c, SPoint p);
  SPoint get_line_closest_point(SPoint P);

  float evaluate(SPoint p);
  float evaluate_x(double x);

};



#endif //CC1PIPROJECT_LINEEQUATION_H
