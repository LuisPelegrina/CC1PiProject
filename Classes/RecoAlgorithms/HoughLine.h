//
// Created by Luis Pelegrina Gutiérrez on 10/10/24.
//
#include "SPoint.h"
#include "LineEquation.h"

#ifndef CC1PIPROJECT_HOUGHLINE_H
#define CC1PIPROJECT_HOUGHLINE_H


class HoughLine {
private:
  LineEquation fEquation;
  float fScore;
  int fNHits;
public:
  HoughLine(LineEquation line = LineEquation(0, 0), float score = -1, int nHits = 0);

  void set_line_equation(LineEquation eq){ fEquation = eq;}
  void set_score(float sc){ fScore = sc;}
  void set_nhits(float nhits){ fNHits = nhits;}

  float get_score(){return fScore;}
  float get_nhits(){return fNHits;}
  LineEquation get_line_equation(){return fEquation;}
};



#endif //CC1PIPROJECT_HOUGHLINE_H
