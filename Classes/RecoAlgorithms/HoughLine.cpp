//
// Created by Luis Pelegrina Gutiérrez on 10/10/24.
//

#include "HoughLine.h"

HoughLine::HoughLine(LineEquation line = LineEquation(0, 0), float score = -1, int nHits = 0):
  fEquation(line),
  fScore(score),
  fNHits(nHits)
{}

LineEquation ComputeLineEquationFromHough(double th, SPoint p) {
  double m = std::tan(th);
  double n = p.get_y() - p.get_x() * m;
  return LineEquation(m, n);
}
