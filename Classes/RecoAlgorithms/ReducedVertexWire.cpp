//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#include "ReducedVertexWire.h"

double ReducedVertexWire::get_distance_to(ReducedVertexWire pnt) {

  double delta_wire = pow(this->drift_t/4 - pnt.drift_t/4, 2);
  double delta_t = pow(this->Wire_ID - pnt.Wire_ID, 2);

  return sqrt(delta_t + delta_wire);
}
