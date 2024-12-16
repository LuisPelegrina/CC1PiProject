//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#include "RecoVertexWire.h"

double RecoVertexWire::get_distance_to(RecoVertexWire pnt) {

  double delta_wire = pow(this->drift_t/4 - pnt.drift_t/4, 2);
  double delta_t = pow(this->Wire_ID - pnt.Wire_ID, 2);

  return sqrt(delta_t + delta_wire);
}