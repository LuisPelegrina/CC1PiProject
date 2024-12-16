//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//
#include "TMath.h"

#ifndef MUPIPROJECT_RECOVERTEXWIRE_H
#define MUPIPROJECT_RECOVERTEXWIRE_H


class RecoVertexWire {
public:
  RecoVertexWire()
    :Wire_ID(-1)
    ,drift_t(-1)
    ,Channel_ID(-1)
    {}

    double Channel_ID;
    double Wire_ID;
    double drift_t;

    double get_distance_to(RecoVertexWire pnt);
};

#endif //MUPIPROJECT_RECOVERTEXWIRE_H
