//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//
#include "TMath.h"
#include "reco_config.h"

#ifndef MUPIPROJECT_ReducedVertexWire_H
#define MUPIPROJECT_ReducedVertexWire_H


class ReducedVertexWire {
public:
  ReducedVertexWire()
    :ID(-1)
    ,Wire_ID(-1)
    ,drift_t(-1)
    ,Channel_ID(-1)
    ,TPC_ID(-1)
    ,Plane_ID(-1)
    {}

    int ID;
    int TPC_ID;
    int Plane_ID;
    int Channel_ID;
    int Wire_ID;
    double drift_t;

    double get_distance_to(ReducedVertexWire pnt);
};

#endif //MUPIPROJECT_ReducedVertexWire_H