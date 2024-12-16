//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_VERTEX_H
#define MUPIPROJECT_VERTEX_H

#include "RecoVertexWire.h"
#include <iostream>
#include "TVector3.h"
using namespace std;

class Vertex {
    public:
      Vertex()
      : vertex_wire({RecoVertexWire(), RecoVertexWire(), RecoVertexWire()})
      , corrected_vertex_wire({RecoVertexWire(), RecoVertexWire(), RecoVertexWire()})
      , vertex_cm(TVector3(-1, -1, -1))
      {}
        std::vector<RecoVertexWire> vertex_wire;
        std::vector<RecoVertexWire> corrected_vertex_wire;
        TVector3 vertex_cm;
};


#endif //MUPIPROJECT_VERTEX_H
