//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_PARTICLE_H
#define MUPIPROJECT_PARTICLE_H

#include "TVector3.h"
#include <iostream>
using namespace std;

class Particle {
public:

  Particle()
    : ID(-1)
    , mother(-1)
    , PDG(-1)
    , mass(-1)
    , E0(-1)
    , X0(TVector3(-1, -1, -1))
    , P0(TVector3(-1, -1, -1))
  { }

    int ID;
    int mother;
    int PDG;
    double mass;
    double E0;
    TVector3 X0;
    TVector3 P0;

    virtual bool is_primary() {return true;};

    virtual void print(bool print_vector) { cout << print_vector << endl;};
};


#endif //MUPIPROJECT_PARTICLE_H
