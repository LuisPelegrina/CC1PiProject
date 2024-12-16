//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_GENERATORPARTICLE_H
#define MUPIPROJECT_GENERATORPARTICLE_H

#include "Particle.h"
#include <iostream>
using namespace std;

class GeneratorParticle : public Particle {
public:

  GeneratorParticle()
    : Particle()
    , status_code(-1)
    {}

    int status_code;

    void print(bool print_vector) override {
        cout << "ID: " << this->ID << " SC: " << this->status_code  << " Mother: " << this->mother
             << " PDG: " << this->PDG  << " Mass: " << this->mass << " E0: "<< this->E0 << endl;
        if(print_vector) {
            cout << "   X0: (" <<  this->X0.X() << ", " <<  this->X0.Y() << ", "<<  this->X0.Z() << ")" << endl;
            cout << "   P0: (" <<  this->P0.X() << ", " <<  this->P0.Y() << ", "<<  this->P0.Z() << ")" << endl;
        }
    };

    bool is_primary() override {
        bool is_primary = (this->status_code == 1);
        if(this->PDG == 22 && this->P0.Mag() < 0.01) is_primary = false;
        if(this->PDG > 1e8) is_primary = false;

        return is_primary;
    };

    void set_gen_particle(int ID, int mother, int PDG, double mass, double E0, TVector3 X0, TVector3 P0, int status_code) {
      this->ID = ID;
      this->mother = mother;
      this->PDG = PDG;
      this->mass = mass;
      this->E0 = E0;
      this->X0 = X0;
      this->P0 = P0;
      this->status_code = status_code;
    };
};

#endif //MUPIPROJECT_GENERATORPARTICLE_H
