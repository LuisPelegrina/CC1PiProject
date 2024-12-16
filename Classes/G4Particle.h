//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//
#ifndef MUPIPROJECT_G4PARTICLE_H
#define MUPIPROJECT_G4PARTICLE_H

#include "TVector3.h"
#include "Particle.h"
class G4Particle : public Particle {
public:
  double Ef;
  double TL;
  double t0;
  double tf;
  TVector3 Xf;
  TVector3 Pf;
  string end_process;
  string start_process;

  G4Particle()
    : Particle()
    , Xf(TVector3(-1, -1, -1))
    {}

    void print(bool print_vector) override {
        cout << "ID: " << this->ID  << " Mother: " << this->mother
             << " PDG: " << this->PDG  << " Mass: " << this->mass << " E0: "<< this->E0  << endl;
        if(print_vector) {
            cout << "   X0: (" <<  this->X0.X() << ", " <<  this->X0.Y() << ", "<<  this->X0.Z() << ")" << endl;
            cout << "   XL: (" <<  this->Xf.X() << ", " <<  this->Xf.Y() << ", "<<  this->Xf.Z() << ")" << endl;
            cout << "End Process: " << this->end_process << endl;
        }
    }

    void set_g4_particle(int ID, int mother, int PDG, double mass, double E0, TVector3 X0, TVector3 P0, double Ef, TVector3 Xf, TVector3 Pf, double t0, double tf, double TL, string start_process, string end_process );

    bool is_primary() override {
        bool is_primary = (this->mother == 0);
        if(this->PDG == 22 && this->P0.Mag() < 0.01) is_primary = false;
        if(this->PDG > 1e8) is_primary = false;

        return is_primary;
    }

    bool is_contained();
    bool is_contained_in_fv(double distance_to_wall);
};



#endif //MUPIPROJECT_G4PARTICLE_H
