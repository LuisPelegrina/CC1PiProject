//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#include "G4Particle.h"

bool G4Particle::is_contained() {
    bool is_contained = true;

    if(this->X0.X() > 200) is_contained = false;
    if(this->X0.X() < -200) is_contained = false;

    if(this->X0.Y() > 200) is_contained = false;
    if(this->X0.Y() < -200) is_contained = false;

    if(this->X0.Z() > 500) is_contained = false;
    if(this->X0.Z() < 0) is_contained = false;

    if(this->Xf.X() > 200) is_contained = false;
    if(this->Xf.X() < -200) is_contained = false;

    if(this->Xf.Y() > 200) is_contained = false;
    if(this->Xf.Y() < -200) is_contained = false;

    if(this->Xf.Z() > 500) is_contained = false;
    if(this->Xf.Z() < 0) is_contained = false;

    return is_contained;
}

bool G4Particle::is_contained_in_fv(double distance_to_wall) {
    bool is_contained = true;

    if(this->X0.X() > 200 - distance_to_wall) is_contained = false;
    if(this->X0.X() < -200 + distance_to_wall) is_contained = false;

    if(this->X0.Y() > 200 - distance_to_wall) is_contained = false;
    if(this->X0.Y() < -200 + distance_to_wall) is_contained = false;

    if(this->X0.Z() > 500 - distance_to_wall) is_contained = false;
    if(this->X0.Z() < 0 + distance_to_wall) is_contained = false;

    if(this->Xf.X() > 200 - distance_to_wall) is_contained = false;
    if(this->Xf.X() < -200 + distance_to_wall) is_contained = false;

    if(this->Xf.Y() > 200 - distance_to_wall) is_contained = false;
    if(this->Xf.Y() < -200 + distance_to_wall) is_contained = false;

    if(this->Xf.Z() > 500 - distance_to_wall) is_contained = false;
    if(this->Xf.Z() < 0 + distance_to_wall) is_contained = false;

    if ((this->X0.X() < 0 + distance_to_wall) && (this->X0.X() > 0 - distance_to_wall)) is_contained = false;
    if ((this->Xf.X() < 0 + distance_to_wall) && (this->Xf.X() > 0 - distance_to_wall)) is_contained = false;

    return is_contained;
}

void G4Particle::set_g4_particle(int ID, int mother, int PDG, double mass, double E0, TVector3 X0, TVector3 P0, double Ef, TVector3 Xf, TVector3 Pf, double t0, double tf, double TL, string start_process, string end_process ) {
  this->ID = ID;
  this->mother = mother;
  this->PDG = PDG;
  this->mass = mass;
  this->E0 = E0;
  this->Ef = Ef;
  this->X0 = X0;
  this->Xf = Xf;
  this->t0 = t0;
  this->tf = tf;
  this->TL = TL;
  this->P0 = P0;
  this->Pf = Pf;
  this->end_process = end_process;
  this->start_process = start_process;
}