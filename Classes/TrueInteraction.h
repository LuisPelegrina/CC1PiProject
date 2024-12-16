//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//
#ifndef MUPIPROJECT_TRUEINTERACTION_H
#define MUPIPROJECT_TRUEINTERACTION_H

#include <iostream>
using namespace std;

#include "GeneratorParticle.h"
#include "G4Particle.h"
#include "Structures.h"



class TrueInteraction {
public:

  TrueInteraction()
  : interaction_mode(-1)
  , interaction_type(-1)
  , E0_nu(-1)
  , PDG_nu(-1)
  , CCNC(-1)
  , target(-1)
  , HitNuc(-1)
  , HitQuark(-1)
  , W(-1)
  , X(-1)
  , Y(-1)
  , QSqr(-1)
  ,primary_vertex(TVector3(-1,-1,-1))
  ,gen_primary_particles({GeneratorParticle()})
  ,gen_particles({GeneratorParticle()})
  ,g4_particles({G4Particle()})
  {}

  int interaction_mode;
  int interaction_type;
  double E0_nu;
  int PDG_nu;
  int CCNC;
  int target;
  int HitNuc;
  int HitQuark;
  double W;
  double X;
  double Y;
  double QSqr;
  TVector3 primary_vertex;
  std::vector<GeneratorParticle> gen_primary_particles;
  std::vector<GeneratorParticle> gen_particles;
  std::vector<G4Particle> g4_particles;

    void print(bool b_print_gen_particles, bool b_print_g4_particles, bool print_vector);
    void print_gen_particles(std::vector<GeneratorParticle> particles_to_print, bool print_vector);
    void print_g4_particles(std::vector<G4Particle> particles_to_print, bool print_vector, int iter,int  mother_ID);

    std::vector<GeneratorParticle>  get_generator_primary_p(int target_PDG);
    std::vector<GeneratorParticle>  get_generator_primary_p_sing(int target_PDG);
    int get_generator_num_primaries();
    int get_generator_num_primary_p(int target_PDG);
    int get_generator_num_primary_p_high_P0(int target_PDG, double P0_threshold);
    int get_generator_num_primary_p_sing(int target_PDG);

    G4Particle get_g4_particle_by_id(int id);

    bool pass_containment_cut(string cut_type, double fv_distance);
    bool is_selected_final_state(string cut_type, double high_p0_threshold);
    bool is_selected_background_production_state(string background_type);
    bool is_selected_background_final_state(string background_type, double high_p0_threshold);
    string get_background_final_state(double high_p0_threshold);
    bool is_vertex_contained(double fv_distance);

    std::vector<G4Particle> get_g4_primary_p(int PDG);
    G4Particle get_g4_particle_mother(G4Particle g4_p);
};

#endif //MUPIPROJECT_TRUEINTERACTION_H
