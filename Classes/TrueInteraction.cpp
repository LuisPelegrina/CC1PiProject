//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#include "TrueInteraction.h"


//PRINTERS
void TrueInteraction::print(bool b_print_gen_particles, bool b_print_g4_particles , bool print_vector) {

  cout << "Interaction Mode: " << this->interaction_mode << "Interaction Type: " << this->interaction_type << " CCNC: " << this->CCNC << " E0_Nu: " << this->E0_nu  << " PDG_Nu: " << this->PDG_nu << endl;
  cout << "Target: " << this->target << "HitNuc: " << this->HitNuc << " HitQuark: " << this->HitQuark << endl;
  cout << " W: " << this->W  << " X: " << this->X << " Y: " << this->Y << " QSqr: " <<this->QSqr<< endl;
  cout << " Primary vertex (x,y,z): " << this->primary_vertex.X()  << " " << this->primary_vertex.Y() << " Y: " << this->primary_vertex.Z()<< endl;
  if(b_print_gen_particles) {
        cout << "Gen particles: " << endl;
        print_gen_particles(this->gen_particles, print_vector);
        cout << endl;
    }
  cout << "Primary gen particles: " << endl;
  print_gen_particles(this->gen_primary_particles, print_vector);
  cout << endl;
    if(b_print_g4_particles) {
        cout << "G4 particles: " << endl;
        print_g4_particles(this->g4_particles, print_vector, 0, 0);
        cout << endl;
    }
}

void TrueInteraction::print_gen_particles(std::vector<GeneratorParticle> particles_to_print, bool print_vector) {
    for(int i_p = 0; i_p < int(particles_to_print.size()); i_p++) {
        particles_to_print.at(i_p).print(print_vector);
    }
}

void TrueInteraction::print_g4_particles(std::vector<G4Particle> particles_to_print, bool print_vector, int iter,int  mother_ID) {
    for(int i_p = 0; i_p < int(particles_to_print.size()); i_p++) {
        if(particles_to_print.at(i_p).mother != mother_ID ) continue;
        for(int i = 0; i< iter;i++) cout << "   ";
        particles_to_print.at(i_p).print(print_vector);
        print_g4_particles(particles_to_print, false, iter+1, particles_to_print.at(i_p).ID );
    }
}

//GETTERS
std::vector<GeneratorParticle> TrueInteraction::get_generator_primary_p(int target_PDG) {

    std::vector<GeneratorParticle> primary_particle;

    for(int i_p = 0; i_p < int(this->gen_primary_particles.size());i_p++) {
        GeneratorParticle gen_particle = this->gen_primary_particles.at(i_p);
        if((abs(gen_particle.PDG) == target_PDG) && (gen_particle.is_primary())) primary_particle.push_back(gen_particle);
    }

    return primary_particle;
}

std::vector<GeneratorParticle> TrueInteraction::get_generator_primary_p_sing(int target_PDG) {

    std::vector<GeneratorParticle> primary_particle;

    for(int i_p = 0; i_p < int(this->gen_primary_particles.size());i_p++) {
        GeneratorParticle gen_particle = this->gen_primary_particles.at(i_p);
        if((gen_particle.PDG == target_PDG) && (gen_particle.is_primary()))  primary_particle.push_back(gen_particle);
    }

    return primary_particle;
}


int TrueInteraction::get_generator_num_primaries() {
    return this->gen_primary_particles.size();
}

int TrueInteraction::get_generator_num_primary_p(int target_PDG) {
    int  num_primary = 0;
    for(int i_p = 0; i_p < this->gen_primary_particles.size(); i_p++) {
        GeneratorParticle primary_particle = this->gen_primary_particles.at(i_p);
        if(abs(primary_particle.PDG) == target_PDG) num_primary++;
    }
    return num_primary;
}

int TrueInteraction::get_generator_num_primary_p_sing(int target_PDG) {
    int  num_primary = 0;
    for(int i_p = 0; i_p < this->gen_primary_particles.size(); i_p++) {
        GeneratorParticle primary_particle = this->gen_primary_particles.at(i_p);
        if(primary_particle.PDG == target_PDG) num_primary++;
    }
    return num_primary;
}

int TrueInteraction::get_generator_num_primary_p_high_P0(int target_PDG, double P0_threshold) {
    int  num_primary = 0;
    for(int i_p = 0; i_p < this->gen_primary_particles.size(); i_p++) {
        GeneratorParticle primary_particle = this->gen_primary_particles.at(i_p);
        if(primary_particle.P0.Mag() < P0_threshold) continue;
        if(abs(primary_particle.PDG) == target_PDG) num_primary++;
    }
    return num_primary;
}

std::vector<G4Particle> TrueInteraction::get_g4_primary_p(int PDG) {

    std::vector<G4Particle> primary_particle;

    for(int i_p = 0; i_p < this->g4_particles.size(); i_p++) {
        G4Particle particle = this->g4_particles.at(i_p);

        if(abs(particle.PDG) != PDG) continue;

        if(particle.is_primary()) primary_particle.push_back(particle);
    }

    return primary_particle;
}

//CUTS
bool TrueInteraction::pass_containment_cut(string cut_type, double fv_distance) {
    bool is_vertex_contained = this->is_vertex_contained(fv_distance);
    bool is_contained = false;

    G4Particle primary_muon;
    G4Particle primary_pion;
    if( this->get_g4_primary_p(13).size() != 0) {
        primary_muon = this->get_g4_primary_p(13).at(0);
    }
    if( this->get_g4_primary_p(211).size() != 0) {
        primary_pion = this->get_g4_primary_p(211).at(0);
    }

     if( cut_type.compare("no_cut") == 0){
        is_contained = true;

    } else if(cut_type.compare("vertex_contained") == 0){
        is_contained = is_vertex_contained;

    } else if(cut_type.compare("both_contained_fv") == 0){
        if(primary_pion.is_contained_in_fv(fv_distance) && primary_muon.is_contained_in_fv(fv_distance)) is_contained = true;
    } else if(cut_type.compare("only_one_contained_fv") == 0){
        if(primary_pion.is_contained_in_fv(fv_distance) && !primary_muon.is_contained_in_fv(fv_distance)) is_contained = true;
        if(!primary_pion.is_contained_in_fv(fv_distance) && primary_muon.is_contained_in_fv(fv_distance)) is_contained = true;
    } else if(cut_type.compare("none_contained_fv") == 0){
        if(!(primary_pion.is_contained_in_fv(fv_distance) || primary_muon.is_contained_in_fv(fv_distance))) is_contained = true;
    } else if(cut_type.compare("muon_contained_fv") == 0){
        if(primary_muon.is_contained_in_fv(fv_distance)) is_contained = true;
    } else if(cut_type.compare("pion_contained_fv") == 0){
        if(primary_pion.is_contained_in_fv(fv_distance)) is_contained = true;
    } else if(cut_type.compare("only_muon_contained_fv") == 0){
        if(!primary_pion.is_contained_in_fv(fv_distance) && primary_muon.is_contained_in_fv(fv_distance)) is_contained = true;
    } else if(cut_type.compare("only_pion_contained_fv") == 0){
        if(primary_pion.is_contained_in_fv(fv_distance) && !primary_muon.is_contained_in_fv(fv_distance)) is_contained = true;
    } else {
      cout << "ERROR FINAL STATE CUT" << endl;
      is_contained = false;

    }

    return is_contained;
}

string TrueInteraction::get_background_final_state(double high_p0_threshold) {
  int num_primaries = this->get_generator_num_primaries();
  int num_mu = this->get_generator_num_primary_p(13);
  int num_e = this->get_generator_num_primary_p(11);
  int num_pi = this->get_generator_num_primary_p(211);
  int num_p = this->get_generator_num_primary_p(2212);
  int num_n = this->get_generator_num_primary_p(2112);
  int num_p_high_p0 = this->get_generator_num_primary_p_high_P0(2212, high_p0_threshold);

  int num_O = num_primaries - num_mu - num_pi;
  string state;
  bool is_CC1Pi = ((num_pi == 1) && (num_mu == 1) && (num_O - num_p - num_n == 0));


  if(is_CC1Pi) {
   state = "CC1Pi";
  } else if(!is_CC1Pi && ((num_pi == 1) && (num_mu == 1))){
   state = "other_CC1pi";
  } else if(num_e == 1)  {
   state = "CC_e";
  } else if ((num_mu != 1) && (num_e != 1)) {
    state = "NC";
  } else  if((num_mu == 1) && (num_pi == 0) && (num_p_high_p0 == 0))  {
   state = "CC_mu_0pi_0p";
  } else if((num_mu == 1) && (num_pi == 0) && (num_p_high_p0 == 1)) {
   state = "CC_mu_0pi_1p";
  } else if((num_mu == 1) && (num_pi == 0) && (num_p_high_p0 > 1)){
   state = "CC_mu_0pi_2p";
  } else if((num_mu == 1) && (num_pi > 1)) {
   state = "CC_mu_2pi";
  } else {
    state = "other";
  }


  return state;


}


bool TrueInteraction::is_selected_final_state(string cut_type, double high_p0_threshold) {
    int num_primaries = this->get_generator_num_primaries();
    int num_mu = this->get_generator_num_primary_p(13);
    int num_pi = this->get_generator_num_primary_p(211);
    int num_pi0 = this->get_generator_num_primary_p(111);
    int num_p = this->get_generator_num_primary_p(2212);
    int num_n = this->get_generator_num_primary_p(2112);
    int num_p_high_p0 = this->get_generator_num_primary_p_high_P0(2212, high_p0_threshold);

    int num_O = num_primaries - num_mu - num_pi;
    bool is_CC1Pi = ((num_pi == 1) && (num_mu == 1) && (num_O - num_p - num_n == 0));
    //bool is_CC1Pi = ((num_pi == 1) && (num_mu == 1) && (num_pi0 == 0));
    //bool is_CC1Pi = ((num_pi == 1) && (num_mu == 1));

    if( cut_type.compare("no_cut") == 0){
        return true;

    } else if( cut_type.compare("CC1Pi") == 0){
      if(is_CC1Pi) return true;

    } else if( cut_type.compare("old_CC1Pi") == 0){
      if((num_pi == 1) && (num_mu == 1)) return true;

    }  else if( cut_type.compare("other_CC1pi") == 0){
      if(!is_CC1Pi && ((num_pi == 1) && (num_mu == 1))) return true;

    } else if (cut_type.compare("N_CC1Pi") == 0) {
      if(!is_CC1Pi) return false;
        if((num_O-num_n) == 0) return true;

    } else if (cut_type.compare("1P_CC1Pi") == 0) {
      if(!is_CC1Pi) return false;
        if((num_p_high_p0 == 1) && ((num_O - num_n - num_p) == 0)) return true;

    } else if (cut_type.compare("plus2P_CC1Pi") == 0) {
      if(!is_CC1Pi) return false;
        if((num_p_high_p0 >= 2) && ((num_O - num_n - num_p) == 0) ) return true;

    } else if (cut_type.compare("1P_XO_CC1Pi") == 0) {
      if(!is_CC1Pi) return false;
        if(num_p_high_p0 == 1) return true;

    } else if (cut_type.compare("plus2P_XO_CC1Pi") == 0) {
      if(!is_CC1Pi) return false;
        if(num_p_high_p0 >= 2) return true;

    } else if (cut_type.compare("hasP_CC1Pi") == 0) {
      if(!is_CC1Pi) return false;
        if(num_p > 0) return true;

    } else if(cut_type.compare("Coherent") == 0){
      if(!is_CC1Pi) return false;
        if(num_O == 0) return true;

    }  else {
      cout << "ERROR FINAL STATE CUT: " << cut_type<< endl;
      return false;

    }

    return false;

}

bool TrueInteraction::is_selected_background_final_state(string background_type, double high_p0_threshold) {
    int num_primaries = this->get_generator_num_primaries();
    int num_mu = this->get_generator_num_primary_p(13);
    int num_e = this->get_generator_num_primary_p(11);
    int num_pi = this->get_generator_num_primary_p(211);
    int num_p = this->get_generator_num_primary_p(2212);
    int num_n = this->get_generator_num_primary_p(2112);
    int num_p_high_p0 = this->get_generator_num_primary_p_high_P0(2212, high_p0_threshold);

    int num_O = num_primaries - num_mu - num_pi;
    bool is_CC1Pi = ((num_pi == 1) && (num_mu == 1) && (num_O - num_p - num_n == 0));

    if(background_type.compare("no_cut") == 0){
        return true;

    }else if( background_type.compare("CC1Pi") == 0){
        if(is_CC1Pi) return true;

    } else if( background_type.compare("other_CC1pi") == 0){
      if(!is_CC1Pi && ((num_pi == 1) && (num_mu == 1))) return true;

    } else if (background_type.compare("CC_mu") == 0) {
        if((num_pi != 1) && (num_mu == 1)) return true;

    } else if (background_type.compare("CC_e") == 0) {
        if(num_e == 1) return true;

    } else if (background_type.compare("NC") == 0) {
        if((num_mu != 1) && (num_e != 1)) return true;

    } else if (background_type.compare("CC_mu_0pi") == 0) {
        if((num_mu == 1) && (num_pi == 0))  return true;

    } else if (background_type.compare("CC_mu_0pi_0p") == 0) {
        if((num_mu == 1) && (num_pi == 0) && (num_p_high_p0 == 0))  return true;

    } else if (background_type.compare("CC_mu_0pi_1p") == 0) {
        if((num_mu == 1) && (num_pi == 0) && (num_p_high_p0 == 1))  return true;

    } else if (background_type.compare("CC_mu_0pi_2p") == 0) {
        if((num_mu == 1) && (num_pi == 0) && (num_p_high_p0 > 1))  return true;

    } else if (background_type.compare("CC_mu_2pi") == 0) {
        if((num_mu == 1) && (num_pi > 1))return true;
    } else {
      cout << "ERROR FINAL STATE CUT" << endl;
      cout << background_type << endl;
      return true;
  }

    return false;

}

bool TrueInteraction::is_selected_background_production_state(string background_type) {
    //cout << background_type << endl;
    int num_primaries = this->get_generator_num_primaries();
    int num_mu = this->get_generator_num_primary_p(13);
    int num_e = this->get_generator_num_primary_p(11);
    int num_pi = this->get_generator_num_primary_p(211);
    int num_p = this->get_generator_num_primary_p(2212);
    int num_n = this->get_generator_num_primary_p(2112);
    int num_p_high_p0 = this->get_generator_num_primary_p_high_P0(2212, 0.25);
    int num_n_high_p0 = this->get_generator_num_primary_p_high_P0(2112, 0.25);

    int num_O = num_primaries - num_mu - num_pi;
  bool is_CC1Pi = ((num_pi == 1) && (num_mu == 1) && (num_O - num_p - num_n == 0));

    if(background_type.compare("no_cut") == 0){
        return true;

    } else if(background_type.compare("CC1Pi") == 0){
        if(is_CC1Pi) return true;

    } else if( background_type.compare("other_CC1pi") == 0){
      if(!is_CC1Pi && ((num_pi == 1) && (num_mu == 1))) return true;

    } else if (background_type.compare("CC_mu") == 0) {
        if((num_pi != 1) && (num_mu == 1)) return true;

    } else if (background_type.compare("CC_e") == 0) {
        if(num_e == 1) return true;

    } else if (background_type.compare("NC") == 0) {
        if((num_mu != 1) && (num_e != 1)) return true;

    } else if (background_type.compare("CC_mu_qe") == 0) {
        if((num_pi != 1) && (num_mu == 1) && this->interaction_mode == 0) return true;

    } else if (background_type.compare("CC_mu_dis") == 0) {
        if((num_pi != 1) && (num_mu == 1) && this->interaction_mode == 2) return true;

    } else if (background_type.compare("CC_mu_res") == 0) {
        if((num_pi != 1) && (num_mu == 1) && this->interaction_mode == 1) return true;

    } else if (background_type.compare("CC_mu_other") == 0) {
        if((num_pi != 1) && (num_mu == 1) && (this->interaction_mode) > 3 && ( this->interaction_mode != 10)) return true;

    } else if (background_type.compare("CC_mu_coh") == 0) {
        if((num_pi != 1) && (num_mu == 1) && this->interaction_mode == 3) return true;

    }else if (background_type.compare("CC_mu_mec") == 0) {
        if((num_pi != 1) && (num_mu == 1) && this->interaction_mode == 10) return true;
    } else {
      cout << "NOT DEF" << endl;
      return false;
    }

    return false;

}

bool TrueInteraction::is_vertex_contained(double fv_distance) {
  bool is_contained = true;

  if(this->primary_vertex.X() > 200 - fv_distance) is_contained = false;
  if(this->primary_vertex.X() < -200 + fv_distance) is_contained = false;

  if(this->primary_vertex.Y() > 200 - fv_distance) is_contained = false;
  if(this->primary_vertex.Y() < -200 + fv_distance) is_contained = false;

  if(this->primary_vertex.Z() > 500 - fv_distance) is_contained = false;
  if(this->primary_vertex.Z() < 0 + fv_distance) is_contained = false;

  if ((this->primary_vertex.X() < 0 + fv_distance) && (this->primary_vertex.X() > 0 - fv_distance)) is_contained = false;

  return is_contained;
}


G4Particle TrueInteraction::get_g4_particle_by_id(int id) {
  G4Particle output_p;
  bool found = false;
  for (int i_pg4 = 0; i_pg4 < this->g4_particles.size(); i_pg4++) {
    if(found) continue;
    G4Particle g4_p = this->g4_particles.at(i_pg4);
    if (id == g4_p.ID) {
      output_p = g4_p;
      found = true;
    }
  }
  return output_p;
}

G4Particle TrueInteraction::get_g4_particle_mother(G4Particle g4_p) {
  G4Particle output_p;
  for (int i_pg4 = 0; i_pg4 < this->g4_particles.size(); i_pg4++) {
    if (this->g4_particles.at(i_pg4).ID == g4_p.mother) output_p = this->g4_particles.at(i_pg4);
  }
  return output_p;
}
