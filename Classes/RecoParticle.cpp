//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#include "RecoParticle.h"




bool Reco_Particle::is_not_razzled_proton(Cut_Parameters cut_p) {
    bool is_proton = true;

    int razzled_pdg = this->razzled_score.razzled_pdg;
    if(razzled_pdg == 2212) is_proton = false;

    return is_proton;
}

bool Reco_Particle::is_razzled_muon_like_candidate(Cut_Parameters cut_p) {
    bool is_candidate = false;

    //double muon_like_razzled_score = this->razzled_score.muon_score +  this->razzled_score.pion_score;
    double muon_like_razzled_score =  this->razzled_score.pion_score;
    if(this->razzled_score.muon_score > muon_like_razzled_score) muon_like_razzled_score =  this->razzled_score.muon_score;

    double shower_score = this->razzled_score.electron_score +  this->razzled_score.photon_score;
    double proton_score = this->razzled_score.proton_score;

    //cout << muon_like_dazzle_score << " " <<  reco_particle.dazzle_score.muon_score<< " " << reco_particle.dazzle_score.pion_score<< endl;

    //if((muon_like_dazzle_score > reco_particle.dazzle_score.proton_score) && (reco_particle.dazzle_score.undef_score < cut_p.dazzle_max_other_score)) is_candidate = true;

    if((shower_score < cut_p.razzled_max_shower_score) && (proton_score < cut_p.razzled_max_proton_score) && ( muon_like_razzled_score > cut_p.razzled_min_muon_like_score)) is_candidate = true;

    //if(muon_like_razzled_score > 1-muon_like_razzled_score) is_candidate = true;
    int razzled_pdg = this->razzled_score.razzled_pdg;
    //if(razzled_pdg == 13) is_candidate = true;
    //if(razzled_pdg == 211) is_candidate = true;

    return is_candidate;
}

bool Reco_Particle::is_razzled_muon_like(Cut_Parameters cut_p) {
    bool is_candidate = false;

    int razzled_pdg = this->razzled_score.razzled_pdg;
    if(razzled_pdg == 13) is_candidate = true;
    if(razzled_pdg == 211) is_candidate = true;

    return is_candidate;
}

/*
bool Reco_Particle::is_muon_like_chi2_candidate(Cut_Parameters cut_p) {
    bool is_candidate = false;
    double muon_chi2_score = this->chi2_score.muon_score;

    if((this->track_lenght > cut_p.chi2_min_TL) &&
       (this->chi2_score.proton_score > cut_p.chi2_min_proton_score) &&
       ( muon_chi2_score < cut_p.chi2_max_muon_score)) is_candidate = true;

    return is_candidate;
}
*/

bool Reco_Particle::is_contained_in_fv(Cut_Parameters cut_p) {
    bool is_contained = true;

    if(this->track_start.X() > 200 - cut_p.min_distance_to_wall_x_y) is_contained = false;
    if(this->track_start.X() < -200 + cut_p.min_distance_to_wall_x_y) is_contained = false;

    if(this->track_start.Y() > 200 - cut_p.min_distance_to_wall_x_y) is_contained = false;
    if(this->track_start.Y() < -200 + cut_p.min_distance_to_wall_x_y) is_contained = false;

    if(this->track_start.Z() > 500 - cut_p.min_distance_to_last_z_wall) is_contained = false;
    if(this->track_start.Z() < 0 + cut_p.min_distance_to_first_z_wall) is_contained = false;

    if(this->track_end.X() > 200 - cut_p.min_distance_to_wall_x_y) is_contained = false;
    if(this->track_end.X() < -200 + cut_p.min_distance_to_wall_x_y) is_contained = false;

    if(this->track_end.Y() > 200 - cut_p.min_distance_to_wall_x_y) is_contained = false;
    if(this->track_end.Y() < -200 +cut_p.min_distance_to_wall_x_y) is_contained = false;

    if(this->track_end.Z() > 500 - cut_p.min_distance_to_last_z_wall) is_contained = false;
    if(this->track_end.Z() < 0 + cut_p.min_distance_to_first_z_wall) is_contained = false;

    if ((this->track_start.X() < 0 + cut_p.min_distance_to_CPA) && (this->track_start.X() > 0 - cut_p.min_distance_to_CPA)) is_contained = false;
    if ((this->track_end.X() < 0 + cut_p.min_distance_to_CPA) && (this->track_end.X() > 0 - cut_p.min_distance_to_CPA)) is_contained = false;


    return is_contained;
}


bool Reco_Particle::pass_quality_cuts(Cut_Parameters cut_p, vector<Hit> slice_hits) {
  bool pass_cuts = true;

  int num_hits_U = 0;
  int num_hits_V = 0;
  int num_hits_C = 0;
  int total_hits = 0;
  for(int i_h = 0; i_h < slice_hits.size(); i_h++){
    Hit hit = slice_hits.at(i_h);
    if(hit.associated_pfp_ID == this->ID) {
      if(hit.Plane_ID == 2) num_hits_C++;
      if(hit.Plane_ID == 1) num_hits_U++;
      if(hit.Plane_ID == 0) num_hits_V++;
      if(hit.Plane_ID != -1) total_hits++;
    }
  }
  int planes_ok = 0;
  if(num_hits_U >= cut_p.quality_cut_min_hits_in_two_planes) planes_ok++;
  if(num_hits_V >= cut_p.quality_cut_min_hits_in_two_planes) planes_ok++;
  if(num_hits_C >= cut_p.quality_cut_min_hits_in_two_planes) planes_ok++;

  if((planes_ok < 2) && (num_hits_C < cut_p.quality_cut_min_total_hits)) pass_cuts = false;
  if(total_hits < cut_p.quality_cut_min_total_hits) pass_cuts = false;
  if(this->track_lenght < cut_p.quality_cut_min_track_lenght) pass_cuts = false;

  if(this->chi2_score.proton_score == -1) pass_cuts = false;
  if(this->track_score == -1) pass_cuts = false;
  //if(this->track_lenght > cut_p.min_track_lenght) pass_cuts = false;

  return pass_cuts;
}


bool Reco_Particle::is_track(Cut_Parameters cut_p) {
    if((this->track_score >= cut_p.min_track_score) && (this->track_lenght > cut_p.min_track_lenght)) {
        return true;
    } else {
        return false;
    }
}

bool Reco_Particle::is_shower(Cut_Parameters cut_p, bool check_energy) {
    bool KE = true;
    if(check_energy) KE = this->shower_energy > cut_p.min_shower_KE;

    if((this->track_score < cut_p.max_shower_track_score) && (KE)) {
        return true;
    } else {
        return false;
    }
}


bool Reco_Particle::is_primary(TVector3 vertex_cm, Cut_Parameters cut_p) {
    bool is_primary = this->is_pandora_primary;
    double start_distance = (this->track_start - vertex_cm).Mag();
    if(start_distance > cut_p.max_primary_distance_to_vertex) is_primary = false;
    return is_primary;
}

void Reco_Particle::print(TVector3 vertex_cm) {

    double delta_x = this->track_start.X() - vertex_cm.X();
    double delta_y = this->track_start.Y() - vertex_cm.Y();
    double delta_z = this->track_start.Z() - vertex_cm.Z();

    double distance = sqrt(pow(delta_x,2) + pow(delta_y,2) + pow(delta_z,2));

    cout << "TS: " << this->track_score << " TL: " <<  this->track_lenght << " DtoV: " << distance << " Is_primary: " << this->is_pandora_primary << endl;
    cout << "matched_pdg: " <<  this->matched_pdg << " purity: " <<  this->purity << " completeness: " <<  this->completeness <<endl;
    cout << "START: " << this->track_start.X()<< " " << this->track_start.Y()<< " " << this->track_start.Z() << endl;
    cout << "DtoSTART: " << this->track_start.X() - vertex_cm.X()<< " " << this->track_start.Y() - vertex_cm.Y()<< " " << this->track_start.Z() - vertex_cm.Z()<< endl;
    cout << "END: " << this->track_end.X()<< " " << this->track_end.Y()<< " " << this->track_end.Z() << endl;
    cout << "Wires: " << endl;
    cout << "UP: " << this->track_start_wires.at(0).Wire_ID << "  " << this->track_start_wires.at(0).drift_t  << endl;
    cout << "VP: " << this->track_start_wires.at(1).Wire_ID << "  " << this->track_start_wires.at(1).drift_t  << endl;
    cout << "CP: " << this->track_start_wires.at(2).Wire_ID << "  " << this->track_start_wires.at(2).drift_t  << endl;

    cout << this->razzled_score.razzled_pdg << endl;
    cout << "Chi2: " << endl;
    cout << "2212: " << this->chi2_score.proton_score << " 321: " << this->chi2_score.kaon_score
         << " 211: " << this->chi2_score.pion_score << " 13: " << this->chi2_score.muon_score << endl;
    cout << "Razzled: " << endl;
    cout << "22: " << this->razzled_score.photon_score << " 11: " << this->razzled_score.electron_score
         << " 2212: " << this->razzled_score.proton_score << " 211: " << this->razzled_score.pion_score
         << " 13: " << this->razzled_score.muon_score << endl;
}
