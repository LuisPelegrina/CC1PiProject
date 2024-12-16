//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#include "Cut_Parameters.h"

void Cut_Parameters::use_default_cut_p() {

  this->min_distance_to_wall_x_y = 20;
  this->min_distance_to_last_z_wall = 20;
  this->min_distance_to_first_z_wall = 30;
  this->min_distance_to_CPA = 5;

  this-> min_crumbs_score = 0;

  this->apply_quality_cuts = true;
  this->quality_cut_min_track_lenght = 3;
  this->quality_cut_min_total_hits = 15;
  this->quality_cut_min_hits_in_two_planes = 5;

  this->max_primary_distance_to_vertex = 10;
  this->min_track_lenght = 3;
  this->min_track_score = 0.45;

  this->min_distance_to_consider_contained = 5;

  this->max_shower_track_score = 0.45;
  this->min_shower_KE = 90;
  this->max_shower_track_score_no_e = 0.45;

  this->razzled_min_muon_like_score = 0;
  this->razzled_max_proton_score = 0.45;
  this->razzled_max_shower_score = 0.65;

  this->proton_rejection_min_proton_score = 80;
  this->proton_rejection_max_muon_score = 60;
  this->proton_rejection_min_TL = 10;

  this->min_muon_candidate_TL = 20;

  this->chi2_min_TL = 12.5;
  this->chi2_min_muon_score = 20;
  this->chi2_max_proton_score = 85;
}