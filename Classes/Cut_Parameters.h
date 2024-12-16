//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_CUT_PARAMETERS_H
#define MUPIPROJECT_CUT_PARAMETERS_H



class Cut_Parameters {
public:
    double min_track_score;
    double quality_cut_min_track_lenght;
    double min_track_lenght;
    double max_primary_distance_to_vertex;

    double min_crumbs_score;

    int quality_cut_min_total_hits;
    int quality_cut_min_hits_in_two_planes;
    bool apply_quality_cuts;

    double min_distance_to_CPA;
    double min_distance_to_wall_x_y;
    double min_distance_to_first_z_wall;
    double min_distance_to_last_z_wall;

    double razzled_min_muon_like_score;
    double razzled_max_proton_score;
    double razzled_max_shower_score;

    double min_distance_to_consider_contained;

    double max_shower_track_score;
    double min_shower_KE;
    double max_shower_track_score_no_e;

    double proton_rejection_min_proton_score;
    double proton_rejection_max_muon_score;
    double proton_rejection_min_TL;

     double min_muon_candidate_TL;

     double chi2_min_TL;
     double chi2_min_muon_score;
     double chi2_max_proton_score;
    void use_default_cut_p();

    Cut_Parameters() { // Constructor with parameters
      min_track_score = 0.5;
      min_track_lenght = 3;

      quality_cut_min_track_lenght = 3;
      quality_cut_min_total_hits = 15;
      quality_cut_min_hits_in_two_planes = 5;
      apply_quality_cuts = true;

      //PIDS
      /*
      chi2_max_muon_like_score = 80;
      chi2_min_proton_score = 45;
      chi2_min_TL = 0;
      */

      razzled_min_muon_like_score = 0;
      razzled_max_proton_score = 0.45;
      razzled_max_shower_score = 0.65;

      max_primary_distance_to_vertex = 10;

      min_distance_to_CPA = 10; // 5 HNL
      min_distance_to_wall_x_y = 20; // 20 HNL
      min_distance_to_first_z_wall = 20; // 10 HNL
      min_distance_to_last_z_wall = 20; //50 HNL
      min_crumbs_score = 0;

      min_distance_to_consider_contained = 5; //50 HNL

      max_shower_track_score = 0.45;
      min_shower_KE = 80;
      max_shower_track_score_no_e = 0.4;

      proton_rejection_min_proton_score = 85;
      proton_rejection_max_muon_score = 60;
      proton_rejection_min_TL = 10;
      min_muon_candidate_TL = 20;

      chi2_min_TL = 10;
      chi2_min_muon_score = 2000;
      chi2_max_proton_score = 70;
    }
};


#endif //MUPIPROJECT_CUT_PARAMETERS_H
