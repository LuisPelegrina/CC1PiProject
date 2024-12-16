//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//


#include "Slice.h"

int Slice::correct_Wire_ID(int input_wire_ID, int Plane_ID, int TPC_ID, int Channel_ID) {
    int corrected_Wire_ID = input_wire_ID;

    if(TPC_ID == 0) {
        if(Plane_ID == 2) {
            if(Channel_ID > 4399) corrected_Wire_ID -= 1;
            if(Channel_ID > 5221) corrected_Wire_ID -= 1;
            if(Channel_ID > 4802) corrected_Wire_ID -= 1;
            if(Channel_ID > 4803) corrected_Wire_ID -= 1;
            if(Channel_ID > 4804) corrected_Wire_ID -= 1;
            if(Channel_ID > 4805) corrected_Wire_ID -= 1;
            if(Channel_ID > 4806) corrected_Wire_ID -= 1;
            if(Channel_ID > 4807) corrected_Wire_ID -= 1;
        }
    } else if(TPC_ID == 1) {
        if(Plane_ID == 2) {
            if(Channel_ID > 10037) corrected_Wire_ID -= 1;
            if(Channel_ID > 10859) corrected_Wire_ID -= 1;
            if(Channel_ID > 10440) corrected_Wire_ID -= 1;
            if(Channel_ID > 10441) corrected_Wire_ID -= 1;
            if(Channel_ID > 10442) corrected_Wire_ID -= 1;
            if(Channel_ID > 10443) corrected_Wire_ID -= 1;
            if(Channel_ID > 10444) corrected_Wire_ID -= 1;
            if(Channel_ID > 10445) corrected_Wire_ID -= 1;
        }
    }
    return corrected_Wire_ID;
}

void Slice::print(Cut_Parameters cut_p, bool b_print_ori_file_name, bool b_print_gen_particles, bool b_print_g4_particles , bool print_vector, bool b_print_hits, bool b_print_reco_vertex) {
    if(b_print_ori_file_name) cout << "File: "<< this->file_name << endl;
    cout << "eventID: " << this->event_ID << " slice_ID: " << this->slice_ID << " Weight: "<<  this-> weight << " gen_index: "<< this->gen_index <<  endl;

    this->true_interaction.print(b_print_gen_particles, b_print_g4_particles, print_vector);

    cout << "IsClearCosmic: " << this->is_clear_cosmic << " IsReconstructed: " << this->is_reconstructed << endl;
    cout << "NuScore: " << this->nu_score << " CrumbsScore: " << this->crumbs_score << endl;

    if(b_print_reco_vertex) {
        cout << "Vertex true: " << endl;
        cout << "time: " << this->t0 << " Score: "<< this->t0_score << endl;
        cout << "cm: " << this->primary_vertex_true.vertex_cm.X() << " " << this->primary_vertex_true.vertex_cm.Y() << " " << this->primary_vertex_true.vertex_cm.Z() << endl;

        cout << "Wires: " << endl;
        cout << "UP: " << this->primary_vertex_true.vertex_wire.at(0).Wire_ID << "  " << this->primary_vertex_true.vertex_wire.at(0).drift_t  << endl;
        cout << "VP: " << this->primary_vertex_true.vertex_wire.at(1).Wire_ID << "  " << this->primary_vertex_true.vertex_wire.at(1).drift_t  << endl;
        cout << "CP: " << this->primary_vertex_true.vertex_wire.at(2).Wire_ID << "  " << this->primary_vertex_true.vertex_wire.at(2).drift_t  << endl;

        cout << "Vertex reco: " << endl;
        cout << "time: " << this->t0 << " Score: "<< this->t0_score << endl;
        cout << "cm: " << this->primary_vertex_reco.vertex_cm.X() << " " << this->primary_vertex_reco.vertex_cm.Y() << " " << this->primary_vertex_reco.vertex_cm.Z() << endl;

        cout << "Wires: " << endl;
        cout << "UP: " << this->primary_vertex_reco.vertex_wire.at(0).Wire_ID << "  " << this->primary_vertex_reco.vertex_wire.at(0).drift_t  << endl;
        cout << "VP: " << this->primary_vertex_reco.vertex_wire.at(1).Wire_ID << "  " << this->primary_vertex_reco.vertex_wire.at(1).drift_t  << endl;
        cout << "CP: " << this->primary_vertex_reco.vertex_wire.at(2).Wire_ID << "  " << this->primary_vertex_reco.vertex_wire.at(2).drift_t  << endl;

        cout << "Wires (corr): " << endl;
        cout << "UP: " << this->primary_vertex_reco.corrected_vertex_wire.at(0).Wire_ID << "  " << this->primary_vertex_reco.corrected_vertex_wire.at(0).drift_t << endl;
        cout << "VP: " << this->primary_vertex_reco.corrected_vertex_wire.at(1).Wire_ID << "  " << this->primary_vertex_reco.corrected_vertex_wire.at(1).drift_t << endl;
        cout << "CP: " << this->primary_vertex_reco.corrected_vertex_wire.at(2).Wire_ID << "  " << this->primary_vertex_reco.corrected_vertex_wire.at(2).drift_t << endl;
    }

    if(b_print_hits) {
        cout << "Hits: " << endl;
        for(int i_h; i_h < this->hits.size() ;i_h++) {
            cout << "Integral: " << hits.at(i_h).integral << " +- " << hits.at(i_h).integral_sigma << " Drift_t: " << hits.at(i_h).drift_t  << " +- "<< hits.at(i_h).drift_t_sigma;
            cout <<" TPCID: "<< hits.at(i_h).TPC_ID << " PlaneID: "<< hits.at(i_h).Plane_ID << " WireID: " << hits.at(i_h).Wire_ID << endl;
        }
    }
    cout << "Particle charac:" << endl;
    for(int i_p = 0; i_p < this->pandora_particle.size(); i_p++) {

        double delta_x = this->pandora_particle.at(i_p).track_start.X() - this->primary_vertex_reco.vertex_cm.X();
        double delta_y = this->pandora_particle.at(i_p).track_start.Y() - this->primary_vertex_reco.vertex_cm.Y();
        double delta_z = this->pandora_particle.at(i_p).track_start.Z() - this->primary_vertex_reco.vertex_cm.Z();

        double distance = sqrt(pow(delta_x,2) + pow(delta_y,2) + pow(delta_z,2));

        cout << "TS: " << this->pandora_particle.at(i_p).track_score << " TL: " <<  this->pandora_particle.at(i_p).track_lenght << " DtoV: " << distance << " Is_primary: " << this->pandora_particle.at(i_p).is_pandora_primary << endl;
        cout << "matched_pdg: " <<  this->pandora_particle.at(i_p).matched_pdg << " purity: " <<  this->pandora_particle.at(i_p).purity << " completeness: " <<  this->pandora_particle.at(i_p).completeness <<endl;
        cout << "START: " << this->pandora_particle.at(i_p).track_start.X()<< " " << this->pandora_particle.at(i_p).track_start.Y()<< " " << this->pandora_particle.at(i_p).track_start.Z() << endl;
        cout << "DtoSTART: " << this->pandora_particle.at(i_p).track_start.X() - this->primary_vertex_reco.vertex_cm.X()<< " " << this->pandora_particle.at(i_p).track_start.Y() - this->primary_vertex_reco.vertex_cm.Y()<< " " << this->pandora_particle.at(i_p).track_start.Z() - this->primary_vertex_reco.vertex_cm.Z()<< endl;
        cout << "END: " << this->pandora_particle.at(i_p).track_end.X()<< " " << this->pandora_particle.at(i_p).track_end.Y()<< " " << this->pandora_particle.at(i_p).track_end.Z() << endl;

        cout << "Chi2: " << endl;
        cout << "2212: " << this->pandora_particle.at(i_p).chi2_score.proton_score << " 321: " << this->pandora_particle.at(i_p).chi2_score.kaon_score
             << " 211: " << this->pandora_particle.at(i_p).chi2_score.pion_score << " 13: " << this->pandora_particle.at(i_p).chi2_score.muon_score << endl;
        cout << "Razzled: " << endl;
        cout << "22: " << this->pandora_particle.at(i_p).razzled_score.photon_score << " 11: " << this->pandora_particle.at(i_p).razzled_score.electron_score
             << " 2212: " << this->pandora_particle.at(i_p).razzled_score.proton_score << " 211: " << this->pandora_particle.at(i_p).razzled_score.pion_score
             << " 13: " << this->pandora_particle.at(i_p).razzled_score.muon_score << endl;
        cout << endl;
    }

    cout << "primary Particle charac:" << endl;
    for(Reco_Particle reco_particle: this->pandora_particle) {
      if((!reco_particle.is_track(cut_p))) continue;
      if(!reco_particle.is_primary(this->primary_vertex_reco.vertex_cm, cut_p)) continue;
      if(!reco_particle.pass_quality_cuts(cut_p, this->hits)) continue;
      if(((reco_particle.chi2_score.proton_score < cut_p.proton_rejection_min_proton_score) || (reco_particle.track_lenght < cut_p.proton_rejection_min_TL))) continue;
      if(reco_particle.track_lenght < cut_p.proton_rejection_min_TL) continue;
      reco_particle.print(this->primary_vertex_reco.vertex_cm);
        cout << endl;
    }
    cout << endl;
}

//SETTERS

void Slice::set_weight(double weight, double expected_POT, double sim_POT) {
    double w = 1;

    if(weight == 0) weight = 1;
    w = expected_POT*weight/sim_POT;

    this->weight = w;
}

void Slice::set_corrected_vertex_wire() {
    int nearest_hit_id[3] = {0, 0, 0};
    double min_distance[3] = {20000,20000,20000};
    int num_planes = 3;

    for(int i_h = 0; i_h < this->hits.size(); i_h++) {
        int hit_plane_ID = this->hits.at(i_h).Plane_ID;
        if(hit_plane_ID == -1) hit_plane_ID = 0;

        double delta_t = pow(this->primary_vertex_reco.vertex_wire.at(hit_plane_ID).drift_t/4 - this->hits.at(i_h).drift_t/4, 2);
        double delta_pos = pow(this->primary_vertex_reco.vertex_wire.at(hit_plane_ID).Wire_ID - this->hits.at(i_h).Wire_ID, 2);

        if(min_distance[hit_plane_ID] > sqrt(delta_pos + delta_t)) {
            min_distance[hit_plane_ID] = sqrt(delta_pos + delta_t);
            nearest_hit_id[hit_plane_ID] = i_h;
        }
    }

    for(int i_p = 0; i_p < num_planes; i_p++) {
        this->primary_vertex_reco.corrected_vertex_wire.at(i_p).Wire_ID = this->hits.at(nearest_hit_id[i_p]).Wire_ID;
        this->primary_vertex_reco.corrected_vertex_wire.at(i_p).drift_t = this->hits.at(nearest_hit_id[i_p]).drift_t;
        this->primary_vertex_reco.corrected_vertex_wire.at(i_p).Channel_ID = this->primary_vertex_reco.vertex_wire.at(i_p).Channel_ID;
    }
}

void Slice::set_corrected_hits() {

  std::vector<Hit> corrected_hits;

  for (int i_h = 0; i_h < this->hits.size(); i_h++) {
    Hit hit = hits.at(i_h);

    int hit_plane_ID = hit.Plane_ID;

    double distance = get_distance_to_vertex(hit);

    if (distance < 50) corrected_hits.push_back(hit);
  }

  this->hits = corrected_hits;
}

//GETTERS
double Slice::get_distance_to_vertex(Hit hit) {

    int hit_plane_ID = hit.Plane_ID;

    double delta_t = pow(this->primary_vertex_reco.corrected_vertex_wire.at(hit_plane_ID).drift_t/4 - hit.drift_t/4, 2);
    double delta_pos = pow(this->primary_vertex_reco.corrected_vertex_wire.at(hit_plane_ID).Wire_ID - hit.Wire_ID, 2);

    double distance = sqrt(delta_t + delta_pos);

    return distance;
}

double Slice::get_corrected_time() {

    double c = 29.9792458;
    double z = this->primary_vertex_reco.vertex_cm.Z();
    double time_of_flight = z/c;

    double t_corr = this->t0 - time_of_flight;

    return t_corr;
}

int Slice::get_num_razzled_particle(int PDG, Cut_Parameters cut_p, bool check_if_primary, bool check_if_track, bool check_quality, bool do_proton_rejection) {
    int num_particle = 0;

    for(int i_p = 0; i_p < this->pandora_particle.size();i_p++) {
      Reco_Particle reco_particle = pandora_particle.at(i_p);
      if((check_if_track) && (!reco_particle.is_track(cut_p))) continue;
      if((check_if_primary) && !reco_particle.is_primary(this->primary_vertex_reco.vertex_cm, cut_p)) continue;
      if((check_quality) && (!reco_particle.pass_quality_cuts(cut_p, this->hits))) continue;
      if(do_proton_rejection && ((reco_particle.chi2_score.proton_score < cut_p.proton_rejection_min_proton_score) || (reco_particle.track_lenght < cut_p.proton_rejection_min_TL))) continue;

      if(reco_particle.razzled_score.razzled_pdg == PDG) num_particle++;;

    }
    return num_particle;
}

double Slice::get_num_showers(Cut_Parameters cut_p, bool check_if_primary, bool check_energy) {
    int num_shower = 0;
    for(int i_p = 0; i_p < this->pandora_particle.size();i_p++) {
        Reco_Particle reco_particle = pandora_particle.at(i_p);

        if(!reco_particle.is_shower(cut_p, check_energy)) continue;
        if((check_if_primary) && !reco_particle.is_primary(this->primary_vertex_reco.vertex_cm, cut_p)) continue;

        num_shower++;
    }
    return num_shower;
}

double Slice::get_num_tracks(Cut_Parameters cut_p, bool check_if_primary, bool check_quality) {
    int num_track = 0;
    for(int i_p = 0; i_p < this->pandora_particle.size();i_p++) {
        Reco_Particle reco_particle = pandora_particle.at(i_p);

        if(!reco_particle.is_track(cut_p)) continue;
        if((check_if_primary) && !reco_particle.is_primary(this->primary_vertex_reco.vertex_cm, cut_p)) continue;
        if((check_quality) && (!reco_particle.pass_quality_cuts(cut_p, this->hits))) continue;
        num_track++;
    }
    return num_track;
}

//CUTS
bool Slice::pass_analysis_cut(Cut_Parameters cut_p, string cut_name) {

    if(cut_name.compare("no_cut") == 0){
      return true;
    } else if (cut_name.compare("reco_cut") == 0) {
        return this->is_reconstructed;
    } else if (cut_name.compare("is_clear_cosmic_cut") == 0) {
        return !this->is_clear_cosmic;
    } else if (cut_name.compare("fv_cut") == 0) {
        return this->is_reco_vertex_inside_FV(cut_p);
    } else if (cut_name.compare("shower_cut") == 0) {
      bool pass = false;
      if(this->get_num_showers(cut_p, false, true) < 1) pass = true;
      return pass;
    } else if (cut_name.compare("shower_cut_energy") == 0) {
      bool pass = true;
      Cut_Parameters cut_p_2 = cut_p;
      cut_p_2.max_shower_track_score = cut_p.max_shower_track_score_no_e;
      int num_showers = this->get_num_showers(cut_p_2, false, false);
      int num_showers_energy =  this->get_num_showers(cut_p, false, true);
      if(num_showers > 1) pass = false;
      if((num_showers < 2) && (num_showers_energy > 0)) pass = false;
      return pass;
    } else if (cut_name.compare("track_cut") == 0) {
        if(this->get_num_tracks(cut_p, true, cut_p.apply_quality_cuts) > 1) return true;
    } else if(cut_name == "at_most_one_pfp_scapping_cut") {
      int num_exiting = 0;
      for(Reco_Particle reco_particle: this->pandora_particle) {
        TVector3 track_start = reco_particle.track_start;
        TVector3 track_end = reco_particle.track_end;

        bool is_contained = true;
        if(track_start.X() > 200 - cut_p.min_distance_to_consider_contained) is_contained = false;
        if(track_start.X() < -200 + cut_p.min_distance_to_consider_contained) is_contained = false;

        if(track_start.Y() > 200 - cut_p.min_distance_to_consider_contained) is_contained = false;
        if(track_start.Y() < -200 + cut_p.min_distance_to_consider_contained) is_contained = false;

        if(track_start.Z() > 500 - cut_p.min_distance_to_consider_contained) is_contained = false;
        if(track_start.Z() < 0 + cut_p.min_distance_to_consider_contained) is_contained = false;

        if(track_end.X() > 200 - cut_p.min_distance_to_consider_contained) is_contained = false;
        if(track_end.X() < -200 + cut_p.min_distance_to_consider_contained) is_contained = false;

        if(track_end.Y() > 200 - cut_p.min_distance_to_consider_contained) is_contained = false;
        if(track_end.Y() < -200 +cut_p.min_distance_to_consider_contained) is_contained = false;

        if(track_end.Z() > 500 - cut_p.min_distance_to_consider_contained) is_contained = false;
        if(track_end.Z() < 0 + cut_p.min_distance_to_consider_contained) is_contained = false;

        if(!is_contained) num_exiting++;
      }
      if(num_exiting < 2) return true;
    } else if (cut_name.compare("crumbs_cut") == 0) {
        if(this->crumbs_score > cut_p.min_crumbs_score) return true;
    } else if (cut_name.compare("razzled_muon_like_cut") == 0) {
      int num_muon_like = this->get_num_razzled_particle(211,cut_p, true, true, cut_p.apply_quality_cuts, false) + this->get_num_razzled_particle(13,cut_p, true, true, cut_p.apply_quality_cuts, false);
      if(num_muon_like == 2)  {
        return true;
      } else {
        return false;
      }
    } else if(cut_name == "proton_rejection_razzled") {
      int num_muon_like = this->get_num_razzled_particle(211,cut_p, true, true, cut_p.apply_quality_cuts, true) + this->get_num_razzled_particle(13,cut_p, true, true, cut_p.apply_quality_cuts, true);
      if(num_muon_like == 2)  {
        return true;
      } else {
        return false;
      }
    } else if(cut_name == "chi2_muon_like_cut_harsh") {
      int num_muon_like = 0;
      for(Reco_Particle reco_particle: this->pandora_particle) {
        if((!reco_particle.is_track(cut_p))) continue;
        if(!reco_particle.is_primary(this->primary_vertex_reco.vertex_cm, cut_p)) continue;
        if(!reco_particle.pass_quality_cuts(cut_p, this->hits)) continue;
        double muon_like_score = max(reco_particle.chi2_score.muon_score, reco_particle.chi2_score.pion_score);
        if((muon_like_score > reco_particle.chi2_score.proton_score) && (muon_like_score > reco_particle.chi2_score.kaon_score)) num_muon_like++;
      }
      if(num_muon_like == 2)  {
        return true;
      } else {
        return false;
      }
    } else if(cut_name == "chi2_muon_like_cut") {
      int num_muon_like = 0;


      for(Reco_Particle reco_particle: this->pandora_particle) {
        if((!reco_particle.is_track(cut_p))) continue;
        if(!reco_particle.is_primary(this->primary_vertex_reco.vertex_cm, cut_p)) continue;
        if(!reco_particle.pass_quality_cuts(cut_p, this->hits)) continue;

        if(reco_particle.track_lenght < cut_p.chi2_min_TL) continue;
        if(reco_particle.chi2_score.muon_score > cut_p.chi2_min_muon_score) continue;
        if(reco_particle.chi2_score.proton_score < cut_p.chi2_max_proton_score) continue;
        num_muon_like++;
      }
      if(num_muon_like == 2)  {
        return true;
      } else {
        return false;
      }
    } else if(cut_name == "proton_rejection_razzled_muon_known") {
      int num_muon_candidate = 0;
      for(Reco_Particle reco_particle: this->pandora_particle) {
        if((!reco_particle.is_track(cut_p))) continue;
        if(!reco_particle.is_primary(this->primary_vertex_reco.vertex_cm, cut_p)) continue;
        if(!reco_particle.pass_quality_cuts(cut_p, this->hits)) continue;
        if(((reco_particle.chi2_score.proton_score < cut_p.proton_rejection_min_proton_score) || (reco_particle.track_lenght < cut_p.proton_rejection_min_TL))) continue;
        if(reco_particle.track_lenght < cut_p.proton_rejection_min_TL) continue;
        num_muon_candidate++;
      }
      if(num_muon_candidate < 1) return false;

      int num_muon_like = this->get_num_razzled_particle(211,cut_p, true, true, cut_p.apply_quality_cuts, true) + this->get_num_razzled_particle(13,cut_p, true, true, cut_p.apply_quality_cuts, true);
        if(num_muon_like == 2)  {
          return true;
        } else {
          return false;
        }

    } else {
        cout << "ERROR FINAL STATE CUT " <<  cut_name.c_str() << endl;
        return false;
    }



    return false;

}

bool Slice::is_reco_vertex_inside_FV(Cut_Parameters cut_p) {
    bool is_inside = true;

    TVector3 vertex = this->primary_vertex_reco.vertex_cm;

    if (vertex.X() > 200 - cut_p.min_distance_to_wall_x_y) is_inside = false;
    if (vertex.X() < -200 + cut_p.min_distance_to_wall_x_y) is_inside = false;

    if ((vertex.X() < 0 + cut_p.min_distance_to_CPA) && (vertex.X() > 0 - cut_p.min_distance_to_CPA)) is_inside = false;

    if (vertex.Y() > 200 - cut_p.min_distance_to_wall_x_y) is_inside = false;
    if (vertex.Y() < -200 + cut_p.min_distance_to_wall_x_y) is_inside = false;

    if (vertex.Z() > 500 - cut_p.min_distance_to_last_z_wall) is_inside = false;
    if (vertex.Z() < 0 + cut_p.min_distance_to_first_z_wall) is_inside = false;

    return is_inside;
}

bool Slice::pass_containment_cut(string cut_type, double fv_distance, Cut_Parameters cut_p) {
    bool is_contained = false;

   if(cut_type.compare("no_cut") == 0){
        is_contained = true;

    } else if(cut_type.compare("vertex_contained") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type , fv_distance);

    } else if(cut_type.compare("all_contained") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type, fv_distance);
    } else if(cut_type.compare("both_contained_fv") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type, fv_distance);
    } else if(cut_type.compare("only_one_contained_fv") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type, fv_distance);
    } else if(cut_type.compare("none_contained_fv") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type, fv_distance);
    } else if(cut_type.compare("muon_contained_fv") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type, fv_distance);
    } else if(cut_type.compare("pion_contained_fv") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type, fv_distance);
    } else if(cut_type.compare("only_muon_contained_fv") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type, fv_distance);
    } else if(cut_type.compare("only_pion_contained_fv") == 0) {
        is_contained = this->true_interaction.pass_containment_cut(cut_type, fv_distance);
    } else if(cut_type.compare("both_candidates_contained_fv") == 0) {
        int num_candidates = 0;
        for(int i_p = 0; i_p < this->pandora_particle.size();i_p++) {
            Reco_Particle reco_particle = pandora_particle.at(i_p);

          //  if(!reco_particle.is_muon_like_candidate(this->primary_vertex_reco.vertex_cm, cut_p)) continue;
            if(reco_particle.is_contained_in_fv(cut_p)) num_candidates++;
        }

        if(num_candidates == 2) is_contained =  true;
    } else if(cut_type.compare("one_candidate_contained_fv") == 0) {
        int num_candidates = 0;
        for(int i_p = 0; i_p < this->pandora_particle.size();i_p++) {
            Reco_Particle reco_particle = pandora_particle.at(i_p);

           // if(!reco_particle.is_muon_like_candidate(this->primary_vertex_reco.vertex_cm,cut_p)) continue;
            if(reco_particle.is_contained_in_fv(cut_p)) num_candidates++;
        }

        if(num_candidates == 1) is_contained =  true;


    } else if(cut_type.compare("no_candidates_contained_fv") == 0) {
        int num_candidates = 0;
        for(int i_p = 0; i_p < this->pandora_particle.size();i_p++) {
            Reco_Particle reco_particle = pandora_particle.at(i_p);

            //if(!reco_particle.is_muon_like_candidate(this->primary_vertex_reco.vertex_cm,cut_p)) continue;
            if(reco_particle.is_contained_in_fv(cut_p)) num_candidates++;
        }

        if(num_candidates == 0) is_contained =  true;
    } else {
       cout << "ERROR FINAL STATE CUT: " << cut_type<< endl;
       is_contained = true;
   }

       return is_contained;
}

TH1* Slice::build_dQdx_hist(std::vector<RecoVertexWire> reco_vertex, int target_plane_ID ,double histogram_low_edge, double histogram_high_edge, int id, int num_signal_points) {

    double delta_t;
    double delta_wire;
    double lower_bound;
    double upper_bound;
    double max_distance;
    double min_distance;
    double distance;

    int num_bins = int(histogram_high_edge - histogram_low_edge);
    TH1* h_dummy = new TH1D(("dQ" +to_string(id) ).c_str(), "dQ", num_bins, histogram_low_edge, histogram_high_edge);

    double delta_x;
    double sigma;
    double mean;
    double area;
    double cumulative_charge;

    //const int num_signal_points = 100;


    double vertex_drift_t = reco_vertex.at(target_plane_ID).drift_t;
    double vertex_Wire_ID = reco_vertex.at(target_plane_ID).Wire_ID;

    for(int i_h = 0; i_h < this->hits.size(); i_h++) {
        //Comprueba la distancia al hit, si es muy grande te lo saltas
        Hit hit = this->hits.at(i_h);
        int hit_plane_ID = hit.Plane_ID;

        if(hit_plane_ID != target_plane_ID) continue;

        lower_bound = hit.drift_t/4 - 3 * hit.drift_t_sigma/4;
        upper_bound = hit.drift_t/4 + 3 * hit.drift_t_sigma/4;

        delta_wire = pow(vertex_Wire_ID - hit.Wire_ID, 2);
        delta_t = pow(vertex_drift_t/4 - lower_bound, 2);
        max_distance = sqrt(delta_t + delta_wire);

        delta_wire = pow(vertex_Wire_ID - hit.Wire_ID, 2);
        delta_t = pow(vertex_drift_t/4 - upper_bound, 2);
        distance = sqrt(delta_t + delta_wire);

        if(distance >= max_distance) {
            min_distance = max_distance;
            max_distance = distance;
        } else {
            min_distance = distance;
        }

        if((min_distance > histogram_high_edge) || (max_distance < histogram_low_edge)) continue;

        //Crea la gaussiana y guarda la información
        delta_x = (upper_bound - lower_bound)/(num_signal_points);
        sigma = hit.drift_t_sigma/4;
        mean = hit.drift_t/4;
        area = hit.integral;

        double value[num_signal_points];
        double x[num_signal_points];

        cumulative_charge = 0;
        for(int i_s = 0; i_s < num_signal_points; i_s++) {
            x[i_s] = lower_bound + delta_x / 2 + i_s * delta_x;
            value[i_s] = area * ROOT::Math::normal_pdf(x[i_s], sigma, mean) * delta_x;
            cumulative_charge += value[i_s];
        }

        //Reescala la gaussiana para no perder la carga de los extremos
        for(int i_s = 0; i_s < num_signal_points; i_s++) {
            value[i_s] *= hit.integral/cumulative_charge;
        }

        //Calcula las distancias y guardalas en el histograma
        for(int i_s = 0; i_s < num_signal_points; i_s++) {
            delta_wire = pow(vertex_Wire_ID - hit.Wire_ID, 2);
            delta_t = pow(vertex_drift_t/4 - x[i_s], 2);

            distance = sqrt(delta_t + delta_wire);
            h_dummy->Fill(distance, value[i_s]);
        }
    }
    return h_dummy;

}

TH1* Slice::build_track_dQdx_hist(std::vector<RecoVertexWire> reco_vertex, int target_plane_ID ,double histogram_low_edge, double histogram_high_edge, int id, int num_signal_points, int pfp_ID) {

  double delta_t;
  double delta_wire;
  double lower_bound;
  double upper_bound;
  double max_distance;
  double min_distance;
  double distance;

  int num_bins = int(histogram_high_edge - histogram_low_edge);
  TH1* h_dummy = new TH1D(("dQ" +to_string(id) ).c_str(), "dQ", num_bins, histogram_low_edge, histogram_high_edge);

  double delta_x;
  double sigma;
  double mean;
  double area;
  double cumulative_charge;

  //const int num_signal_points = 100;


  double vertex_drift_t = reco_vertex.at(target_plane_ID).drift_t;
  double vertex_Wire_ID = reco_vertex.at(target_plane_ID).Wire_ID;

  for(int i_h = 0; i_h < this->hits.size(); i_h++) {
    //Comprueba la distancia al hit, si es muy grande te lo saltas
    Hit hit = this->hits.at(i_h);
    if(hit.associated_pfp_ID != pfp_ID) continue;
    int hit_plane_ID = hit.Plane_ID;

    if(hit_plane_ID != target_plane_ID) continue;

    lower_bound = hit.drift_t/4 - 3 * hit.drift_t_sigma/4;
    upper_bound = hit.drift_t/4 + 3 * hit.drift_t_sigma/4;

    delta_wire = pow(vertex_Wire_ID - hit.Wire_ID, 2);
    delta_t = pow(vertex_drift_t/4 - lower_bound, 2);
    max_distance = sqrt(delta_t + delta_wire);

    delta_wire = pow(vertex_Wire_ID - hit.Wire_ID, 2);
    delta_t = pow(vertex_drift_t/4 - upper_bound, 2);
    distance = sqrt(delta_t + delta_wire);

    if(distance >= max_distance) {
      min_distance = max_distance;
      max_distance = distance;
    } else {
      min_distance = distance;
    }

    if((min_distance > histogram_high_edge) || (max_distance < histogram_low_edge)) continue;

    //Crea la gaussiana y guarda la información
    delta_x = (upper_bound - lower_bound)/(num_signal_points);
    sigma = hit.drift_t_sigma/4;
    mean = hit.drift_t/4;
    area = hit.integral;

    double value[num_signal_points];
    double x[num_signal_points];

    cumulative_charge = 0;
    for(int i_s = 0; i_s < num_signal_points; i_s++) {
      x[i_s] = lower_bound + delta_x / 2 + i_s * delta_x;
      value[i_s] = area * ROOT::Math::normal_pdf(x[i_s], sigma, mean) * delta_x;
      cumulative_charge += value[i_s];
    }

    //Reescala la gaussiana para no perder la carga de los extremos
    for(int i_s = 0; i_s < num_signal_points; i_s++) {
      value[i_s] *= hit.integral/cumulative_charge;
    }

    //Calcula las distancias y guardalas en el histograma
    for(int i_s = 0; i_s < num_signal_points; i_s++) {
      delta_wire = pow(vertex_Wire_ID - hit.Wire_ID, 2);
      delta_t = pow(vertex_drift_t/4 - x[i_s], 2);

      distance = sqrt(delta_t + delta_wire);
      h_dummy->Fill(distance, value[i_s]);
    }
  }
  return h_dummy;
}

dQdxSegmentInfo Slice::get_dQdx_segment_info(TH1* h_dQdx, double segment_min_distance, double segment_max_distance) {

    dQdxSegmentInfo segment_info;

    TH1* h_segment =  new TH1D("a", " ", 10000, 0, 500000);


    int low_bin = int(segment_min_distance)+1;
    int high_bin = int(segment_max_distance)+1;

    double max = 0;
    double min = 1e10;
    bool first = true;
    for(int j = low_bin; j < high_bin;j++) {
        //cout << j << " "<< h_segment->GetBinLowEdge(j) << " " << h_segment->GetBinWidth(j) << endl;
        h_segment->Fill(h_dQdx->GetBinContent(j));

        if(h_dQdx->GetBinContent(j) > max) max = h_dQdx->GetBinContent(j);
        if(h_dQdx->GetBinContent(j) < min) min = h_dQdx->GetBinContent(j);
    }

    segment_info.mean = h_segment->GetMean();
    segment_info.max = max;
    segment_info.min = min;
    segment_info.std_dev = h_segment->GetStdDev();

    h_segment->Delete();
    return segment_info;
}

Reco_Particle Slice::get_pandora_p_by_ID(int id) {
  Reco_Particle pandora_p;
  for(int i_p = 0; i_p < this->pandora_particle.size(); i_p++) {
    if(this->pandora_particle.at(i_p).ID == id) pandora_p = this->pandora_particle.at(i_p);
  }
  return pandora_p;
}