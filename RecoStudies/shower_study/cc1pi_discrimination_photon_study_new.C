#include "../../Includes.h"
double get_min_distance(TVector3 P, TVector3 P_v, TVector3 v){
    TVector3 PQ = P - P_v;
    TVector3 crossProduct = PQ.Cross(v);
    return  crossProduct.Mag() / v.Mag();
}
G4Particle get_g4_mother(int g4_id, vector<G4Particle> g4_part_vec) {

  for (G4Particle g4_p: g4_part_vec) {
    if (g4_p.ID != g4_id) continue;
    for (G4Particle g4_p_mother: g4_part_vec) {
      if (g4_p_mother.ID != g4_p.mother) continue;
      return g4_p_mother;
    }
  }

  return G4Particle();

}

void cc1pi_discrimination_photon_study_new() {
  GenerateDictionaries();
  TTree *tree;
  TFile *input;

  Cut_Parameters cut_p;
  cut_p.min_track_score = 0.5;
  gStyle->SetPalette(56);

  //Declare the variables
  Slice *slice = 0;

  TH1* h_mother_common = new TH1D("h1", "h", 100, 0, 300);
  TH1* h_mother_not_common = new TH1D("h2", "h", 100, 0, 300);
  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_82p7k.root";
  input = new TFile(strRuta.c_str());
  tree = (TTree *) input->Get("tree");

  tree->SetBranchAddress("slice", &slice);

  vector<struct LocalBinInformation> bin_info = {
    {"min_distance_to_vtx",          "direction distance to vtx [cm]",           50, 0, 100},
    {"min_distance_to_vtx_true",     "direction distance to vtx true [cm]",      50, 0, 100},
    {"min_distance_to_vtx_all_true", "true direction distance to vtx true [cm]", 50, 0, 100},
    {"num_pho", "# #gamma", 5, 0, 5},
    //{"num_showers", "# Showers", 10, 0, 10}
  };
  const int num_hist = bin_info.size();

  vector<string> classification_names = {"CC1Pi_promt", "CC1Pi_non_prompt", "old_CC1Pi_mother_non_primary", "old_CC1Pi_mother_primary"};

  map<string, vector<TH1 *>> TH1_vec_map;
  for (int i_v = 0; i_v < classification_names.size(); i_v++) {
    vector<TH1 *> dummy_vec;
    for (int i_hist = 0; i_hist < num_hist; i_hist++) {
      string title = classification_names.at(i_v) + to_string(i_hist);
      dummy_vec.push_back(new TH1D(title.c_str(), " ", bin_info[i_hist].nBins, bin_info[i_hist].LowBin,
                                   bin_info[i_hist].UpBin));
    }
    TH1_vec_map[classification_names.at(i_v)] = dummy_vec;
  }



  map<string, double> Data;
  int nEntries = tree->GetEntries();
  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;



    int num_mu = slice->true_interaction.get_generator_num_primary_p(13);
    int num_pi = slice->true_interaction.get_generator_num_primary_p(211);
    int num_pi_0 = slice->true_interaction.get_generator_num_primary_p(111);
    bool is_CC1Pi = ((num_pi == 1) && (num_mu == 1));



    if (!is_CC1Pi) continue;
    vector<Reco_Particle> matched_pho;
    if(!slice->is_reconstructed) continue;
    if(slice->is_clear_cosmic) continue;

    G4Particle primary_muon = slice->true_interaction.get_g4_primary_p(13).at(0);
    G4Particle primary_pion = slice->true_interaction.get_g4_primary_p(211).at(0);

    if(primary_muon.TL < 3) continue;
    if(primary_pion.TL < 3) continue;


    for (Reco_Particle reco_p: slice->pandora_particle) {
      if (reco_p.track_score > 0.5) continue;
      if (abs(reco_p.matched_pdg) != 22) continue;
      matched_pho.push_back(reco_p);
    }
    Data["num_pho"]  = matched_pho.size();

    if (matched_pho.size() == 1) {
      for (Reco_Particle reco_p: matched_pho) {
        Data["min_distance_to_vtx"] = get_min_distance(slice->primary_vertex_reco.vertex_cm, reco_p.track_start, reco_p.shower_direction);
        Data["min_distance_to_vtx_true"]  = get_min_distance(slice->primary_vertex_true.vertex_cm, reco_p.track_start, reco_p.shower_direction);

        G4Particle mother = get_g4_mother(reco_p.true_track_id, slice->true_interaction.g4_particles);
        int mother_ID = mother.ID;
        TVector3 mother_dir = mother.P0;
        TVector3 mother_xf = mother.Xf;
        bool is_mother_primary = (mother.mother == 0);
        double mother_start_distance_to_vtx =  (mother.X0 - slice->primary_vertex_true.vertex_cm).Mag();

        TVector3 PQ_all_true = slice->primary_vertex_true.vertex_cm - mother_xf;
        TVector3 crossProduct_all_true = PQ_all_true.Cross(mother_dir);
        Data["min_distance_to_vtx_all_true"] = get_min_distance(slice->primary_vertex_true.vertex_cm, mother_xf, mother_dir);

        if(slice->true_interaction.is_selected_final_state("CC1Pi",0.325)) {
          for (int i_hist = 0; i_hist < num_hist; i_hist++) {
            if(mother_start_distance_to_vtx < 2) {
              TH1_vec_map["CC1Pi_promt"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
            } else {
              TH1_vec_map["CC1Pi_non_prompt"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
            }
          }
          //h_CC1Pi->Fill(reco_p.track_score, reco_p.track_kinetic_energy);
        } else {
          if(!slice->true_interaction.is_selected_final_state("old_CC1Pi",0.325)) continue;
          if(slice->true_interaction.is_selected_final_state("CC1Pi",0.325)) continue;
            for (int i_hist = 0; i_hist < num_hist; i_hist++) {
              if(is_mother_primary) {
                TH1_vec_map["old_CC1Pi_mother_primary"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
              } else {
                TH1_vec_map["old_CC1Pi_mother_non_primary"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
              }
            }
        }
      }
    } else if(matched_pho.size() > 0) {
      vector<int> used_mother_id;
      //True fill and #num_pho
      for (Reco_Particle reco_p: matched_pho) {
        //cout << parent_pdg << " " << reco_p.matched_pdg << endl;
        G4Particle mother = get_g4_mother(reco_p.true_track_id, slice->true_interaction.g4_particles);
        int mother_ID = mother.ID;
        TVector3 mother_dir = mother.P0;
        TVector3 mother_xf = mother.Xf;
        bool is_mother_primary = (mother.mother == 0);
        double mother_start_distance_to_vtx =  (mother.X0 - slice->primary_vertex_true.vertex_cm).Mag();

        if (std::find(used_mother_id.begin(), used_mother_id.end(), mother_ID) == used_mother_id.end()) {
          used_mother_id.push_back(mother_ID);
          //cout << mother_ID << " " << mother_PDG << endl;
          Data["min_distance_to_vtx_all_true"] = get_min_distance(slice->primary_vertex_true.vertex_cm, mother_xf, mother_dir);

          if(slice->true_interaction.is_selected_final_state("CC1Pi",0.325)) {
            for (int i_hist = 0; i_hist < num_hist; i_hist++) {
              if((bin_info[i_hist].FillDataType == "min_distance_to_vtx_all_true")||(bin_info[i_hist].FillDataType == "num_pho")) {
                if(mother_start_distance_to_vtx < 2) {
                  TH1_vec_map["CC1Pi_promt"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
                } else {
                  TH1_vec_map["CC1Pi_non_prompt"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
                }
              }
            }
            //h_CC1Pi->Fill(reco_p.track_score, reco_p.track_kinetic_energy);
          } else {
            if(!slice->true_interaction.is_selected_final_state("old_CC1Pi",0.325)) continue;
            if(slice->true_interaction.is_selected_final_state("CC1Pi",0.325)) continue;
            for (int i_hist = 0; i_hist < num_hist; i_hist++) {
              if((bin_info[i_hist].FillDataType == "min_distance_to_vtx_all_true")||(bin_info[i_hist].FillDataType == "num_pho")) {
                if (is_mother_primary) {
                  TH1_vec_map["old_CC1Pi_mother_primary"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
                } else {
                  TH1_vec_map["old_CC1Pi_mother_non_primary"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
                }
              }
            }
          }
        }
      } //END TRUE FILL


      //START RECO FILL

      //get posible phopair
      cout << slice->event_ID << " " << slice->subrun_ID << " " << slice->run_ID << " " << i_e << " " << slice->slice_ID <<  endl;
      cout << "crumbs score: " << slice->crumbs_score << endl;
      cout << "Pfparicles: " << endl;
      for (Reco_Particle reco_p: slice->pandora_particle) {
        cout << reco_p.ID << " " << reco_p.matched_pdg <<  endl;
      }
      cout << "Photons: "<< endl;
      for (Reco_Particle reco_p: matched_pho) {
        cout << reco_p.ID << " " << reco_p.matched_pdg <<  endl;
      }

      cout << "Pairs: "<< endl;
      std::vector<pair<Reco_Particle, Reco_Particle>> photon_pair_vec;
      for (Reco_Particle reco_p: matched_pho) {
        for (Reco_Particle reco_p_2: matched_pho) {
          if(reco_p.ID == reco_p_2.ID) continue;
          if(reco_p.ID < reco_p_2.ID) continue;
          pair<Reco_Particle, Reco_Particle> p_pair = make_pair(reco_p,reco_p_2);
          photon_pair_vec.push_back(p_pair);
        }
      }

      for(pair<Reco_Particle, Reco_Particle> photon_pair: photon_pair_vec) {
        cout << photon_pair.first.ID << " " << photon_pair.second.ID << endl;
      }

      /*
      //Get only the best pairs
      vector<int> used_pfp_ids;
      for (Reco_Particle reco_p: matched_pho) {
        double best_mass = 0;
        pair<Reco_Particle, Reco_Particle> best_pair;
        for (pair<Reco_Particle, Reco_Particle> photon_pair: photon_pair_vec) {
          if((photon_pair.first.ID != reco_p.ID) && (photon_pair.second.ID != reco_p.ID)) continue;
          if(find(used_pfp_ids.begin(), used_pfp_ids.end(), photon_pair.first.ID) != used_pfp_ids.end() ) continue;
          if(find(used_pfp_ids.begin(), used_pfp_ids.end(), photon_pair.second.ID) != used_pfp_ids.end() ) continue;

          TVector3 pho_1_P = photon_pair.first.shower_direction * photon_pair.first.shower_energy;
          TVector3 pho_2_P =  photon_pair.second.shower_direction * photon_pair.second.shower_energy;
          double inv_mass = sqrt(2*(pho_1_P.Mag()*pho_2_P.Mag() - pho_1_P*pho_2_P));
          if(abs(inv_mass - 134.9768) < abs(best_mass - 134.9768)) {
            best_mass = inv_mass;
            best_pair = photon_pair;
          }
        }
      }
       */

        for(pair<Reco_Particle, Reco_Particle> photon_pair: photon_pair_vec) {
        TVector3 pho_1_P = photon_pair.first.shower_direction * photon_pair.first.shower_energy;
        TVector3 pho_2_P =  photon_pair.second.shower_direction * photon_pair.second.shower_energy;
        double inv_mass = sqrt(2*(pho_1_P.Mag()*pho_2_P.Mag() - pho_1_P*pho_2_P));
        cout << photon_pair.first.ID << " " << photon_pair.second.ID << " invariant mass: " << inv_mass;

        G4Particle mother_1 = get_g4_mother(photon_pair.first.true_track_id, slice->true_interaction.g4_particles);
        G4Particle mother_2 = get_g4_mother(photon_pair.second.true_track_id, slice->true_interaction.g4_particles);




        if((mother_1.ID == mother_2.ID) && (mother_1.PDG == 111)) {
          cout << " common mother" << endl;
          h_mother_common->Fill(inv_mass);
        } else {
          cout << endl;
          h_mother_not_common->Fill(inv_mass);
        }

        bool is_mother_primary = (mother_1.mother == 0);
        double mother_start_distance_to_vtx =  (mother_1.X0 - slice->primary_vertex_true.vertex_cm).Mag();
        if(abs(inv_mass - 134.9778) < 50) {
          TVector3 position = photon_pair.first.track_start + photon_pair.second.track_start;
          position.SetX(position.X()/2);
          position.SetY(position.Y()/2);
          position.SetZ(position.Z()/2);
          Data["min_distance_to_vtx"] = get_min_distance(slice->primary_vertex_reco.vertex_cm, position, pho_1_P+pho_2_P);
          Data["min_distance_to_vtx_true"]  = get_min_distance(slice->primary_vertex_true.vertex_cm, position, pho_1_P+pho_2_P);

          if(slice->true_interaction.is_selected_final_state("CC1Pi",0.325)) {
            for (int i_hist = 0; i_hist < num_hist; i_hist++) {
              if((bin_info[i_hist].FillDataType == "min_distance_to_vtx_true")||(bin_info[i_hist].FillDataType == "min_distance_to_vtx")) {
                if(mother_start_distance_to_vtx < 2) {
                  TH1_vec_map["CC1Pi_promt"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
                } else {
                  TH1_vec_map["CC1Pi_non_prompt"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
                }
              }
            }
            //h_CC1Pi->Fill(reco_p.track_score, reco_p.track_kinetic_energy);
          } else {
            if(!slice->true_interaction.is_selected_final_state("old_CC1Pi",0.325)) continue;
            if(slice->true_interaction.is_selected_final_state("CC1Pi",0.325)) continue;
            for (int i_hist = 0; i_hist < num_hist; i_hist++) {
              if((bin_info[i_hist].FillDataType == "min_distance_to_vtx_true")||(bin_info[i_hist].FillDataType == "min_distance_to_vtx")) {
                if (is_mother_primary) {
                  TH1_vec_map["old_CC1Pi_mother_primary"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
                } else {
                  TH1_vec_map["old_CC1Pi_mother_non_primary"][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
                }
              }
            }
          }

        }
      }


    }
  }









  TLegend *leg = new TLegend(0.55, 0.6, .87, .87);

  //leg->AddEntry(h_old_CC1Pi_pi0_mother[0], "Old CC1pi distribution with pi0", "l");
  //leg->AddEntry(h_old_CC1Pi_others[0], "Old CC1pi distribution without pi0", "l");
  //leg->AddEntry(h_old_CC1Pi_pi0_mother_primary[0], "Old CC1pi distribution with pi0 mother primary", "l");


  string folder = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/shower_studies/pho_dir_study/";
  for (int i_v = 0; i_v < classification_names.size(); i_v++) {
    leg->AddEntry(TH1_vec_map[classification_names.at(i_v)][0], classification_names.at(i_v).c_str(), "l");
  }
  for (int i_hist = 0; i_hist < num_hist; i_hist++) {
    TCanvas *ci = new TCanvas();
    double Max = 0;
    for (int i_v = 0; i_v < classification_names.size(); i_v++) {
      TH1_vec_map[classification_names.at(i_v)][i_hist]->Scale(1. / TH1_vec_map[classification_names.at(i_v)][i_hist]->Integral());
      if(TH1_vec_map[classification_names.at(i_v)][i_hist]->GetMaximum() > Max) Max = TH1_vec_map[classification_names.at(i_v)][i_hist]->GetMaximum();
      TH1_vec_map[classification_names.at(0)][i_hist]->GetYaxis()->SetRangeUser(0,Max*1.2);
      TH1_vec_map[classification_names.at(i_v)][i_hist]->SetStats(0);

      if(i_v ==  0) {
        TH1_vec_map[classification_names.at(i_v)][i_hist]->Draw("hist");
      } else {
        TH1_vec_map[classification_names.at(i_v)][i_hist]->Draw("hist same");
      }
      TH1_vec_map[classification_names.at(i_v)][i_hist]->SetTitle((";" + bin_info[i_hist].Title + "; Events").c_str());
      TH1_vec_map[classification_names.at(i_v)][i_hist]->SetLineColor(colors.at(i_v));
    }
    leg->Draw();
    ci->SaveAs((folder + "/normalzed_mo_shot_muon_like_"+ bin_info.at(i_hist).FillDataType +  ".pdf").c_str());
    ci->Close();
    delete ci;
  }

  TCanvas* c3 = new TCanvas();
  h_mother_common->Draw("hist");
  h_mother_not_common->SetLineColor(kRed + 2);
  h_mother_not_common->Draw("hist same");

}