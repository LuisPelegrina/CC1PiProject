#include "../../Includes.h"
#include "Graphs_utils.cpp"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}

void Energy_comp()
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
  GenerateDictionaries();
  TTree *tree;
  TFile *input;
    
  Cut_Parameters cut_p;
  cut_p.min_track_lenght = 0;
  cut_p.min_track_score = 0;


  //Declare the variables
  Slice *slice = 0;

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
  input = new TFile(strRuta.c_str());
  tree =(TTree*)input->Get("tree");
  tree->SetBranchAddress("slice", &slice);
  int nEntries = tree->GetEntries();

  std::vector<TGraph*> gr_general_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_muon_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_pion_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_proton_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_other_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_muon_primary_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_pion_primary_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_proton_primary_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_other_primary_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_primary_vertex = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_primary_vertex_truth = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_secondary_vertex = {new TGraph(), new TGraph(), new TGraph()};


  TGraph2D* gr_general_sp_2D = new TGraph2D();
  TGraph2D* gr_muon_track_sp_2D = new TGraph2D();
  TGraph2D* gr_pion_track_sp_2D = new TGraph2D();
  TGraph2D* gr_proton_track_sp_2D = new TGraph2D();
  TGraph2D* gr_other_track_sp_2D = new TGraph2D();
  TGraph2D* gr_muon_primary_track_sp_2D = new TGraph2D();
  TGraph2D* gr_pion_primary_track_sp_2D = new TGraph2D();
  TGraph2D* gr_proton_primary_track_sp_2D = new TGraph2D();
  TGraph2D* gr_other_primary_track_sp_2D = new TGraph2D();
  TGraph2D* gr_primary_vertex_2D = new TGraph2D();
  TGraph2D* gr_primary_vertex_truth_2D = new TGraph2D();
  TGraph2D* gr_secondary_vertex_2D = new TGraph2D();

  std::vector<TCanvas*> c = {new TCanvas("c1", "My Canvas", 500, 400), new TCanvas("c2", "My Canvas", 500, 400), new TCanvas("c3", "My Canvas", 500, 400)};
  c[0]->SetWindowPosition(10, 50);
  c[1]->SetWindowPosition(700, 50);
  c[2]->SetWindowPosition(10, 500);

  TCanvas* c2D = new TCanvas("c4", "My Canvas", 500, 400);
  c2D->SetWindowPosition(700, 500);

  TCanvas* button = new TCanvas("c12", "My Canvas", 200, 200);
  button->SetWindowPosition(1400, 50);
  button->SetFillColor(kRed);
  button->Draw();

  cut_p.min_distance_to_CPA = 20;
  cut_p.min_distance_to_first_z_wall = 10;
  cut_p.min_distance_to_last_z_wall = 25;
  cut_p.min_distance_to_wall_x_y = 20;

  const int num_cuts = 4;
  string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut"};

  //pp = Primary Particle
  map<string,double> cases_cont;

  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

    if(slice->slice_ID != 0) continue;
    if(!slice->true_interaction.is_selected_final_state("N_CC1Pi")) continue;

    bool pass_cut = true;
    //Plot only if correctly reconstructed
    for(int i_cut = 0; i_cut < num_cuts; i_cut++) {
      if(!pass_cut) continue;
      if(!slice->pass_analysis_cut(cut_p, cut_name[i_cut])) {
        pass_cut = false;
      }
    }
    if(!pass_cut) continue;

    G4Particle primary_muon = slice->true_interaction.get_g4_primary_p(13).at(0);
    G4Particle primary_pion = slice->true_interaction.get_g4_primary_p(211).at(0);

    int cont_primary_p = 0;
    int cont_p = slice->pandora_particle.size();
    map<int, int> cont_by_pdg_pur_80;
    map<int, int> cont_by_pdg_pur_50;
    map<int, int> cont_by_pdg;

    for(int i_p = 0; i_p < slice->pandora_particle.size() ;i_p++) {
      Reco_Particle pandora_p = slice->pandora_particle.at(i_p);
        if(slice->pandora_particle.at(i_p).is_pandora_primary){
          cont_by_pdg[11111]++;
          cont_by_pdg[abs(pandora_p.matched_pdg)]++;
          if((pandora_p.completeness > 0.8) && (pandora_p.purity > 0.8)) cont_by_pdg_pur_80[abs(pandora_p.matched_pdg)]++;
          if((pandora_p.completeness > 0.5) && (pandora_p.purity > 0.5)) cont_by_pdg_pur_50[abs(pandora_p.matched_pdg)]++;
          if((abs(pandora_p.matched_pdg) == 211) || (abs(pandora_p.matched_pdg) == 13)) continue;
          cont_by_pdg[abs(pandora_p.matched_pdg)]++;
          cont_by_pdg[-1]++;
        }
    }

    int num_extra = cont_by_pdg[11111] - 2;

    double distance_to_true_vertex = (slice->primary_vertex_true.vertex_cm - slice->primary_vertex_reco.vertex_cm).Mag();
    bool show = false;
    bool classified = false;

    cases_cont["tot"] += 1;
    if(num_extra == 0) {
      if((cont_by_pdg_pur_80[211] == 1) && (cont_by_pdg_pur_80[13] == 1)) {
        cases_cont["good_reco_80"]++;
        show = false;
      } else if ((cont_by_pdg_pur_50[211] == 1) && (cont_by_pdg_pur_50[13] == 1)){
        cases_cont["good_reco_50"]++;
        show = false;
      } else if((cont_by_pdg[211] == 1) && (cont_by_pdg[13] == 1)) {
        cases_cont["miss_reco"] += 1;
        cases_cont["miss_reco_mu_pi_purlow"] += 1;
        show = false;
      } else {
        cases_cont["miss_reco"] += 1;
        string miss_reco_name = "miss_reco_2_pfp_missing_mu_pi";
        cases_cont[miss_reco_name] += 1;

        if(cont_by_pdg[13] > 0) {
          miss_reco_name += "_only_mu";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        } else if(cont_by_pdg[211] > 0) {
          miss_reco_name += "_only_pi";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }

        } else {
          miss_reco_name += "_no_mu_pi";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        }

      }
    } else if (num_extra > 0) {
      cases_cont["miss_reco"] += 1;
      string miss_reco_name = "miss_reco";

      if ((cont_by_pdg_pur_80[211] == 1) && (cont_by_pdg_pur_80[13] == 1)) {
        miss_reco_name += "_mu_pi_pur80_others";
        cases_cont[miss_reco_name]++;
        double distance_to_closest_non_80 = 1212212;
        int pandora_extra_p_id = 0;
        for (int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
          Reco_Particle pandora_p = slice->pandora_particle.at(i_p);
          if (slice->pandora_particle.at(i_p).is_pandora_primary) {
            if(( (pandora_p.matched_pdg == 211) || (pandora_p.matched_pdg = 13)) &&(pandora_p.completeness > 0.8) && (pandora_p.purity > 0.8)) continue;
            pandora_extra_p_id = i_p;
            if(distance_to_closest_non_80 > (pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag()) distance_to_closest_non_80 = (pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag();
          }
        }
        if(num_extra == 1)  {
          miss_reco_name += "_1";
          cases_cont[miss_reco_name]++;

          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;


            if(distance_to_closest_non_80 > 8) {
              miss_reco_name += "_p_far_from_vertex";
              cases_cont[miss_reco_name]++;
              show = false;
            } else {
              miss_reco_name += "_p_close_to_vertex";
              show = false;
              cases_cont[miss_reco_name]++;
              if((slice->true_interaction.get_g4_primary_p(211).at(0).TL < 3) && ((slice->true_interaction.get_g4_primary_p(211).at(0).end_process.compare("pi+Inelastic") == 0) || (slice->true_interaction.get_g4_primary_p(211).at(0).end_process.compare("pi-Inelastic") == 0))) {
                miss_reco_name += "_promt_pion_interaction";
                show = false;
                cases_cont[miss_reco_name]++;
              } else if((slice->pandora_particle.at(pandora_extra_p_id).completeness < 0.1)|| (slice->pandora_particle.at(pandora_extra_p_id).track_lenght < 3)) {
                miss_reco_name += "_really_low_completenesss_or_TL";
                show = false;
                cases_cont[miss_reco_name]++;

              } else {
                miss_reco_name += "_other";
                show = false;
                cases_cont[miss_reco_name]++;
              }

            }
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }

        } else if(num_extra > 1)  {
          miss_reco_name += "_2+";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;

            if(distance_to_closest_non_80 > 8) {
              miss_reco_name += "_p_far_from_vertex";
              cases_cont[miss_reco_name]++;
              show = false;
            } else {
              miss_reco_name += "_p_close_to_vertex";
              show = false;
              cases_cont[miss_reco_name]++;
            }
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        }


      } else if ((cont_by_pdg_pur_50[211] == 1) && (cont_by_pdg_pur_50[13] == 1)) {
        miss_reco_name += "_mu_pi_pur50_others";
        cases_cont[miss_reco_name]++;
        show = false;

        double distance_to_closest_non_50 = 1212212;
        int pandora_extra_p_id = 0;
        for (int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
          Reco_Particle pandora_p = slice->pandora_particle.at(i_p);
          if (slice->pandora_particle.at(i_p).is_pandora_primary) {
            if(( (pandora_p.matched_pdg == 211) || (pandora_p.matched_pdg = 13)) &&(pandora_p.completeness > 0.5) && (pandora_p.purity > 0.5)) continue;
            pandora_extra_p_id = i_p;
            if(distance_to_closest_non_50 > (pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag()) distance_to_closest_non_50 = (pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag();
          }
        }
        if(num_extra == 1)  {
          miss_reco_name += "_1";
          cases_cont[miss_reco_name]++;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
            if(distance_to_closest_non_50 > 8) {
              miss_reco_name += "_p_far_from_vertex";
              cases_cont[miss_reco_name]++;
              show = false;
            } else {
              miss_reco_name += "_p_close_to_vertex";
              show = false;
              cases_cont[miss_reco_name]++;
            }
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        } else if(num_extra > 1)  {
          miss_reco_name += "_2+";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
            if(distance_to_closest_non_50 > 8) {
              miss_reco_name += "_p_far_from_vertex";
              cases_cont[miss_reco_name]++;
              show = false;
            } else {
              miss_reco_name += "_p_close_to_vertex";
              show = false;
              cases_cont[miss_reco_name]++;
            }
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        }


      } else if ((cont_by_pdg[211] == 1) && (cont_by_pdg[13] == 1)){
        miss_reco_name += "_mu_pi_purlow_others";
        cases_cont[miss_reco_name]++;
        show = false;
        if(num_extra == 1)  {
          miss_reco_name += "_1";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        } else if(num_extra > 1)  {
          miss_reco_name += "_2+";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        }
      } else {
        //CÓDIGO PARA CLASIFICAR
        miss_reco_name += "_2+_pfp_missing_mu_pi";
        cases_cont[miss_reco_name]++;
        if(cont_by_pdg[13] > 0) {
          miss_reco_name += "_only_mu";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
            if(slice->true_interaction.get_g4_primary_p(211).at(0).TL < 3) {
              miss_reco_name += "_low_tl_pion";
              cases_cont[miss_reco_name]++;
            } else {
              show = false;
            }
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            if(slice->true_interaction.get_g4_primary_p(211).at(0).TL < 3) {
              miss_reco_name += "_low_tl_pion";
              cases_cont[miss_reco_name]++;
            } else {
              show = false;
            }
          }

        } else if(cont_by_pdg[211] > 0) {
          miss_reco_name += "_only_pi";
          cases_cont[miss_reco_name]++;
          show = false;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        } else {
          miss_reco_name += "_no_mu_pi";
          cases_cont[miss_reco_name]++;
          if(distance_to_true_vertex < 2) {
            miss_reco_name += "_vertex_ok";
            cases_cont[miss_reco_name]++;
            show = false;
          } else {
            miss_reco_name += "_miss_vertex";
            cases_cont[miss_reco_name]++;
            show = false;
          }
        }
      }
    } if(num_extra < 0) {
      cases_cont["miss_reco"] += 1;
      string miss_reco_name = "miss_reco_<1_pfp_missing_mu_pi";
      cases_cont[miss_reco_name]++;
      if(num_extra == -2) {
        miss_reco_name += "_no_p";
        cases_cont[miss_reco_name]++;
        show = false;
      } else {
        if(cont_by_pdg[13] > 0) {
          miss_reco_name += "_only_mu";
          cases_cont[miss_reco_name]++;
          show = false;

        } else if(cont_by_pdg[211] > 0) {
          miss_reco_name += "_only_pi";
          cases_cont[miss_reco_name]++;
          show = false;

        } else {
          miss_reco_name += "_no_mu_pi";
          cases_cont[miss_reco_name]++;
          show = false;
        }
      }

    }

    if(!show) continue;
    slice->print(cut_p, true, true, true, true,false, true);
    cout << "DISTANCE TO TRUE VERTEX" << distance_to_true_vertex << endl;
    cout << "PRINT ENDED" << endl;


    int TPC = 0;
    if(slice->primary_vertex_true.vertex_cm.X() > 0) TPC = 1;
    cout << "TPC Chosen" << endl;

    //bucle para los hits
    for(int i_h = 0; i_h < slice->hits.size(); i_h++) {
      if (slice->hits.at(i_h).TPC_ID != TPC) continue;
      Hit hit = slice->hits.at(i_h);
     // if(slice->get_distance_to_vertex(hit) > 200 ) continue;
      gr_general_hits[hit.Plane_ID]->SetPoint(gr_general_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
      if(hit.associated_pfp_ID == -1) continue;

      //Search for the pfp that the hit is associated to:
      int particle_index = 0;
      for(int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
        if(slice->pandora_particle.at(i_p).ID == hit.associated_pfp_ID) particle_index = i_p;
      }
      if(slice->pandora_particle.at( particle_index).is_pandora_primary) {
        if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 211) {
          gr_pion_primary_track_hits[hit.Plane_ID]->SetPoint(gr_pion_primary_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 13) {
          gr_muon_primary_track_hits[hit.Plane_ID]->SetPoint(gr_muon_primary_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 2212) {
          gr_proton_primary_track_hits[hit.Plane_ID]->SetPoint(gr_proton_primary_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else {
          gr_other_primary_track_hits[hit.Plane_ID]->SetPoint(gr_other_primary_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        }
      } else {
        if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 211) {
          gr_pion_track_hits[hit.Plane_ID]->SetPoint(gr_pion_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 13) {
          gr_muon_track_hits[hit.Plane_ID]->SetPoint(gr_muon_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 2212) {
          gr_proton_track_hits[hit.Plane_ID]->SetPoint(gr_proton_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else {
          gr_other_track_hits[hit.Plane_ID]->SetPoint(gr_other_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        }
      }

    }
    cout << "hits done" << endl;

    //bucle en las 3 vistas para el vértice primario
    for(int i_p = 0; i_p < 3; i_p++) {
      double vertex_drift_t = slice->primary_vertex_reco.vertex_wire.at(i_p).drift_t;
      int vertex_wire_ID = slice->primary_vertex_reco.vertex_wire.at(i_p).Wire_ID;

      gr_primary_vertex[i_p]->SetPoint(gr_primary_vertex[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
      gr_primary_vertex[i_p]->SetPoint(gr_primary_vertex[i_p]->GetN(), vertex_wire_ID, vertex_drift_t);
      gr_primary_vertex[i_p]->SetPoint(gr_primary_vertex[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
    }
    cout << "Primary vertex done" << endl;

    //bucle en las 3 vistas para el vértice true
    for(int i_p = 0; i_p < 3; i_p++) {
      double vertex_drift_t = slice->primary_vertex_true.vertex_wire.at(i_p).drift_t;
      int vertex_wire_ID = slice->primary_vertex_true.vertex_wire.at(i_p).Wire_ID;
      gr_primary_vertex_truth[i_p]->SetPoint(gr_primary_vertex_truth[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
      gr_primary_vertex_truth[i_p]->SetPoint(gr_primary_vertex_truth[i_p]->GetN(), vertex_wire_ID, vertex_drift_t);
      gr_primary_vertex_truth[i_p]->SetPoint(gr_primary_vertex_truth[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
    }
    cout << "true vertex done" << endl;

    //bucle en las 3 vistas para vertices secundarios true
    for(int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      for(int i_p = 0; i_p < 3; i_p++) {
        double vertex_drift_t = slice->vertex_vec.at(i_v).vertex_wire.at(i_p).drift_t;
        int vertex_wire_ID = slice->vertex_vec.at(i_v).vertex_wire.at(i_p).Wire_ID;
        gr_secondary_vertex[i_p]->SetPoint(gr_secondary_vertex[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
        gr_secondary_vertex[i_p]->SetPoint(gr_secondary_vertex[i_p]->GetN(), vertex_wire_ID, vertex_drift_t);
        gr_secondary_vertex[i_p]->SetPoint(gr_secondary_vertex[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
      }
    }
    cout << "secondary vertex done" << endl;

    //Dibuja las 3 vistas
    for(int i_p = 0; i_p < 3; i_p++) {
      c[i_p]->cd();

      gr_general_hits[i_p]->SetTitle(";Wire ID; Time [#mu s]");
      gr_general_hits[i_p]->SetMarkerSize(0.5);
      gr_general_hits[i_p]->SetMarkerColor(kGray+2);

      gr_muon_track_hits[i_p]->SetMarkerSize(0.5);
      gr_muon_track_hits[i_p]->SetMarkerColor(kViolet+2);
      gr_muon_track_hits[i_p]->SetMarkerStyle(8);

      gr_pion_track_hits[i_p]->SetMarkerSize(0.5);
      gr_pion_track_hits[i_p]->SetMarkerColor(kViolet+2);
      gr_pion_track_hits[i_p]->SetMarkerStyle(8);

      gr_proton_track_hits[i_p]->SetMarkerSize(0.5);
      gr_proton_track_hits[i_p]->SetMarkerColor(kViolet+2);
      gr_proton_track_hits[i_p]->SetMarkerStyle(8);

      gr_other_track_hits[i_p]->SetMarkerSize(0.5);
      gr_other_track_hits[i_p]->SetMarkerColor(kViolet+2);
      gr_other_track_hits[i_p]->SetMarkerStyle(8);


      gr_muon_primary_track_hits[i_p]->SetMarkerSize(0.5);
      gr_muon_primary_track_hits[i_p]->SetMarkerColor(kGreen+2);
      gr_muon_primary_track_hits[i_p]->SetMarkerStyle(8);

      gr_pion_primary_track_hits[i_p]->SetMarkerSize(0.5);
      gr_pion_primary_track_hits[i_p]->SetMarkerColor(kRed+2);
      gr_pion_primary_track_hits[i_p]->SetMarkerStyle(8);

      gr_proton_primary_track_hits[i_p]->SetMarkerSize(0.5);
      gr_proton_primary_track_hits[i_p]->SetMarkerColor(kBlue+2);
      gr_proton_primary_track_hits[i_p]->SetMarkerStyle(8);

      gr_other_primary_track_hits[i_p]->SetMarkerSize(0.5);
      gr_other_primary_track_hits[i_p]->SetMarkerColor(kOrange+2);
      gr_other_primary_track_hits[i_p]->SetMarkerStyle(8);

      gr_primary_vertex[i_p]->SetMarkerSize(1);
      gr_primary_vertex[i_p]->SetMarkerColor(kBlack);
      gr_primary_vertex[i_p]->SetMarkerStyle(8);

      gr_primary_vertex_truth[i_p]->SetMarkerSize(1);
      gr_primary_vertex_truth[i_p]->SetMarkerColor(kBlue +2 );
      gr_primary_vertex_truth[i_p]->SetMarkerStyle(22);

      gr_secondary_vertex[i_p]->SetMarkerSize(1);
      gr_secondary_vertex[i_p]->SetMarkerColor(kBlack);
      gr_secondary_vertex[i_p]->SetMarkerStyle(34);



      TMultiGraph *mg = new TMultiGraph();
      if(gr_general_hits[i_p]->GetN() > 0)mg->Add(gr_general_hits[i_p],"AP");
      if(gr_muon_track_hits[i_p]->GetN() > 0)mg->Add(gr_muon_track_hits[i_p],"AP");
      if(gr_pion_track_hits[i_p]->GetN() > 0)mg->Add(gr_pion_track_hits[i_p],"AP");
      if(gr_proton_track_hits[i_p]->GetN() > 0)mg->Add(gr_proton_track_hits[i_p],"AP");
      if(gr_other_track_hits[i_p]->GetN() > 0)mg->Add(gr_other_track_hits[i_p],"AP");
      if(gr_muon_primary_track_hits[i_p]->GetN() > 0)mg->Add(gr_muon_primary_track_hits[i_p],"AP");
      if(gr_pion_primary_track_hits[i_p]->GetN() > 0)mg->Add(gr_pion_primary_track_hits[i_p],"AP");
      if(gr_proton_primary_track_hits[i_p]->GetN() > 0)mg->Add(gr_proton_primary_track_hits[i_p],"AP");
      if(gr_other_primary_track_hits[i_p]->GetN() > 0)mg->Add(gr_other_primary_track_hits[i_p],"AP");
      if(gr_primary_vertex_truth[i_p]->GetN() > 0)mg->Add(gr_primary_vertex_truth[i_p],"AP");
      if(gr_primary_vertex[i_p]->GetN() > 0)mg->Add(gr_primary_vertex[i_p],"AP");
      if(gr_secondary_vertex[i_p]->GetN() > 0)mg->Add(gr_secondary_vertex[i_p],"AP");

      mg->SetTitle(";Wire ID; Time [#mu s]");
      mg->Draw("A");

      c[i_p]->Update();
    }

    //Lo mimso pero en 3D

    double max_sp_x=0;
    double max_sp_y=0;
    double max_sp_z=0;
    double min_sp_x=2000;
    double min_sp_y=2000;
    double min_sp_z=2000;


    for(int i_sp = 0; i_sp < slice->space_point_vec.size(); i_sp++) {
      SpacePoint sp = slice->space_point_vec.at(i_sp);
      gr_general_sp_2D->SetPoint(gr_general_sp_2D->GetN(), sp.z, sp.x, sp.y);
      if(sp.associated_pfp_ID == -1) continue;

      if(max_sp_x > sp.x) max_sp_x = sp.x;
      if(max_sp_y > sp.y) max_sp_y = sp.y;
      if(max_sp_z > sp.z) max_sp_z = sp.z;
      if(min_sp_x < sp.x) min_sp_x = sp.x;
      if(min_sp_y < sp.y) min_sp_y = sp.y;
      if(min_sp_z < sp.z) min_sp_z = sp.z;

      //Search for the pfp that the hit is associated to:
      int particle_index = 0;
      for(int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
        if(slice->pandora_particle.at(i_p).ID == sp.associated_pfp_ID) particle_index = i_p;
      }
      if(slice->pandora_particle.at( particle_index).is_pandora_primary) {
        if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 211) {
          gr_pion_primary_track_sp_2D->SetPoint(gr_pion_primary_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 13) {
          gr_muon_primary_track_sp_2D->SetPoint(gr_muon_primary_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 2212) {
          gr_proton_primary_track_sp_2D->SetPoint(gr_proton_primary_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else {
          gr_other_primary_track_sp_2D->SetPoint(gr_other_primary_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        }
      } else {
        if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 211) {
          gr_pion_track_sp_2D->SetPoint(gr_pion_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 13) {
          gr_muon_track_sp_2D->SetPoint(gr_muon_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 2212) {
          gr_proton_track_sp_2D->SetPoint(gr_proton_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else {
          gr_other_track_sp_2D->SetPoint(gr_other_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        }
      }

    }

    cout << "sp done" << endl;
    //bucle en las 3 vistas para el vértice primario
    gr_primary_vertex_2D->SetPoint(gr_primary_vertex_2D->GetN(),slice->primary_vertex_reco.vertex_cm.Z(),slice->primary_vertex_reco.vertex_cm.X(),slice->primary_vertex_reco.vertex_cm.Y());
    cout << "Primary vertex done" << endl;

    gr_primary_vertex_truth_2D->SetPoint(gr_primary_vertex_2D->GetN(),slice->primary_vertex_true.vertex_cm.Z(),slice->primary_vertex_true.vertex_cm.X(),slice->primary_vertex_true.vertex_cm.Y());
    cout << "True vertex done" << endl;

    //bucle en las 3 vistas para vertices secundarios true
    for(int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      gr_secondary_vertex_2D->SetPoint(gr_primary_vertex_2D->GetN(),slice->vertex_vec.at(i_v).vertex_cm.Z(),slice->vertex_vec.at(i_v).vertex_cm.X(),slice->vertex_vec.at(i_v).vertex_cm.Y());
    }
    cout << "secondary vertex done" << endl;

    c2D->cd();
    gr_general_sp_2D->SetTitle(";z [cm]; y [cm];x [cm]");
    gr_general_sp_2D->SetMarkerSize(0.5);
    gr_general_sp_2D->SetMarkerColor(kViolet+2);

    gr_muon_track_sp_2D->SetMarkerSize(0.5);
    gr_muon_track_sp_2D->SetMarkerColor(kViolet+2);
    gr_muon_track_sp_2D->SetMarkerStyle(8);

    gr_pion_track_sp_2D->SetMarkerSize(0.5);
    gr_pion_track_sp_2D->SetMarkerColor(kViolet+2);
    gr_pion_track_sp_2D->SetMarkerStyle(8);

    gr_proton_track_sp_2D->SetMarkerSize(0.5);
    gr_proton_track_sp_2D->SetMarkerColor(kViolet+2);
    gr_proton_track_sp_2D->SetMarkerStyle(8);

    gr_other_track_sp_2D->SetMarkerSize(0.5);
    gr_other_track_sp_2D->SetMarkerColor(kViolet+2);
    gr_other_track_sp_2D->SetMarkerStyle(8);

    gr_muon_primary_track_sp_2D->SetMarkerSize(0.5);
    gr_muon_primary_track_sp_2D->SetMarkerColor(kGreen+2);
    gr_muon_primary_track_sp_2D->SetMarkerStyle(8);

    gr_pion_primary_track_sp_2D->SetMarkerSize(0.5);
    gr_pion_primary_track_sp_2D->SetMarkerColor(kRed+2);
    gr_pion_primary_track_sp_2D->SetMarkerStyle(8);

    gr_proton_primary_track_sp_2D->SetMarkerSize(0.5);
    gr_proton_primary_track_sp_2D->SetMarkerColor(kBlue+2);
    gr_proton_primary_track_sp_2D->SetMarkerStyle(8);

    gr_other_primary_track_sp_2D->SetMarkerSize(0.5);
    gr_other_primary_track_sp_2D->SetMarkerColor(kOrange+2);
    gr_other_primary_track_sp_2D->SetMarkerStyle(8);

    gr_primary_vertex_2D->SetMarkerSize(1);
    gr_primary_vertex_2D->SetMarkerColor(kBlack);
    gr_primary_vertex_2D->SetMarkerStyle(8);

    gr_primary_vertex_truth_2D->SetMarkerSize(1);
    gr_primary_vertex_truth_2D->SetMarkerColor(kBlue +2 );
    gr_primary_vertex_truth_2D->SetMarkerStyle(22);

    gr_secondary_vertex_2D->SetMarkerSize(1);
    gr_secondary_vertex_2D->SetMarkerColor(kBlack);
    gr_secondary_vertex_2D->SetMarkerStyle(34);

    //SetPoints2D(gr_general_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_muon_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_pion_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_proton_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_other_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_muon_primary_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_pion_primary_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_proton_primary_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_other_primary_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_primary_vertex_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_primary_vertex_truth_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_secondary_vertex_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);

    if(gr_general_sp_2D->GetN() > 0) gr_general_sp_2D->Draw("P");
    if(gr_muon_track_sp_2D->GetN() > 0) gr_muon_track_sp_2D->Draw("SAME P");
    if(gr_pion_track_sp_2D->GetN() > 0) gr_pion_track_sp_2D->Draw("SAME P");
    if(gr_proton_track_sp_2D->GetN() > 0) gr_proton_track_sp_2D->Draw("SAME P");
    if(gr_other_track_sp_2D->GetN() > 0) gr_other_track_sp_2D->Draw("SAME P");
    if(gr_muon_primary_track_sp_2D->GetN() > 0) gr_muon_primary_track_sp_2D->Draw("SAME P");
    if(gr_pion_primary_track_sp_2D->GetN() > 0) gr_pion_primary_track_sp_2D->Draw("SAME P");
    if(gr_proton_primary_track_sp_2D->GetN() > 0) gr_proton_primary_track_sp_2D->Draw("SAME P");
    if(gr_other_primary_track_sp_2D->GetN() > 0) gr_other_primary_track_sp_2D->Draw("SAME P");
    if(gr_primary_vertex_2D->GetN() > 0) gr_primary_vertex_2D->Draw("SAME P");
    if(gr_primary_vertex_truth_2D->GetN() > 0) gr_primary_vertex_truth_2D->Draw("SAME P");
    if(gr_secondary_vertex_2D->GetN() > 0) gr_secondary_vertex_2D->Draw("SAME P");
    c2D->Update();

    button->WaitPrimitive();

    bool clicked = false;
    while(!clicked) {
      if(button->WaitPrimitive() == nullptr) {
        std::cout << "Canvas 1 clicked! Proceeding..." << std::endl;
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }
    for(int i_p = 0; i_p < 3; i_p++) {
      int N_points = gr_general_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_general_hits[i_p]->RemovePoint(0);
      N_points = gr_pion_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_pion_track_hits[i_p]->RemovePoint(0);
      N_points = gr_muon_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_muon_track_hits[i_p]->RemovePoint(0);
      N_points = gr_proton_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_proton_track_hits[i_p]->RemovePoint(0);
      N_points = gr_other_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_other_track_hits[i_p]->RemovePoint(0);

      N_points = gr_pion_primary_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_pion_primary_track_hits[i_p]->RemovePoint(0);
      N_points = gr_muon_primary_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_muon_primary_track_hits[i_p]->RemovePoint(0);
      N_points = gr_proton_primary_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_proton_primary_track_hits[i_p]->RemovePoint(0);
      N_points = gr_other_primary_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_other_primary_track_hits[i_p]->RemovePoint(0);


      N_points = gr_primary_vertex[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_primary_vertex[i_p]->RemovePoint(0);
      N_points = gr_primary_vertex_truth[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_primary_vertex_truth[i_p]->RemovePoint(0);
      N_points = gr_secondary_vertex[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_secondary_vertex[i_p]->RemovePoint(0);
    }


    int N_points = gr_general_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_general_sp_2D->RemovePoint(0);
    N_points = gr_pion_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_pion_track_sp_2D->RemovePoint(0);
    N_points = gr_muon_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_muon_track_sp_2D->RemovePoint(0);
    N_points = gr_proton_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_proton_track_sp_2D->RemovePoint(0);
    N_points = gr_other_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_other_track_sp_2D->RemovePoint(0);

    N_points = gr_pion_primary_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_pion_primary_track_sp_2D->RemovePoint(0);
    N_points = gr_muon_primary_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_muon_primary_track_sp_2D->RemovePoint(0);
    N_points =  gr_proton_primary_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_proton_primary_track_sp_2D->RemovePoint(0);
    N_points = gr_other_primary_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_other_primary_track_sp_2D->RemovePoint(0);

    N_points = gr_primary_vertex_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_primary_vertex_2D->RemovePoint(0);
    N_points = gr_primary_vertex_truth_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_primary_vertex_truth_2D->RemovePoint(0);
    N_points = gr_secondary_vertex_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_secondary_vertex_2D->RemovePoint(0);
    /*
    if(cont_primary_p < 2) {
        miss_reco_name += "nopp";
        cases_cont[miss_reco_name]++;
        if(cont_p < 2) {
            miss_reco_name += "_nop_" + to_string(cont_primary_pion) +"pi_" + to_string(cont_primary_muon) + "mu";
            cases_cont[miss_reco_name]++;
        } else {
            miss_reco_name += "_hasp_" + to_string(cont_primary_pion) +"pi_" + to_string(cont_primary_muon) + "mu";
            cases_cont[miss_reco_name]++;
        }
        classified = true;
    } else {
        miss_reco_name += "haspp";
        */

       /*
        for (std::map<int, int>::iterator it = cont_by_pdg.begin(); it != cont_by_pdg.end(); ++it) {
            if((it->first == 211) || (it->first == 13)) continue;
            miss_reco_name += to_string(it->second) + to_string(it->first) +"_";
        }
        */

        /*
        int closet_PDG = 0;
        string start_or_end = "start";
        string prim_or_sec = "primary";
        double closest_distance = (slice->primary_vertex_true.vertex_cm - slice->primary_vertex_reco.vertex_cm).Mag();
        for(int i_p = 0; i_p < slice->true_interaction.g4_particles.size(); i_p++) {
            G4Particle g4_part = slice->true_interaction.g4_particles.at(i_p);

            double start_distance_to_vertex = (g4_part.X0 - slice->primary_vertex_reco.vertex_cm).Mag();
            double end_distance_to_vertex = (g4_part.XL - slice->primary_vertex_reco.vertex_cm).Mag();
            if(!g4_part.is_primary()) continue;
            if(end_distance_to_vertex < closest_distance - 0.01) {
                closest_distance = end_distance_to_vertex;
                closet_PDG = g4_part.PDG;
                start_or_end = "end";
                prim_or_sec = "primary";
            }

        }

        if ((slice->primary_vertex_true.vertex_cm - slice->primary_vertex_reco.vertex_cm).Mag() > 3){
           if(closest_distance < 3) {
                miss_reco_name += "vertex_in_" + start_or_end + "_of_" + prim_or_sec + "_" + to_string(abs(closet_PDG));
            } else {
                miss_reco_name += "miss_vertex_unknown";
            }
        } else {
            miss_reco_name += "other";
        }
        cases_cont[miss_reco_name]++;

        miss_reco_name += "_" + to_string(cont_by_pdg[211]) + "pi_" + to_string(cont_by_pdg[13]) +"mu_"+to_string(cont_by_pdg[0]) +"o";
        cases_cont[miss_reco_name]++;

        if((miss_reco_name.find("0pi") != std::string::npos) && (miss_reco_name.find("211_") != std::string::npos)) classified = true;
        cont_reconstructed++;
        */


        
    }


    for (std::map<string,double>::iterator it=cases_cont.begin(); it!=cases_cont.end(); ++it)
        std::cout << it->first << ": " << it->second << '\n';


}