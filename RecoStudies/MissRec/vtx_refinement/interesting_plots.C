#include "../../../Includes.h"
#include "../Graphs_utils.cpp"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}

void interesting_plots()
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

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_83k.root";
  input = new TFile(strRuta.c_str());
  tree =(TTree*)input->Get("tree");
  tree->SetBranchAddress("slice", &slice);
  int nEntries = tree->GetEntries();

  vector<struct LocalBinInformation> bin_info = {
    {"m_0", "m_0", 50, 0, 10000},
    {"min_0", "min_0", 50, 0, 10000},
    {"max_0", "max_0", 50, 0, 10000},
    {"n_pfps", "n_pfps", 6, 0, 6},
    {"n_good_pfps", "n_good_pfps", 6, 0, 6},
    {"longest_TL", "longest_TL", 20, 0, 100},
    {"longest_track_theta_z", "longest_track_theta_z", 20, 0, 180}
  };

  TH1 *h_good_vtx = new TH1D("h_good_vtx", "h_good_vtx",20,0,1000);
  TH1 *h_miss_vtx = new TH1D("h_miss_vtx", "h_miss_comp",20,0,1000);
  TH1 *h_random_slice_comp = new TH1D("h_random_slice_comp", "h_miss_comp",20,0,180);

cout << bin_info.size()<< endl << endl;
  int num_hist = bin_info.size();
  TH1 *h_good[num_hist];
  TH1 *h_miss[num_hist];
  TH1 *h_random_slice[num_hist];
  for(int i_hist = 0; i_hist < num_hist; i_hist++) {
    string title = "h_" + bin_info.at(i_hist).FillDataType;
    h_good[i_hist] = new TH1D((title + "_good").c_str(),(title + "_good").c_str(), bin_info[i_hist].nBins, bin_info[i_hist].LowBin, bin_info[i_hist].UpBin);
    h_miss[i_hist] = new TH1D((title + "_miss").c_str(),(title + "_good").c_str(), bin_info[i_hist].nBins, bin_info[i_hist].LowBin, bin_info[i_hist].UpBin);
    h_random_slice[i_hist] = new TH1D((title + "_random_slice").c_str(),(title + "_good").c_str(), bin_info[i_hist].nBins, bin_info[i_hist].LowBin, bin_info[i_hist].UpBin);
  }


  const int num_cuts = 3;
  string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut"};

  map<string, double> Data;
  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;

    //if(slice->slice_ID != 0) continue;
    if (!slice->true_interaction.is_selected_final_state("N_CC1Pi")) continue;
    if (!slice->true_interaction.is_vertex_contained(10)) continue;

    bool pass_cut = true;
    //Plot only if correctly reconstructed
    for (int i_cut = 0; i_cut < num_cuts; i_cut++) {
      if (!pass_cut) continue;
      if (!slice->pass_analysis_cut(cut_p, cut_name[i_cut])) {
        pass_cut = false;
      }
    }
    if (!pass_cut) continue;

    G4Particle primary_muon = slice->true_interaction.get_g4_primary_p(13).at(0);
    G4Particle primary_pion = slice->true_interaction.get_g4_primary_p(211).at(0);

    int cont_primary_p = 0;
    int cont_p = slice->pandora_particle.size();
    map<int, int> cont_by_pdg_pur_80;
    map<int, int> cont_by_pdg_pur_50;
    map<int, int> cont_by_pdg;

    int close_muon = 0;
    int close_pion = 0;
    for (int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
      Reco_Particle pandora_p = slice->pandora_particle.at(i_p);
      if (slice->pandora_particle.at(i_p).is_pandora_primary) {
        cont_by_pdg[11111]++;
        cont_by_pdg[abs(pandora_p.matched_pdg)]++;
        if ((pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag() < 10) {
          cont_by_pdg[22222]++;
          if (pandora_p.matched_pdg == 211) {
            close_pion++;
          } else if (pandora_p.matched_pdg == 13) {
            close_muon++;
          }
        }
        if ((pandora_p.completeness > 0.8) && (pandora_p.purity > 0.8))
          cont_by_pdg_pur_80[abs(pandora_p.matched_pdg)]++;
        if ((pandora_p.completeness > 0.5) && (pandora_p.purity > 0.5))
          cont_by_pdg_pur_50[abs(pandora_p.matched_pdg)]++;
        if ((abs(pandora_p.matched_pdg) == 211) || (abs(pandora_p.matched_pdg) == 13)) continue;
        cont_by_pdg[-1]++;
      }
    }

    double distance_to_true_vertex = (slice->primary_vertex_true.vertex_cm -
                                      slice->primary_vertex_reco.vertex_cm).Mag();
    bool show = true;


    if ((primary_muon.TL < 3) || (primary_pion.TL < 3)) continue;
    if (distance_to_true_vertex < 2) continue;

    double run_Id = slice->run_ID;
    double subrun_Id = slice->subrun_ID;
    double event_Id = slice->event_ID;

    //Check if there is an extra slice or not
    bool extra_slice = false;

    tree->GetEntry(i_e - 1);
    if ((run_Id = slice->run_ID) && (subrun_Id = slice->subrun_ID) && (event_Id = slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e + 1);
    if ((run_Id = slice->run_ID) && (subrun_Id = slice->subrun_ID) && (event_Id = slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e - 2);
    if ((run_Id = slice->run_ID) && (subrun_Id = slice->subrun_ID) && (event_Id = slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e + 2);
    if ((run_Id = slice->run_ID) && (subrun_Id = slice->subrun_ID) && (event_Id = slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e);

    int mother_PDG = -1;
    bool first = false;
    if (extra_slice) {
      for (int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
        if ((slice->pandora_particle.at(i_p).track_start - slice->primary_vertex_reco.vertex_cm).Mag() < 2) {
          if (first) continue;
          first = true;
          for (int i_pg4 = 0; i_pg4 < slice->true_interaction.g4_particles.size(); i_pg4++) {
            G4Particle g4_p = slice->true_interaction.g4_particles.at(i_pg4);
            if (slice->pandora_particle.at(i_p).true_track_id == g4_p.ID) {
              for (int i_pg4_2 = 0; i_pg4_2 < slice->true_interaction.g4_particles.size(); i_pg4_2++) {
                if (slice->true_interaction.g4_particles.at(i_pg4_2).ID == g4_p.mother) {
                  mother_PDG = slice->true_interaction.g4_particles.at(i_pg4_2).PDG;
                }
              }
            }
          }
        }
      }
    }

    bool has_possible_vtx = false;
    for (int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      if ((slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag() < 4) has_possible_vtx = true;
    }

    for (int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      TH1 *h = slice->build_dQdx_hist(slice->vertex_vec.at(i_v).vertex_wire, 2, 0, 100, i_v, 100);
      Data["m_0"] = slice->get_dQdx_segment_info(h, 1, 10).mean;
      Data["min_0"] = slice->get_dQdx_segment_info(h, 1, 10).min;
      Data["max_0"] = slice->get_dQdx_segment_info(h, 1, 10).max;

      int n_pfps = 0;
      int n_good_pfps = 0;
      double longest_TL = 0;
      double longest_track_theta_z = 0;
      for(int i_p = 0; i_p < slice->pandora_particle.size();i_p++) {
        Reco_Particle pandora_p = slice->pandora_particle.at(i_p);

        int num_hits_U = 0;
        int num_hits_V = 0;
        int num_hits_C = 0;
        int total_hits = 0;
        for (int i_h = 0; i_h < slice->hits.size(); i_h++) {
          Hit hit = slice->hits.at(i_h);
          if (hit.associated_pfp_ID == pandora_p.ID) {
            if (hit.Plane_ID == 2) num_hits_C++;
            if (hit.Plane_ID == 1) num_hits_U++;
            if (hit.Plane_ID == 0) num_hits_V++;
            if (hit.Plane_ID != -1) total_hits++;
          }
        }
        int planes_ok = 0;
        if (num_hits_U >= 5) planes_ok++;
        if (num_hits_V >= 5) planes_ok++;
        if (num_hits_C >= 5) planes_ok++;



        if (((pandora_p.track_start - slice->vertex_vec.at(i_v).vertex_cm).Mag() < 2) ||
        (pandora_p.track_end - slice->vertex_vec.at(i_v).vertex_cm).Mag() < 2) {
          n_pfps++;
          if( pandora_p.track_lenght > longest_TL)  {
            longest_TL = pandora_p.track_lenght;
            longest_track_theta_z = pandora_p.start_direction.Angle(TVector3(0,0,1))*360/(2*TMath::Pi());
            //cout << pandora_p.track_theta<< " " << pandora_p.start_direction.Angle(TVector3(0,0,1)) << endl;
            //cout << longest_track_theta_z << endl;
          }
          if ((total_hits > 15) && ((planes_ok >= 2) || num_hits_C > 15) && (pandora_p.track_lenght > 3) &&
              (pandora_p.chi2_score.muon_score != -1)) {
                 n_good_pfps++;
          }
        }
      }
      Data["n_pfps"] = n_pfps;
      Data["n_good_pfps"] = n_good_pfps;
      Data["longest_TL"] = longest_TL;
      Data["longest_track_theta_z"] = longest_track_theta_z;

      //if(n_pfps != 2) continue;

      if(extra_slice && ((mother_PDG == 2112) || (mother_PDG == 111)) && !has_possible_vtx) {
        for(int i_hist = 0; i_hist < num_hist; i_hist++) h_random_slice[i_hist]->Fill(Data[bin_info.at(i_hist).FillDataType]);
      } else if ((slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag() < 4) {
        for(int i_hist = 0; i_hist < num_hist; i_hist++) h_good[i_hist]->Fill(Data[bin_info.at(i_hist).FillDataType]);
      } else {
        for(int i_hist = 0; i_hist < num_hist; i_hist++) h_miss[i_hist]->Fill(Data[bin_info.at(i_hist).FillDataType]);
      }
      delete h;
    }

  }

  cout << "NICE" << endl;
  for(int i_hist = 0; i_hist < num_hist; i_hist++) {
    TCanvas *ci = new TCanvas();
    h_miss[i_hist]->SetLineColor(kRed +2);
    h_random_slice[i_hist]->SetLineColor(kGreen +2);
    double Max =h_good[i_hist]->GetMaximum();
    if(h_miss[i_hist]->GetMaximum()> Max) Max =h_miss[i_hist]->GetMaximum() ;
    if(h_random_slice[i_hist]->GetMaximum()> Max) Max =h_random_slice[i_hist]->GetMaximum() ;
    h_good[i_hist]->SetMaximum(Max*1.2);
    h_good[i_hist]->Draw("hist");
    h_miss[i_hist]->Draw("hist same");
    h_random_slice[i_hist]->Draw("hist same");
   }


}