#include "../../Includes.h"
#include "Graphs_utils.cpp"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}

void TL_missing_graphs()
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

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_84k.root";
  input = new TFile(strRuta.c_str());
  tree =(TTree*)input->Get("tree");
  tree->SetBranchAddress("slice", &slice);
  int nEntries = tree->GetEntries();

  TH2* missing_muon_TL_MissV = new TH2D("missing_muon_TL_MissV", "missing_muon_TL_MissV", 20, 0, 20, 20, 0, 20);
  TH2* missing_pion_TL_MissV = new TH2D("missing_pion_TL_MissV", "missing_pion_TL_MissV", 20, 0, 20, 20, 0, 20);
  TH2* missing_muon_TL_hist = new TH2D("missing_muon_TL_hist", "missing_muon_TL_hist", 20, 0, 20,2, 0, 2);
  TH2* missing_pion_TL_hist = new TH2D("missing_pion_TL_hist", "missing_pion_TL_hist", 20, 0, 20,2, 0, 2);

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
    if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;

    if (slice->slice_ID != 0) continue;
    if (!slice->true_interaction.is_selected_final_state("CC1Pi")) continue;

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
        if ((abs(pandora_p.matched_pdg) == 211) || (abs(pandora_p.matched_pdg) == 13)) continue;
        cont_by_pdg[abs(pandora_p.matched_pdg)]++;
        cont_by_pdg[-1]++;
      }
    }

    double distance_to_true_vertex = (slice->primary_vertex_true.vertex_cm -
                                      slice->primary_vertex_reco.vertex_cm).Mag();
    cases_cont["tot"] += 1;
    string case_name = "";

    if (cont_by_pdg[13] == 0)  {
      missing_muon_TL_MissV->Fill(primary_muon.TL, distance_to_true_vertex);
      missing_muon_TL_hist->Fill(primary_muon.TL, 1);
    } else {
      missing_muon_TL_hist->Fill(primary_muon.TL, 0);

    }

    if (cont_by_pdg[211] == 0)  {
      missing_pion_TL_MissV->Fill(primary_pion.TL, distance_to_true_vertex);
      missing_pion_TL_hist->Fill(primary_pion.TL, 1);
    } else {
      missing_pion_TL_hist->Fill(primary_pion.TL, 0);

    }
  }


 TCanvas *c1 = new TCanvas();
  c1->cd();
  missing_muon_TL_MissV->Draw("colz");

  for(int i_b = 1; i_b <= missing_muon_TL_hist->GetNbinsX(); i_b++) {
    double colum_down_value = missing_muon_TL_hist->GetBinContent(missing_muon_TL_hist->GetBin(i_b, 1));
    double colum_up_value = missing_muon_TL_hist->GetBinContent(missing_muon_TL_hist->GetBin(i_b, 2));
    double colum_total = colum_up_value + colum_down_value;
    missing_muon_TL_hist->SetBinContent(missing_muon_TL_hist->GetBin(i_b, 1),colum_down_value/colum_total);
    missing_muon_TL_hist->SetBinContent(missing_muon_TL_hist->GetBin(i_b, 2),colum_up_value/colum_total);
  }


  TCanvas *c2 = new TCanvas();
  c2->cd();
  missing_muon_TL_hist->Draw("colz");

  for(int i_b = 1; i_b <= missing_pion_TL_hist->GetNbinsX(); i_b++) {
    double colum_down_value = missing_pion_TL_hist->GetBinContent(missing_pion_TL_hist->GetBin(i_b, 1));
    double colum_up_value = missing_pion_TL_hist->GetBinContent(missing_pion_TL_hist->GetBin(i_b, 2));
    double colum_total = colum_up_value + colum_down_value;
    missing_pion_TL_hist->SetBinContent(missing_pion_TL_hist->GetBin(i_b, 1),colum_down_value/colum_total);
    missing_pion_TL_hist->SetBinContent(missing_pion_TL_hist->GetBin(i_b, 2),colum_up_value/colum_total);
  }
  TCanvas *c3 = new TCanvas();
  c3->cd();
  missing_pion_TL_MissV->Draw("colz");

  TCanvas *c4 = new TCanvas();
  c4->cd();
  missing_pion_TL_hist->Draw("colz");
}