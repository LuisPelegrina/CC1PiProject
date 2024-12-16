#include "../../Includes.h"
#include "Graphs_utils.cpp"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}

void min_TL_dtov_TH2D()
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

  double max_TL = 6;
  double min_TL = 1;
  int n_TL = 20;
  double pas_TL = (max_TL - min_TL)/n_TL;


  double max_dtv = 10;
  double min_dtv = 5;
  int n_dtv = 20;
  double pas_dtv = (max_dtv - min_dtv)/n_dtv;
  TH2 *h_lost = new TH2D("h_lost","h_lost", n_TL, min_TL, max_TL,n_dtv, min_dtv, max_dtv);
  TH2 *h_recovered = new TH2D("h_recovered","h_recovered", n_TL, min_TL, max_TL,n_dtv, min_dtv, max_dtv);

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_83k.root";
  input = new TFile(strRuta.c_str());
  tree =(TTree*)input->Get("tree");
  tree->SetBranchAddress("slice", &slice);
  int nEntries = tree->GetEntries();

  const int num_cuts = 4;
  string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut"};


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

    if ((primary_muon.TL < 3) || (primary_pion.TL < 3)) {
      continue;
    }


    for(int i_tl = 0; i_tl < n_TL; i_tl++) {
      double TL = i_tl*pas_TL + min_TL;
      for(int i_dtv = 0; i_dtv< n_dtv; i_dtv++) {
        double dtv = i_dtv * pas_dtv + min_dtv;
        map<int, int> cont_by_pdg;
        map<int, int> cont_by_pdg_after_corr;

        for(int i_p = 0; i_p < slice->pandora_particle.size() ;i_p++) {
          Reco_Particle pandora_p = slice->pandora_particle.at(i_p);
          if(slice->pandora_particle.at(i_p).is_pandora_primary){
            cont_by_pdg[11111]++;
            cont_by_pdg[abs(pandora_p.matched_pdg)]++;

            int num_hits_U = 0;
            int num_hits_V = 0;
            int num_hits_C = 0;
            int total_hits = 0;
            for(int i_h = 0; i_h < slice->hits.size(); i_h++){
              Hit hit = slice->hits.at(i_h);
              if(hit.associated_pfp_ID == pandora_p.ID) {
                if(hit.Plane_ID == 2) num_hits_C++;
                if(hit.Plane_ID == 1) num_hits_U++;
                if(hit.Plane_ID == 0) num_hits_V++;
                if(hit.Plane_ID != -1) total_hits++;
              }
            }
            int planes_ok = 0;
            if(num_hits_U >= 5) planes_ok++;
            if(num_hits_V >= 5) planes_ok++;
            if(num_hits_C >= 5) planes_ok++;

            if((total_hits > 15) && ((planes_ok >= 2)|| num_hits_C > 15) && (pandora_p.track_lenght > TL) && (pandora_p.chi2_score.muon_score != -1)
               && (pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag() < dtv) {
              cont_by_pdg_after_corr[11111]++;
              cont_by_pdg_after_corr[abs(pandora_p.matched_pdg)]++;
            }
           }
        }
        bool good_reco_start = false;

        if(((cont_by_pdg[211] == 1) && (cont_by_pdg[13] == 1)) && (cont_by_pdg[11111] == 2)) {
          good_reco_start = true;
        } else {
          good_reco_start = false;
        }

        if(((cont_by_pdg_after_corr[211] == 1) && (cont_by_pdg_after_corr[13] == 1)) && (cont_by_pdg_after_corr[11111] == 2)) {
          if(!good_reco_start) h_recovered->Fill(TL, dtv, 1);
        } else {
          if(good_reco_start) h_lost->Fill(TL, dtv, 1);
        }
      }
    }


  }

  string path = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/TLvsDTVRecoLostComparison/";


  TCanvas *c1 = new TCanvas();
  h_recovered->SetTitle("Recovered;TL [cm]; Distance to vertex [cm]");
  h_recovered->SetStats(0);
  h_recovered->Draw("colz");
  c1->SaveAs((path+ "Recovered.pdf").c_str());

  TCanvas *c2 = new TCanvas();
  h_lost->SetTitle("Lost;TL [cm]; Distance to vertex [cm]");
  h_lost->SetStats(0);
  h_lost->Draw("colz");
  c2->SaveAs((path+ "Lost.pdf").c_str());


  TCanvas *c3 = new TCanvas();
  TH2* h_difference = (TH2F*)h_recovered->Clone("h_difference");
  h_difference->Add(h_lost, -1);
  h_difference->SetTitle("Diference;TL [cm]; Distance to vertex [cm]");
  h_difference->Draw("colz");
  c3->SaveAs((path+ "Difference.pdf").c_str());

}