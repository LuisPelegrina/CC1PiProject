#include "../../Includes.h"

void shower_study_tl_e_razzled_score()
{
    GenerateDictionaries();
    TTree *tree;
    TFile *input;
    
    Cut_Parameters cut_p;
    cut_p.min_track_score = 0.5;
    
    gStyle->SetPalette(56);

    //Declare the variables
    Slice *slice = 0;

    string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
    input = new TFile(strRuta.c_str());
    tree =(TTree*)input->Get("tree");

    tree->SetBranchAddress("slice", &slice);


    TH1* h_pi0_shower_TL = new TH1D("h_pi0_shower_TL", "h_pi0_shower_KE", 100, 0, 100);
    TH1* h_other_shower_TL = new TH1D("h_other_shower_TL", "h_other_shower_KE", 100, 0, 100);
  TH1* h_pi0_shower_ts = new TH1D("h_pi0_shower_ts", "h_pi0_shower_ts", 100, 0, 1);
  TH1* h_other_shower_ts = new TH1D("h_other_shower_ts", "h_other_shower_ts", 100, 0, 1);

  TH1* h_pi0_shower_e_score = new TH1D("h_pi0_shower_e", "h_pi0_shower_e", 100, 0, 1);
  TH1* h_other_shower_e_score = new TH1D("h_other_shower_e", "h_other_shower_e", 100, 0, 1);

  TH1* h_pi0_shower_pho_score = new TH1D("h_pi0_shower_pho", "h_pi0_shower_pho", 100, 0, 1);
  TH1* h_other_shower_pho_score = new TH1D("h_other_shower_pho", "h_other_shower_pho", 100, 0, 1);


  TH2* h_2D_pi0 = new TH2D("h_pi0_shower_ts", "h_pi0_shower_ts", 100, 0, 1, 100, 0, 200);
  TH2* h_2D_other = new TH2D("h_other_shower_ts", "h_other_shower_ts", 100, 0, 1, 100, 0, 200);

  int nEntries = tree->GetEntries();
  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;

    int num_mu = slice->true_interaction.get_generator_num_primary_p(13);
    int num_pi = slice->true_interaction.get_generator_num_primary_p(211);
    int num_pi_0 = slice->true_interaction.get_generator_num_primary_p(111);
    bool is_CC1Pi = ((num_pi == 1) && (num_mu == 1));

    if(!is_CC1Pi) continue;
    for(Reco_Particle reco_p: slice->pandora_particle) {
      if((abs(reco_p.matched_pdg) != 11) && (abs(reco_p.matched_pdg) != 22)) continue;

      int parent_pdg = 0;
      for(G4Particle g4_p: slice->true_interaction.g4_particles) {
        if(g4_p.ID == reco_p.parent_ID) parent_pdg = abs(g4_p.PDG);
      }
      //cout << parent_pdg << " " << reco_p.matched_pdg << endl;


      if((abs(reco_p.matched_pdg == 22))) {
        h_pi0_shower_TL->Fill(reco_p.track_kinetic_energy);
        h_pi0_shower_ts->Fill(reco_p.track_score);
        h_pi0_shower_e_score->Fill(reco_p.razzled_score.electron_score);
        h_pi0_shower_pho_score->Fill(reco_p.razzled_score.photon_score);
        h_2D_pi0->Fill(reco_p.track_score, reco_p.track_kinetic_energy);

      } else {
        cout << reco_p.matched_pdg << endl;
        h_other_shower_TL->Fill(reco_p.track_kinetic_energy);
        h_other_shower_ts->Fill(reco_p.track_score);
        h_other_shower_e_score->Fill(reco_p.razzled_score.electron_score);
        h_other_shower_pho_score->Fill(reco_p.razzled_score.photon_score);
        h_2D_other->Fill(reco_p.track_score, reco_p.track_kinetic_energy);
      }

    }


  }

  TCanvas *c1 = new TCanvas();
  h_other_shower_TL->Draw("hist");
  h_pi0_shower_TL->Draw("hist same");
  h_pi0_shower_TL->SetLineColor(kRed +2);
  h_other_shower_TL->SetLineColor(kBlue +2);

  TCanvas *c2 = new TCanvas();
  h_other_shower_ts->Draw("hist");
  h_pi0_shower_ts->Draw("hist same");
  h_pi0_shower_ts->SetLineColor(kRed +2);
  h_other_shower_ts->SetLineColor(kBlue +2);


  TCanvas *c5 = new TCanvas();
  h_other_shower_e_score->Draw("hist");
  h_pi0_shower_e_score->Draw("hist same");
  h_pi0_shower_e_score->SetLineColor(kRed +2);
  h_other_shower_e_score->SetLineColor(kBlue +2);


  TCanvas *c6 = new TCanvas();
  h_other_shower_pho_score->Draw("hist");
  h_pi0_shower_pho_score->Draw("hist same");
  h_pi0_shower_pho_score->SetLineColor(kRed +2);
  h_other_shower_pho_score->SetLineColor(kBlue +2);

  TCanvas *c3 = new TCanvas();
  h_2D_other->Draw("hist");

  TCanvas *c4 = new TCanvas();
  h_2D_pi0->Draw("hist");
}