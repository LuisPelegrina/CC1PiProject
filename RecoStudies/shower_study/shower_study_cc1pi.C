#include "../../Includes.h"

void shower_study()
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


    TH1* h_cc1pi_shower_KE = new TH1D("h_muon_P0_TL", "h_muon_P0_TL", 100, 0, 1000);
    TH1* h_other_shower_KE = new TH1D("h_pion_P0_TL", "h_pion_P0_TL", 100, 0, 1000);

  int nEntries = tree->GetEntries();
  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;
    if (!slice->true_interaction.is_selected_final_state("CC1Pi", 0.325)) {
      for (int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
        Reco_Particle reco_particle = slice->pandora_particle.at(i_p);

        if (!reco_particle.is_shower(cut_p)) continue;


        h_cc1pi_shower_KE->Fill(reco_particle.track_kinetic_energy);
      }

    } else {
      for (int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
        Reco_Particle reco_particle = slice->pandora_particle.at(i_p);

        if (!reco_particle.is_shower(cut_p)) continue;
        h_other_shower_KE->Fill(reco_particle.track_kinetic_energy);
      }
    }
  }

  TCanvas *c1 = new TCanvas();
  h_cc1pi_shower_KE->Draw("hist");
  h_other_shower_KE->Draw("hist same ");
  h_other_shower_KE->SetLineColor(kRed +2);
  h_cc1pi_shower_KE->SetLineColor(kBlue +2);

}