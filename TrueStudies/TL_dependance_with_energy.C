#include "../Includes.h"

void TL_dependance_with_energy()
{
    GenerateDictionaries();
    TTree *tree;
    TFile *input;
    
    Cut_Parameters cut_p;
    cut_p.min_track_lenght = 0;
    cut_p.min_track_score = 0;
    
    gStyle->SetPalette(56);

    //Declare the variables
    Slice *slice = 0;

    string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
    input = new TFile(strRuta.c_str());
    tree =(TTree*)input->Get("tree");

    tree->SetBranchAddress("slice", &slice);
    int nEntries = tree->GetEntries();


    TH2* h_muon_P0_TL = new TH2D("h_muon_P0_TL", "h_muon_P0_TL", 100, 0, 0.5, 50, 0, 50);
    TH2* h_pion_P0_TL = new TH2D("h_pion_P0_TL", "h_pion_P0_TL", 100, 0, 0.5, 50, 0, 50);
  TH2* h_muon_KE_TL = new TH2D("h_muon_KE_TL", "h_muon_KE_TL", 100, 0, 0.2, 50, 0, 50);
  TH2* h_pion_KE_TL = new TH2D("h_muon_KE_TL", "h_muon_KE_TL", 100, 0, 0.2, 50, 0, 50);

    int cont_total = 0;
    for(int i_e = 0; i_e < nEntries; ++i_e) {
        tree->GetEntry(i_e);
        if(i_e%100 == 0) cout << "Entry:" << i_e << endl;



        cont_total++;
        if(slice->slice_ID != 0) continue;
        if(!slice->true_interaction.is_selected_final_state("CC1Pi",0.325)) continue;


        G4Particle primary_muon = slice->true_interaction.get_g4_primary_p(13).at(0);
        G4Particle primary_pion = slice->true_interaction.get_g4_primary_p(211).at(0);

        h_muon_P0_TL->Fill(primary_muon.P0.Mag(), primary_muon.TL);
        h_pion_P0_TL->Fill(primary_pion.P0.Mag(), primary_pion.TL);

      vector<G4Particle> primary_proton = slice->true_interaction.get_g4_primary_p(2212);
      for(G4Particle g4_p: primary_proton) {
        h_muon_KE_TL->Fill(sqrt(pow(g4_p.P0.Mag(),2) + pow(0.938272,2)) - 0.938272 , g4_p.TL);

      }
      //h_muon_KE_TL->Fill(sqrt(pow(primary_muon.P0.Mag(),2) + pow(0.13957,2)) - 0.13957 , primary_muon.TL);
      h_pion_KE_TL->Fill( sqrt(pow(primary_pion.P0.Mag(),2) + pow(0.10566,2)) - 0.10566, primary_pion.TL);

    }
    TCanvas *c1 = new TCanvas();
    h_muon_P0_TL->SetTitle("Missed Muon TL; P_{0} [GeV]; TL [cm]");
    h_muon_P0_TL->Draw("hist");

    TCanvas *c2 = new TCanvas();
    h_pion_P0_TL->SetTitle("Missed Pion TL; P_{0} [GeV]; TL [cm]");
    h_pion_P0_TL->Draw("hist");

    TCanvas *c3 = new TCanvas();
    h_muon_KE_TL->SetTitle("Missed Muon TL; Kinetic Energy [GeV]; TL [cm]");
    h_muon_KE_TL->Draw("hist");

    TCanvas *c4 = new TCanvas();
    h_pion_KE_TL->SetTitle("Missed Pion TL; Kinetic Energy [GeV]; TL [cm]");
    h_pion_KE_TL->Draw("hist");


}