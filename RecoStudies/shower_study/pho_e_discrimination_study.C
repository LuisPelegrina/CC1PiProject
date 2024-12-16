#include "../../Includes.h"

void pho_e_discrimination_study()
{
    GenerateDictionaries();
    TTree *tree;
    TFile *input;
    
    Cut_Parameters cut_p;
    cut_p.min_track_score = 0.5;
    
    gStyle->SetPalette(56);

    //Declare the variables
    Slice *slice = 0;

    string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_82p7k.root";
    input = new TFile(strRuta.c_str());
    tree =(TTree*)input->Get("tree");

    tree->SetBranchAddress("slice", &slice);

    const int num_hist = bin_info.size();

    vector<vector<TH2*>> TH2_vec_e;
    vector<vector<TH2*>> TH2_vec_pho;
    for(int i_x = 0; i_x < num_hist; i_x++) {
      vector<TH2*> TH2_e_dummy;
      vector<TH2*> TH2_pho_dummy;
      for(int i_y = i_x+1; i_y < num_hist; i_y++) {
        string name = "h2_e_" + to_string(i_x) + "_" + to_string(i_y);
        TH2_e_dummy.push_back(new TH2D(name.c_str(), name.c_str(),
                                   bin_info[i_x].nBins, bin_info[i_x].LowBin, bin_info[i_x].UpBin,
                                   bin_info[i_y].nBins, bin_info[i_y].LowBin, bin_info[i_y].UpBin));

        name = "h2_pho_" + to_string(i_x) + "_" + to_string(i_y);
        TH2_pho_dummy.push_back(new TH2D(name.c_str(), name.c_str(),
                                        bin_info[i_x].nBins, bin_info[i_x].LowBin, bin_info[i_x].UpBin,
                                        bin_info[i_y].nBins, bin_info[i_y].LowBin, bin_info[i_y].UpBin));
      }
      TH2_vec_e.push_back(TH2_e_dummy);
      TH2_vec_pho.push_back(TH2_pho_dummy);
    }


  TH1 *h_e[num_hist];
  TH1 *h_pho[num_hist];
  for (int i_hist = 0; i_hist < num_hist; i_hist++) {
    string title = "h_e" + to_string(i_hist);
    h_e[i_hist] = new TH1D(title.c_str(), " ", bin_info[i_hist].nBins, bin_info[i_hist].LowBin,
                                          bin_info[i_hist].UpBin);
    title = "h_pho" + to_string(i_hist);
    h_pho[i_hist] = new TH1D(title.c_str(), " ", bin_info[i_hist].nBins, bin_info[i_hist].LowBin,
                                bin_info[i_hist].UpBin);
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

    if(!is_CC1Pi) continue;
    for(Reco_Particle reco_p: slice->pandora_particle) {
      if((abs(reco_p.matched_pdg) != 22) && (abs(reco_p.matched_pdg) != 11)) continue;
      //if(reco_p.track_score > 0.5) continue;
      //if(reco_p.purity < 0.8) continue;
      //if(reco_p.completeness < 0.8) continue;

      int parent_pdg = 0;
      for(G4Particle g4_p: slice->true_interaction.g4_particles) {
        if(g4_p.ID == reco_p.parent_ID) parent_pdg = abs(g4_p.PDG);
      }
      //cout << parent_pdg << " " << reco_p.matched_pdg << endl;

      Data["track_length"] = reco_p.track_lenght;
      Data["track_score"] = reco_p.track_score;
      Data["e_score"] = reco_p.razzled_score.electron_score;
      Data["pho_score"] = reco_p.razzled_score.photon_score;
      Data["dEdx"] = reco_p.shower_dEdx;
      Data["opening_angle"] = reco_p.shower_opening_angle;
      Data["shower_length"] = reco_p.shower_length;
      Data["shower_energy"] = reco_p.shower_energy;

      if(abs(reco_p.matched_pdg) == 11) {
        for (int i_hist = 0; i_hist < num_hist; i_hist++) {
          h_e[i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
        }
        //h_2D_shower->Fill(reco_p.track_score, reco_p.track_kinetic_energy);

        for(int i_x = 0; i_x < num_hist; i_x++) {
          for(int i_y = 0; i_y < num_hist-(i_x+1); i_y++) {
            TH2_vec_e.at(i_x).at(i_y)->Fill(Data[bin_info[i_x].FillDataType], Data[bin_info[i_y + i_x +1].FillDataType]);
          }
        }
      } else {
        cout << reco_p.matched_pdg << endl;
        for (int i_hist = 0; i_hist < num_hist; i_hist++) {
          h_pho[i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
        }

        for(int i_x = 0; i_x < num_hist; i_x++) {
          for(int i_y = 0; i_y < num_hist-(i_x+1); i_y++) {
            TH2_vec_pho.at(i_x).at(i_y)->Fill(Data[bin_info[i_x].FillDataType], Data[bin_info[i_y + i_x +1].FillDataType]);
          }
        }
        //h_2D_miss_shower->Fill(reco_p.track_score, reco_p.track_kinetic_energy);
      }

    }


  }


  string folder = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/shower_studies/e_vs_pho/";
  TLegend *leg = new TLegend(0.55, 0.6, .87, .87);

  leg->AddEntry(h_e[0], "e^{-}  distribution", "l");
  leg->AddEntry(h_pho[0], "#gamma distribution", "l");
  for (int i_hist = 0; i_hist < num_hist; i_hist++) {
    TCanvas *ci = new TCanvas();
    h_e[i_hist]->SetStats(0);
    h_e[i_hist]->Draw("hist");
    h_pho[i_hist]->Draw("hist same");
    h_e[i_hist]->SetTitle((";" + bin_info[i_hist].Title + "; Events").c_str());

    h_pho[i_hist]->SetLineColor(kRed +2);
    h_e[i_hist]->SetLineColor(kBlue +2);
    leg->Draw();
    ci->SaveAs((folder + "TH1/ "+ bin_info.at(i_hist).FillDataType +  ".pdf").c_str());
    ci->Close();
    delete ci;
  }

  gStyle->SetPadRightMargin(0.175);
  for(int i_x = 0; i_x < num_hist; i_x++) {
    for(int i_y = 0; i_y < num_hist-(i_x+1); i_y++) {
      string name = "ci" + to_string(i_x) + "_" + to_string(i_y);
      TCanvas *ci = new TCanvas(name.c_str(), "ci", 2000, 800);
      ci->Divide(2,1);
      ci->cd(1);
      TH2_vec_e.at(i_x).at(i_y)->SetTitle(("e;" + bin_info[i_x].Title + ";" + bin_info[i_y + i_x +1].Title + "; Events").c_str());
      TH2_vec_e.at(i_x).at(i_y)->SetStats(0);
      TH2_vec_e.at(i_x).at(i_y)->Draw("colz");
      ci->cd(2);
      TH2_vec_pho.at(i_x).at(i_y)->SetTitle(("#gamma;" + bin_info[i_x].Title + ";" + bin_info[i_y + i_x +1].Title + "; Events").c_str());
      TH2_vec_pho.at(i_x).at(i_y)->SetStats(0);
      TH2_vec_pho.at(i_x).at(i_y)->Draw("colz");

      ci->SaveAs((folder + "TH2/" + bin_info.at(i_x).FillDataType + "_" + bin_info.at(i_y + i_x +1).FillDataType + ".pdf").c_str());
      ci->Close();
      delete ci;
    }
  }

  /*
  TCanvas *c3 = new TCanvas();
  h_2D_shower->Draw("hist");

  TCanvas *c4 = new TCanvas();
  h_2D_miss_shower->Draw("hist");
*/
}