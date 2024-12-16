#include "../../Includes.h"

void Chi2_PID_p_study()
{
    GenerateDictionaries();
    TTree *tree;
    TFile *input;

    gStyle->SetPalette(56);

    //Declare the variables
    Slice *slice = 0;

    string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_82p7k.root";
    input = new TFile(strRuta.c_str());
    tree =(TTree*)input->Get("tree");

    tree->SetBranchAddress("slice", &slice);



  vector<struct LocalBinInformation> bin_info = {
    {"track_length", "Track length [cm]", 50, 0, 100},
    {"p_chi2", "proton #chi^{2}", 50, 0, 200},
    {"mu_chi2", "#muon #chi^{2}", 50, 0, 200},
    {"pi_chi2", "#pion #chi^{2}", 50, 0, 200},
    {"chi2_diff", "#pion #chi^{2} - #muon #chi^{2}", 50, -100, 100},
  };
  const int num_hist = bin_info.size();

  vector<int> pdg_vec = {2212,211,13};
  int num_pdg = pdg_vec.size();

  map<int,vector<TH1*>> h_proton_id;
  map<int,vector<TH1*>> h_muon_id;
  vector<TH1*> dummy_vec;
  for (int i_hist = 0; i_hist < num_hist; i_hist++) {
    dummy_vec.push_back(new TH1D());
  }

  for (int i_pdg = 0; i_pdg < num_pdg; i_pdg++) {
    h_proton_id[pdg_vec.at(i_pdg)] = dummy_vec;
    h_muon_id[pdg_vec.at(i_pdg)] = dummy_vec;
    for (int i_hist = 0; i_hist < num_hist; i_hist++) {
      cout << i_hist << " " << i_pdg << endl;
      string title = "h_proton" + to_string(i_hist)+ "_" +  to_string(i_pdg);
      h_proton_id[pdg_vec.at(i_pdg)][i_hist] = new TH1D(title.c_str(), " ", bin_info[i_hist].nBins, bin_info[i_hist].LowBin,
                                  bin_info[i_hist].UpBin);
       title = "h_muon" + to_string(i_hist)+ "_" +  to_string(i_pdg);
      h_muon_id[pdg_vec.at(i_pdg)][i_hist] = new TH1D(title.c_str(), " ", bin_info[i_hist].nBins, bin_info[i_hist].LowBin,
                                                      bin_info[i_hist].UpBin);
    }
  }
  cout << "NICE" << endl;

  map<string, double> Data;

  //Only select events that pass the cuts
  vector<string> cut_name = {"no_cut", "is_clear_cosmic_cut", "reco_cut" , "fv_cut",  "crumbs_cut", "track_cut",
                             "chi2_muon_like_cut", "shower_cut", "at_most_one_pfp_scapping_cut"};

  const int num_cuts = cut_name.size();
  Cut_Parameters cut_p;
  cut_p.use_default_cut_p();



  int nEntries = tree->GetEntries();
  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;

    bool pass_cut = true;
    for (int i_pcut = 0; i_pcut < num_cuts; i_pcut++) {
      if (!pass_cut) continue;
      if (!slice->pass_analysis_cut(cut_p, cut_name[i_pcut])) {
        pass_cut = false;
      }
    }
    if (!pass_cut) continue;

    for(Reco_Particle reco_p: slice->pandora_particle) {
      if((abs(reco_p.matched_pdg) != 211) && (abs(reco_p.matched_pdg) != 2212) && (abs(reco_p.matched_pdg) != 13)) continue;

      Data["track_length"] = reco_p.track_lenght;
      Data["p_chi2"] = reco_p.chi2_score.proton_score;
      Data["pi_chi2"] = reco_p.chi2_score.pion_score;
      Data["mu_chi2"] = reco_p.chi2_score.muon_score;
      Data["chi2_diff"] =  reco_p.chi2_score.muon_score -  reco_p.chi2_score.pion_score;
      double muon_like_score = min(reco_p.chi2_score.pion_score, reco_p.chi2_score.muon_score);

      bool is_candidate = false;
      if((reco_p.track_lenght > cut_p.chi2_min_TL) && (reco_p.chi2_score.proton_score > cut_p.chi2_max_proton_score) && ( reco_p.chi2_score.muon_score < cut_p.chi2_min_muon_score)) is_candidate = true;

      if(is_candidate) {
        for (int i_hist = 0; i_hist < num_hist; i_hist++) {
          h_muon_id[abs(reco_p.matched_pdg)][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
        }
      } else {
        for (int i_hist = 0; i_hist < num_hist; i_hist++) {
          h_proton_id[abs(reco_p.matched_pdg)][i_hist]->Fill(Data[bin_info[i_hist].FillDataType]);
        }
      }
    }
  }


  string folder = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/VariablesAfterChi2/";
  TLegend *leg = new TLegend(0.55, 0.6, .87, .87);
  TLegend *leg_2 = new TLegend(0.2, 0.7, .42, .87);
  TIter nextEntry(leg->GetListOfPrimitives());
  TLegendEntry *entry;
  while ((entry = (TLegendEntry*)nextEntry())) {
    leg_2->AddEntry(entry->GetObject(), entry->GetLabel(), entry->GetOption());
  }

  leg->AddEntry( h_proton_id[2212][0], "p distribution", "l");
  leg->AddEntry( h_proton_id[211][0], "pi distribution", "l");
  leg->AddEntry( h_proton_id[13][0], "mu distribution", "l");
  for (int i_hist = 0; i_hist < num_hist; i_hist++) {
    TCanvas *ci = new TCanvas();
    double max = 0;
    for (int i_pdg = 0; i_pdg < num_pdg; i_pdg++) {
      if(i_pdg == 0) h_proton_id[pdg_vec.at(i_pdg)][i_hist]->Draw("hist");
      if(i_pdg != 0) h_proton_id[pdg_vec.at(i_pdg)][i_hist]->Draw("hist same");
      h_proton_id[pdg_vec.at(i_pdg)][i_hist]->SetLineColor(colors_alpha.at(i_pdg+1));
      h_proton_id[pdg_vec.at(i_pdg)][i_hist]->SetStats(0);
      if(h_proton_id[pdg_vec.at(i_pdg)][i_hist]->GetMaximum() > max) max = h_proton_id[pdg_vec.at(i_pdg)][i_hist]->GetMaximum();
    }
    h_proton_id[pdg_vec.at(0)][i_hist]->SetTitle(("Identified as p;" + bin_info[i_hist].Title + "; Events").c_str());
    h_proton_id[pdg_vec.at(0)][i_hist]->GetYaxis()->SetRangeUser(0, max*1.2);
    if(bin_info[i_hist].FillDataType.compare("p_chi2") == 0) {
      leg_2->Draw();
    } else {
      leg->Draw();
    }
    //ci->SaveAs((folder + "proton/ "+ bin_info.at(i_hist).FillDataType +  ".pdf").c_str());
    //ci->Close();
    //delete ci;
  }

  for (int i_hist = 0; i_hist < num_hist; i_hist++) {
    double max = 0;
    TCanvas *ci = new TCanvas();
    for (int i_pdg = 0; i_pdg < num_pdg; i_pdg++) {
      if(i_pdg == 0) h_muon_id[pdg_vec.at(i_pdg)][i_hist]->Draw("hist");
      if(i_pdg != 0) h_muon_id[pdg_vec.at(i_pdg)][i_hist]->Draw("hist same");
      h_muon_id[pdg_vec.at(i_pdg)][i_hist]->SetLineColor(colors_alpha.at(i_pdg+1));
      h_muon_id[pdg_vec.at(i_pdg)][i_hist]->SetStats(0);
      if(h_muon_id[pdg_vec.at(i_pdg)][i_hist]->GetMaximum() > max) max = h_muon_id[pdg_vec.at(i_pdg)][i_hist]->GetMaximum();
    }
    h_muon_id[pdg_vec.at(0)][i_hist]->SetTitle(("Identified as muon-like;" + bin_info[i_hist].Title + "; Events").c_str());
    h_muon_id[pdg_vec.at(0)][i_hist]->GetYaxis()->SetRangeUser(0, max*1.2);
    if(bin_info[i_hist].FillDataType.compare("p_chi2") == 0) {
      leg_2->Draw();
    } else {
      leg->Draw();
    }
    //ci->SaveAs((folder + "muon/ "+ bin_info.at(i_hist).FillDataType +  ".pdf").c_str());
    //ci->Close();
    //delete ci;
  }


}