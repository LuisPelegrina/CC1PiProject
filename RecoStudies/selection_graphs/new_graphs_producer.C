#include "../../Includes.h"

bool Normalize = true;
bool save_graphs = true;
bool check_miss_slice = true;
void new_graphs_producer() {
  GenerateDictionaries();

  string strRuta;
  int num_trees = 1;

  TTree *tree[num_trees];
  TFile *input[num_trees];


  //Declare the variables
  Slice *slice = nullptr;

  //VERTEX CUT AND HIT CT
  vector<string> cut_name = {"no_cut", "is_clear_cosmic_cut", "reco_cut" , "fv_cut",  "crumbs_cut", "track_cut",
                             "razzled_muon_like_cut", "shower_cut", "at_most_one_pfp_scapping_cut", "proton_rejection_razzled_muon_known"};

  const int num_cuts = cut_name.size();
  Cut_Parameters cut_p;
  cut_p.use_default_cut_p();

  const int num_final_states = 4;
  string final_states[num_final_states] = {"no_cut", "CC1Pi", "old_CC1Pi", "N_CC1Pi"};

  vector<struct LocalBinInformation> bin_info = {
    {"E_nu", "E_{#nu} [GeV]", 50, 0, 6},
    {"P_mu", "P_{#mu} [GeV]", 50, 0, 2},
    {"P_pi", "P_{#pi} [GeV]", 50, 0, 2},
    //{"num_razzled_muon_like_candidates", "# Primary tracks", 10, 0, 10},
    //{"num_showers", "# Showers", 10, 0, 10}
  };
  const int num_hist = bin_info.size();


  TH1 *h[num_cuts][num_hist][num_final_states];
  for (int i_cut = 0; i_cut < num_cuts; i_cut++) {
    for (int i_hist = 0; i_hist < num_hist; i_hist++) {
      for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
        string title = "hi" + to_string(i_cut) + "j" + to_string(i_hist) + "k" + to_string(i_fs);
        h[i_cut][i_hist][i_fs] = new TH1D(title.c_str(), " ", bin_info[i_hist].nBins, bin_info[i_hist].LowBin,
                                          bin_info[i_hist].UpBin);
      }
    }
  }

  map<string, double> Data;

  for (int i_t = 0; i_t < num_trees; i_t++) {
    if (i_t == 0)
      strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_82p7k.root";
    //if (i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/Data/processed_data/processed_data_nu_single.root";
    input[i_t] = new TFile(strRuta.c_str());
    tree[i_t] = (TTree *) input[i_t]->Get("tree");
    tree[i_t]->SetBranchAddress("slice", &slice);
    int nEntries = tree[i_t]->GetEntries();


    double prev_run_Id;
    double prev_subrun_Id;
    double prev_event_Id;
    for (int i_e = 0; i_e < nEntries; ++i_e) {
      tree[i_t]->GetEntry(i_e);
      if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;
      //if (slice->slice_ID != 0) continue;
      //Double counting is enabled

      double w = 1;
      if (Normalize) w = slice->weight;

      //Check it is not a "phantom slice"

      double run_Id = slice->run_ID;
      double subrun_Id = slice->subrun_ID;
      double event_Id = slice->event_ID;

      bool extra_slice = false;
      bool has_activity_near_true_vertex = false;

      if(check_miss_slice) {
        if((run_Id == prev_run_Id) && (subrun_Id == prev_subrun_Id) && (event_Id == prev_event_Id))extra_slice = true;
        tree[i_t]->GetEntry(i_e+1);
        if((run_Id == slice->run_ID) && (subrun_Id == slice->subrun_ID) && (event_Id == slice->event_ID))extra_slice = true;
        tree[i_t]->GetEntry(i_e);

        if(extra_slice) {
          vector<int> pfp_ID_vector;
          for(Reco_Particle slice_pfp: slice->pandora_particle) {
            pfp_ID_vector.push_back(slice_pfp.ID);
          }

          for(int i_sp = 0; i_sp < slice->space_point_vec.size(); i_sp++) {
            SpacePoint sp = slice->space_point_vec.at(i_sp);
            if(find(pfp_ID_vector.begin(), pfp_ID_vector.end(), sp.associated_pfp_ID) == pfp_ID_vector.end()) continue;
            double delta_x = sp.x - slice->primary_vertex_true.vertex_cm.x();
            double delta_y = sp.y - slice->primary_vertex_true.vertex_cm.y();
            double delta_z = sp.z - slice->primary_vertex_true.vertex_cm.z();

            if(sqrt(pow(delta_x,2) + sqrt(pow(delta_y,2)) +sqrt(pow(delta_z,2))) < 10) has_activity_near_true_vertex = true;
          }
        }
      }
      bool miss_slice = false;
      if(extra_slice && !has_activity_near_true_vertex) {
        miss_slice = true;
      }

      for (int i_cut = 0; i_cut < num_cuts; i_cut++) {
        bool pass_cut = true;

        for (int i_pcut = 0; i_pcut <= i_cut; i_pcut++) {
          if (!pass_cut) continue;
          if (!slice->pass_analysis_cut(cut_p, cut_name[i_pcut])) {
            pass_cut = false;
          }
        }
        if (!pass_cut) continue;

        Data["E_nu"] = slice->true_interaction.E0_nu;
        vector<GeneratorParticle> muon_vec =  slice->true_interaction.get_generator_primary_p(13);
        vector<GeneratorParticle> pion_vec =  slice->true_interaction.get_generator_primary_p(211);
        Data["P_mu"] = -1;
        Data["P_pi"] = -1;
        if(muon_vec.size() == 1) {
          Data["P_mu"] = muon_vec.at(0).P0.Mag();
        }
        if(pion_vec.size() == 1) {
          Data["P_pi"] = pion_vec.at(0).P0.Mag();
        }

        Data["num_razzled_muon_like_candidates"] =slice->get_num_razzled_particle(211, cut_p,true,true,cut_p.apply_quality_cuts, false) + slice->get_num_razzled_particle(13, cut_p,true,true,cut_p.apply_quality_cuts, false);
        if(miss_slice) {
          for (int i_hist = 0; i_hist < num_hist; i_hist++) {
            h[i_cut][i_hist][0]->Fill(Data[bin_info[i_hist].FillDataType], w);
          }
        } else {
          for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
            if (!slice->true_interaction.is_selected_final_state(final_states[i_fs], 0.325)) continue;

            for (int i_hist = 0; i_hist < num_hist; i_hist++) {
              h[i_cut][i_hist][i_fs]->Fill(Data[bin_info[i_hist].FillDataType], w);
            }
          }
        }
      }

      prev_run_Id = slice->run_ID;
      prev_subrun_Id = slice->subrun_ID;
      prev_event_Id = slice->event_ID;
    }
    input[i_t]->Close();
  }

  //Cout the result
  cout << endl;
  double num_events_0[num_final_states];
  for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
    num_events_0[i_fs] = h[3][0][i_fs]->Integral(0, h[3][0][i_fs]->GetNbinsX());
  }

  for (int i_cut = 0; i_cut < num_cuts; i_cut++) {
    cout << endl;
    cout << cut_name[i_cut] << endl;

    for (int i_fs = 1; i_fs < num_final_states; i_fs++) {
      cout << final_states[i_fs] << endl;

      double num_signal = h[i_cut][0][i_fs]->Integral(0, h[i_cut][0][i_fs]->GetNbinsX());
      double num_background = h[i_cut][0][0]->Integral(0, h[i_cut][0][0]->GetNbinsX()) - num_signal;

      cout << "BackGround: " << num_background << " Signal: " << num_signal << endl;
      cout << "Purity: " << 100 * num_signal / (num_background + num_signal) << endl;

      if (i_cut > 0) {
        cout << "Efficiency ->" << " BackGround: " << 100 * num_background / (num_events_0[0] - num_events_0[i_fs])
             << " Signal: " << 100 * num_signal / num_events_0[i_fs] << endl;
        cout << "signal/srqt(background) ->" << num_signal / sqrt(num_background) << endl;
        cout << "EFF * PUR ->" << num_signal / num_events_0[i_fs] * num_signal / (num_background + num_signal) << endl;

        double num_p_signal = h[i_cut - 1][0][i_fs]->Integral(0, h[i_cut - 1][0][i_fs]->GetNbinsX());
        double num_p_background = h[i_cut - 1][0][0]->Integral(0, h[i_cut - 1][0][0]->GetNbinsX()) - num_p_signal;

        cout << "Efficiency (last cut) ->" << " BackGround: " << 100 * num_background / num_p_background << " Signal: "
             << 100 * num_signal / num_p_signal << endl;
      }
    }
  }

  string folder = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/new_selection_graphs/";
  if (save_graphs) {
    gSystem->Exec(("rm -fr "+folder +"*").c_str());
    int normalize_index = 3;
    TH1 *h_0[num_final_states];
    for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
      h_0[i_fs] = (TH1D *) h[normalize_index][0][i_fs]->Clone();
    }
    //Make the bakcground everything but 1
    int signal_index = 1;
    h_0[0]->Add(h_0[signal_index], -1);

    TH1 *h_0_clone[num_final_states];
    for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
      h_0_clone[i_fs] = (TH1D *) h_0[i_fs]->Clone();
      h_0_clone[i_fs]->Divide(h_0[i_fs]);
    }

    for (int i_cut = 0; i_cut < num_cuts; i_cut++) {
      gSystem->Exec(("mkdir " + folder +cut_name[i_cut]).c_str());
      //prepare the background
      for (int i_hist = 0; i_hist < num_hist; i_hist++) {
        TCanvas *ci = new TCanvas();
        h[i_cut][i_hist][0]->SetTitle((";" + bin_info[i_hist].Title + "; Events").c_str());
        double Max = 0;

        //Substract the signal from the background
        h[i_cut][i_hist][0]->Add(h[i_cut][i_hist][1], -1);
        for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
          if (i_fs == 0) h[i_cut][i_hist][i_fs]->Draw("hist");
          if (i_fs != 0) h[i_cut][i_hist][i_fs]->Draw("hist same");

          if (i_fs == 0) {
            h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kBlack + 2, 0.1);
            h[i_cut][i_hist][i_fs]->SetLineColor(kBlack + 2);
          } else if (i_fs == 1) {
            h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kBlue + 2, 0.1);
            h[i_cut][i_hist][i_fs]->SetLineColor(kBlue + 2);
          } else if (i_fs == 2) {
            h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kRed + 2, 0.1);
            h[i_cut][i_hist][i_fs]->SetLineColor(kRed + 2);
          }

          h[i_cut][i_hist][i_fs]->SetStats(0);
          if (Max < h[i_cut][i_hist][i_fs]->GetMaximum()) Max = h[i_cut][i_hist][i_fs]->GetMaximum();
        }
        h[i_cut][i_hist][0]->SetMaximum(1.2 * Max);
        h[i_cut][i_hist][0]->SetMinimum(0);

        TLegend *leg = new TLegend(0.6, 0.65, .85, .85);
        for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
          if (i_fs == 0) {
            leg->AddEntry(h[i_cut][i_hist][i_fs], "background CC1Pi", "lf");
          } else {
            leg->AddEntry(h[i_cut][i_hist][i_fs], final_states[i_fs].c_str(), "lf");
          }
        }
        leg->Draw();
        ci->SaveAs((folder +cut_name[i_cut] + "/" + bin_info[i_hist].FillDataType + ".pdf").c_str());
        ci->Close();

        TH1 *h_clone[num_final_states];
        for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
          h_clone[i_fs] = (TH1F *) h[i_cut][i_hist][i_fs]->Clone();
        }

        for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
          TLegend *leg_2 = new TLegend(0.6, 0.65, .85, .85);
          leg_2->AddEntry(h_0[i_fs], "Original variable", "lf");
          leg_2->AddEntry(h[i_cut][i_hist][i_fs], "Variable after cuts", "lf");


          TCanvas *cj = new TCanvas();
          // Define two pads with unequal space
          TPad *pad1 = new TPad("pad1", "Pad 1", 0.0, 0.3, 1.0, 1.0); // Larger pad (top)
          TPad *pad2 = new TPad("pad2", "Pad 2", 0.0, 0.0, 1.0, 0.3); // Smaller pad (bottom)

          // Set margins for each pad if needed
          pad2->SetTopMargin(0.02);    // Reduces gap between pads
          pad2->SetBottomMargin(0.3);    // Reduces gap between pads
          pad1->SetBottomMargin(0.05);    // Reduces gap between pads

          // Draw pads on the canvas
          pad1->Draw();
          pad2->Draw();

          pad1->cd();

          TH1* h_in= (TH1D *) h[normalize_index][i_hist][i_fs]->Clone();

          //Make the bakcground everything but 1
          TH1 *h_in_clone = (TH1D *) h_in->Clone();
          h_in_clone->Divide(h_in);


          h_in->SetStats(0);
          h_in->SetLineColor(kBlack);
          h_in->SetFillColorAlpha(kBlack,0);
          h[i_cut][i_hist][i_fs]->SetLineColor(kBlue);
          h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kBlue,0);
          h_in->SetTitle((";" + bin_info[i_hist].Title + "; Events").c_str());
          h_in->Draw("hist");
          h[i_cut][i_hist][i_fs]->Draw("hist same");

          leg_2->Draw();
          pad2->cd();
          h_clone[i_fs]->Divide(h_in);
          h_in_clone->SetTitle((";" + bin_info[i_hist].Title + "; % initial events").c_str());
          h_in_clone->SetStats(0);
          h_in_clone->GetXaxis()->SetTitleSize(0.13); // Adjust these values as needed
          h_in_clone->GetXaxis()->SetLabelSize(0.1); // Adjust these values as needed
          h_in_clone->GetYaxis()->SetTitleSize(0.13);
          h_in_clone->GetYaxis()->SetTitleOffset(0.45);
          h_in_clone->GetYaxis()->SetLabelSize(0.1); // Adjust these values as needed

          h_in_clone->SetLineColor(kBlack);
          h_in_clone->SetFillColorAlpha(kBlack,0);
          h_clone[i_fs]->SetLineColor(kBlue);
          h_clone[i_fs]->SetFillColorAlpha(kBlue,0);

          h_in_clone->GetYaxis()->SetRangeUser(0,1.1);
          h_in_clone->Draw("hist");
          h_clone[i_fs]->Draw("hist same");

          cj->SaveAs((folder + cut_name[i_cut] + "/" + bin_info[i_hist].FillDataType + "_"+ final_states[i_fs] + "comparison.pdf").c_str());
          cj->Close();
        }

      }
    }
  }
}