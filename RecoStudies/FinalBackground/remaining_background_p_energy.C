#include "../../Includes.h"

bool Normalize = false;
bool save_graphs = true;
double high_p0_threshold = 0.325;

void remaining_background_p_energy()
{
  GenerateDictionaries();

  Cut_Parameters cut_p;
  cut_p.min_distance_to_wall_x_y = 20;
  cut_p.min_distance_to_last_z_wall = 20;
  cut_p.min_distance_to_first_z_wall = 30;
  cut_p.min_distance_to_CPA = 5;
  cut_p.min_crumbs_score = 0;

  cut_p.apply_quality_cuts = false;
  cut_p.max_primary_distance_to_vertex = 15;
  cut_p.min_track_lenght = 10;
  cut_p.min_track_score = 0.5;

  Slice *slice = 0;

  const int num_trees = 1;

  TTree *tree[num_trees];
  TFile *input[num_trees];
  string strRuta;

  for(int i_t = 0; i_t < num_trees; i_t++) {
    if(i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
    //if(i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_nu_cosmics.root";
    //if(i_t == 1) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/Data/processed_data/processed_data_in_time_cosmics.root";
    //if(i_t == 2) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/Data/processed_data/processed_data_nu.root";
    //if(i_t == 3) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/Data/processed_data/processed_data_nu_ncc1pi.root";
    //if(i_t == 4) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/Data/processed_data/processed_data_nu_coh.root";
    input[i_t] = new TFile(strRuta.c_str());

    tree[i_t] =(TTree*)input[i_t]->Get("tree");
    tree[i_t]->SetBranchAddress("slice", &slice);
  }

  const int num_cuts = 7;
  string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut", "crumbs_cut", "track_cut", "razzled_muon_like_cut"};


  vector<string> final_states_particle = { "no_cut", "out_AV_nu", "cosmic", "CC1Pi", "NC", "CC_e", "CC_mu_0pi_0p", "CC_mu_0pi_1p", "CC_mu_0pi_2p", "CC_mu_2pi", "other_CC1pi"};
  vector<string> final_states_particle_title = { "no_cut", "out_AV_nu", "cosmic","CC1#pi", "NC", "CCe", "CC#mu0#pi0p", "CC#mu0#pi1p", "CC#mu0#pi>2p", "CC#mu>2#pi","other CC1#pi"};

  TH1 *h_p_energy =

  const int num_hist = 1;
  struct LocalBinInformation bin_info[num_hist] = {
    {"E_nu", "E_{#nu} [GeV]", 35, 0, 7},
  };

  int num_final_states = final_states_particle.size();
  TH1 *h[num_hist][num_final_states];
  for(int i_hist = 0; i_hist < num_hist; i_hist++) {
    for(int i_fs = 0; i_fs < num_final_states; i_fs++) {
      string title = "hj"+ to_string(i_hist)+"k" + to_string(i_fs);
      h[i_hist][i_fs] = new TH1D(title.c_str()," ", bin_info[i_hist].nBins, bin_info[i_hist].LowBin, bin_info[i_hist].UpBin);
    }
  }

  map<string, double> Data;


  for(int i_t = 0; i_t < num_trees; i_t++) {
    int num_entries = tree[i_t]->GetEntries();


    for (int i_e = 0; i_e < num_entries; ++i_e) {
      tree[i_t]->GetEntry(i_e);
      if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;

      double w = 1;
      if (Normalize) w = slice->weight;

      bool pass_cut = true;
      for (int i_cut = 0; i_cut < num_cuts; i_cut++) {
        if (!pass_cut) continue;
        if (!slice->pass_analysis_cut(cut_p, cut_name[i_cut])) {
          pass_cut = false;
        }
      }
      if (!pass_cut) continue;
      Data["E_nu"] = slice->true_interaction.E0_nu;
      Data["crumbs_score"] = slice->true_interaction.E0_nu;

      for (int i_hist = 0; i_hist < num_hist; i_hist++) {
        h[i_hist][0]->Fill(Data[bin_info[i_hist].FillDataType], w);
      }

      string background_state = slice->true_interaction.get_background_final_state(0.325);
      for (int i_fs = 3; i_fs < num_final_states; i_fs++) {
        string current_final_state = final_states_particle.at(i_fs);
        if(background_state != current_final_state) continue;
        for (int i_hist = 0; i_hist < num_hist; i_hist++) {
          h[i_hist][i_fs]->Fill(Data[bin_info[i_hist].FillDataType], w);
        }
      }

      /*
      for (int i_fs = 3; i_fs < num_final_states; i_fs++) {
        string current_final_state = final_states_particle.at(i_fs);
        if(!slice->true_interaction.is_selected_background_final_state(current_final_state, high_p0_threshold)) continue;
        for (int i_hist = 0; i_hist < num_hist; i_hist++) {
          h[i_hist][i_fs]->Fill(Data[bin_info[i_hist].FillDataType], w);
        }
      }
       */

      for (int i_fs = 3; i_fs < num_final_states; i_fs++) {
        string current_final_state = final_states_particle.at(i_fs);
        if(!slice->true_interaction.is_selected_background_final_state(current_final_state, high_p0_threshold)) continue;
        for (int i_hist = 3; i_hist < num_hist; i_hist++) {
          h[i_hist][i_fs]->Fill(Data[bin_info[i_hist].FillDataType], w);
        }
      }

      /*
      if(!is_inside_AV(slice->true_interaction.primary_vertex.X(), slice->true_interaction.primary_vertex.Y(), slice->true_interaction.primary_vertex.Z())) {
        for (int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][1]->Fill(Data[bin_info[i_hist].FillDataType], w);
      } else {
        for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
          string current_final_state = final_states_particle.at(i_fs);
            if(!slice->true_interaction.is_selected_background_final_state(current_final_state, high_p0_threshold)) continue;
            for (int i_hist = 0; i_hist < num_hist; i_hist++) {
              h[i_hist][i_fs]->Fill(Data[bin_info[i_hist].FillDataType], w);
            }
          }
        }
      */



    }
  }


  cout << endl;
  double cont = 0;
  for(int i_fs =0; i_fs < num_final_states; i_fs++) {

    double num_signal = h[0][i_fs]->Integral();
    cout << final_states_particle[i_fs] << " events: " << num_signal << " percent: " << num_signal/h[0][0]->Integral() << " percent back: " << num_signal/(h[0][0]->Integral() - h[0][1]->Integral()) << endl;

    if(i_fs > 0) cont += num_signal ;
  }
  cout << cont << endl;


  gSystem->Exec("mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/background_composition");
  gSystem->Exec("rm -fr /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/background_composition/*");

  //prepare the background


  for(int i_hist = 0; i_hist < num_hist; i_hist++) {
    TCanvas *ci = new TCanvas();
    double Max = 0;
    auto hs = new THStack("hs","");

    //for (int i_fs = 1; i_fs < num_final_states_particle; i_fs++) {
    for (int i_fs = num_final_states-1; i_fs > 0; i_fs--) {
      h[i_hist][i_fs]->SetFillColorAlpha(i_fs, 0.1);
      h[i_hist][i_fs]->SetLineColor(i_fs);
      string current_background = final_states_particle.at(i_fs);

      if(current_background.compare("out_av_nu") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kBlue, 0.1);
        h[i_hist][i_fs]->SetLineColor(kBlue+2 );
      } else if(current_background.compare("cosmic") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kRed, 0.1);
        h[i_hist][i_fs]->SetLineColor(kRed + 2);
      } else if(current_background.compare("CC1Pi") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kBlack, 0.1);
        h[i_hist][i_fs]->SetLineColor(kBlack);
      } else if(current_background.compare("other_CC1pi") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kGray, 0.1);
        h[i_hist][i_fs]->SetLineColor(kGray);
      }  else if(current_background.compare("NC") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kGreen + 2, 0.1);
        h[i_hist][i_fs]->SetLineColor(kGreen + 2);
      } else if(current_background.compare("CC_e") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kCyan, 0.1);
        h[i_hist][i_fs]->SetLineColor(kCyan);
      } else if(current_background.compare("CC_mu_0pi_0p") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kViolet-5, 0.1);
        h[i_hist][i_fs]->SetLineColor(kViolet-5);
      } else if(current_background.compare("CC_mu_0pi_1p") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kOrange+2 , 0.1);
        h[i_hist][i_fs]->SetLineColor(kOrange+2 );
      } else if(current_background.compare("CC_mu_0pi_2p") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kMagenta +2, 0.1);
        h[i_hist][i_fs]->SetLineColor(kMagenta +2);
      } else if(current_background.compare("CC_mu_2pi") == 0) {
        h[i_hist][i_fs]->SetFillColorAlpha(kOrange, 0.1);
        h[i_hist][i_fs]->SetLineColor(kOrange);
      }

      h[i_hist][i_fs]->SetStats(0);
    }

    cout << "NICE" << endl;

    hs->Add(h[i_hist][7]);
    hs->Add(h[i_hist][8]);
    hs->Add(h[i_hist][9]);
    hs->Add(h[i_hist][10]);
    //hs->Add(h[i_hist][1]);
    hs->Add(h[i_hist][4]);
    //hs->Add(h[i_hist][2]);
    hs->Add(h[i_hist][6]);
    hs->Add(h[i_hist][5]);
    hs->Add(h[i_hist][3]);


    cout << "NICE 2" << endl;
    hs->SetTitle( (";" + bin_info[i_hist].Title + "; Events").c_str());

    hs->Draw("hist");

    cout << "NICE 3" << endl;
    TLegend *leg = new TLegend(0.6, 0.4, .88, .85);

    leg->AddEntry(h[i_hist][3], final_states_particle_title[3].c_str(),"lf");
    leg->AddEntry(h[i_hist][10], final_states_particle_title[10].c_str(),"lf");
    leg->AddEntry(h[i_hist][7], final_states_particle_title[7].c_str(),"lf");
    leg->AddEntry(h[i_hist][8], final_states_particle_title[8].c_str(),"lf");
    leg->AddEntry(h[i_hist][9], final_states_particle_title[9].c_str(),"lf");
    leg->AddEntry(h[i_hist][1], final_states_particle_title[1].c_str(),"lf");
    leg->AddEntry(h[i_hist][4], final_states_particle_title[4].c_str(),"lf");
    leg->AddEntry(h[i_hist][2], final_states_particle_title[2].c_str(),"lf");
    leg->AddEntry(h[i_hist][6], final_states_particle_title[6].c_str(),"lf");
    leg->AddEntry(h[i_hist][5], final_states_particle_title[5].c_str(),"lf");



    cout << "NICE 4" << endl;
    leg->Draw();
    //ci->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/background_composition/" + bin_info[i_hist].FillDataType + ".pdf").c_str());
    //ci->Close();

    TCanvas *c1 = new TCanvas();
    h_p_energy->Draw("hist");

    TCanvas *c2 = new TCanvas();
    h_p_tl->Draw("hist");

    if(slice->true_interaction.is_selected_background_final_state("CC_mu_0pi_0p", 0.325) ||  slice->true_interaction.is_selected_background_final_state("CC_mu_0pi_1p", 0.325) || slice->true_interaction.is_selected_background_final_state("CC_mu_0pi_2p", 0.325)) {
      vector<GeneratorParticle> gen_proton_vec = slice->true_interaction.get_generator_primary_p(2212);
      for(GeneratorParticle proton: gen_proton_vec) {
        h_p_energy->Fill(sqrt(proton.P0*proton.P0 + 0.938272 ));
      }
      vector<G4Particle> g4_proton_vec = slice->true_interaction.get_g4_primary_p(2212);
      for(G4Particle proton: g4_proton_vec) {
        h_p_tl->Fill(proton.TL);
      }

    }
  }








