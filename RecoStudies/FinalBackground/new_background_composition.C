#include "../../Includes.h"
#include "../tree_utils.cpp"

bool Normalize = true;
double high_p0_threshold = 0.325;

bool is_inside_AV(double v_x, double v_y, double v_z){
  bool is_inside = true;

  if (v_x > 200) is_inside = false;
  if (v_x < -200) is_inside = false;

  if (v_y > 200) is_inside = false;
  if (v_y < -200) is_inside = false;

  if (v_z > 500) is_inside = false;
  if (v_z < 0) is_inside = false;

  return is_inside;
}

void new_background_composition()
{
  GenerateDictionaries();

    Slice *slice = 0;

    const int num_trees = 1;

    TTree *tree[num_trees];
    TFile *input[num_trees];
    string strRuta;

  for(int i_t = 0; i_t < num_trees; i_t++) {
    if(i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_82p7k.root";
    input[i_t] = new TFile(strRuta.c_str());

    tree[i_t] =(TTree*)input[i_t]->Get("tree");
    tree[i_t]->SetBranchAddress("slice", &slice);
  }

  vector<string> cut_name = {"no_cut", "is_clear_cosmic_cut", "reco_cut" , "fv_cut",  "crumbs_cut", "track_cut",
                             "razzled_muon_like_cut", "shower_cut", "at_most_one_pfp_scapping_cut", "proton_rejection_razzled_muon_known"};

/*
  vector<string> cut_name = {"no_cut", "is_clear_cosmic_cut", "reco_cut" , "fv_cut",  "crumbs_cut", "track_cut",
                             "razzled_muon_like_cut", "shower_cut", "at_most_one_pfp_scapping_cut"};
*/
  const int num_cuts = cut_name.size();
  Cut_Parameters cut_p;
  cut_p.use_default_cut_p();

  vector<string> final_states_particle = {"no_cut", "out_AV_nu", "cosmic","CC1Pi", "CC_mu_0pi_1p", "CC_mu_0pi_2p", "CC_mu_0pi_0p", "CC_mu_2pi", "NC", "other_CC1pi", "CC_e"};
  vector<string> final_states_particle_title = {"no_cut", "out_AV_nu", "cosmic","#nu_{#mu}CC1#pi", "#nu_{#mu}CC0#pi1p", "#nu_{#mu}CC0#pi2^{+}p", "#nu_{#mu}CC0#pi0p", "#nu_{#mu}CC2^{+}#pi", "NC", "Excluded 1#mu1#pi events", "#nu_{e}CC"};

    vector<struct LocalBinInformation> bin_info= {
	    {"E_nu", "E_{#nu} [GeV]", 50, 0, 7},
      {"crumbs_score", "CRUMBS Score", 50, 0, 1},
      {"angle_candidates", "Angle between candidates [Rad]", 50, 0, 4},
  	};
  const int num_hist = bin_info.size();

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
        //num_entries = 1000;

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
            Data["crumbs_score"] = slice->crumbs_score;

            bool first = true;
            Reco_Particle candidate_1;
            Reco_Particle candidate_2;
            for(Reco_Particle reco_particle: slice->pandora_particle) {
              if((!reco_particle.is_track(cut_p))) continue;
              if(!reco_particle.is_primary(slice->primary_vertex_reco.vertex_cm, cut_p)) continue;
              if(!reco_particle.pass_quality_cuts(cut_p, slice->hits)) continue;
              if(((reco_particle.chi2_score.proton_score < cut_p.proton_rejection_min_proton_score) || (reco_particle.track_lenght < cut_p.proton_rejection_min_TL))) continue;
              if(first) {
                candidate_1 = reco_particle;
                first = false;
              } else {
                candidate_2 = reco_particle;
              }

            }
            Data["angle_candidates"] = candidate_1.start_direction.Angle(candidate_2.start_direction);

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
    for (int i_fs = 3; i_fs < num_final_states; i_fs++) {
      //h[i_hist][i_fs]->SetFillColorAlpha(i_fs, 0.1);
      //h[i_hist][i_fs]->SetLineColor(i_fs);
      string current_background = final_states_particle.at(i_fs-3);
      h[i_hist][i_fs]->SetFillColorAlpha(colors_alpha.at(i_fs-3), 0.1);
      h[i_hist][i_fs]->SetLineColor(colors_alpha.at(i_fs-3));

      h[i_hist][i_fs]->SetStats(0);
    }

    cout << "NICE" << endl;

    hs->Add(h[i_hist][4]); //0pi1p
    hs->Add(h[i_hist][5]);//0pi2p
    hs->Add(h[i_hist][6]); //0pi0p
    hs->Add(h[i_hist][7]);//2pi
    hs->Add(h[i_hist][8]);//NC
    hs->Add(h[i_hist][9]);//other
    //hs->Add(h[i_hist][1]); //cosmics

    //hs->Add(h[i_hist][2]); //out_av_nu
    hs->Add(h[i_hist][10]);//cce
    hs->Add(h[i_hist][3]);//CC1π


    cout << "NICE 2" << endl;
    hs->SetTitle( (";" + bin_info[i_hist].Title + "; Events").c_str());

    hs->Draw("hist");

    cout << "NICE 3" << endl;
    TLegend *leg = new TLegend(0.2, 0.5, .4, .85);
    leg->AddEntry(h[i_hist][3], final_states_particle_title[3].c_str(),"lf");
    leg->AddEntry(h[i_hist][4], final_states_particle_title[4].c_str(),"lf");
    leg->AddEntry(h[i_hist][5], final_states_particle_title[5].c_str(),"lf");
    leg->AddEntry(h[i_hist][6], final_states_particle_title[6].c_str(),"lf");
    leg->AddEntry(h[i_hist][7], final_states_particle_title[7].c_str(),"lf");
    leg->AddEntry(h[i_hist][8], final_states_particle_title[8].c_str(),"lf");
    leg->AddEntry(h[i_hist][9], final_states_particle_title[9].c_str(),"lf");
    //leg->AddEntry(h[i_hist][1], final_states_particle_title[1].c_str(),"lf");
    //leg->AddEntry(h[i_hist][2], final_states_particle_title[2].c_str(),"lf");
    leg->AddEntry(h[i_hist][10], final_states_particle_title[10].c_str(),"lf");

    TLegend *leg_2 = new TLegend(0.68, 0.5, .88, .85);
    TIter nextEntry(leg->GetListOfPrimitives());
    TLegendEntry *entry;
    while ((entry = (TLegendEntry*)nextEntry())) {
      leg_2->AddEntry(entry->GetObject(), entry->GetLabel(), entry->GetOption());
    }

    gPad->Modified();
    gPad->Update();
    cout << "NICE 4" << endl;
    if(bin_info[i_hist].FillDataType.compare("angle_candidates") == 0) {
      leg_2->Draw();
    } else {
      leg->Draw();
    }
    ci->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/background_composition/" + bin_info[i_hist].FillDataType +"_" +  cut_name[num_cuts-1] + ".pdf").c_str());
    //ci->Close();


  }




}