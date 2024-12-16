#include "../../Includes.h"

bool Normalize = true;
bool save_graphs = false;
bool check_miss_slice = true;
double high_p0_threshold = 0.325;

void old_graphs_producer()
{
    GenerateDictionaries();

    string strRuta;
    int num_trees = 1;

    TTree *tree[num_trees];
    TFile *input[num_trees];
    
    Cut_Parameters cut_p;

    //Declare the variables
    Slice *slice = 0;

    const int num_cuts = 7;
    //VERTEX CUT AND HIT CT
    vector<string> cut_name = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut", "crumbs_cut", "track_cut", "razzled_muon_like_cut"};

    cut_p.min_distance_to_wall_x_y = 20;
    cut_p.min_distance_to_last_z_wall = 20;
    cut_p.min_distance_to_first_z_wall = 30;
    cut_p.min_distance_to_CPA = 5;
    cut_p.min_crumbs_score = 0;

    cut_p.apply_quality_cuts = false;
    cut_p.max_primary_distance_to_vertex = 15;
    cut_p.min_track_lenght = 10;
    cut_p.min_track_score = 0.5;


  cut_p.chi2_min_TL = 0;
  cut_p.chi2_max_muon_score = 60;
  cut_p.chi2_max_muon_like_score = 60;
  cut_p.chi2_min_proton_score = 85;

    const int num_final_states = 4;
    string final_states[num_final_states] = {"no_cut", "CC1Pi", "old_CC1Pi", "N_CC1Pi"};

    const int num_hist = 1;
    struct LocalBinInformation bin_info[num_hist] = {
      {"E_nu", "E_{#nu} [GeV]", 20, 0, 7},

      //{"num_primary_tracks", "# Primary Tracks", 10, 0, 10},
	    //{"num_tracks", "# Tracks", 20, 0, 20},
      //{"num_primary_showers", "# Primary Showers", 10, 0, 10},
	    //{"num_showers", "# Showers", 10, 0, 10},
	    //{"num_primary_muon", "# Primary muons", 10, 0, 10},
	    //{"num_primary_pion", "# Primary pions", 10, 0, 10},
      //{"num_primary_muon_like", "# Primary muon like", 10, 0, 10},
	    //{"num_complete_primary_tracks", "# Complete Primary Tracks", 10, 0, 10},	    
        //{"num_hits", "# Hits", 100, 0, 1000}
        
	     //{"num_primary_dazzle_undef", "# dazzle undef", 10, 0, 10},
	     //{"num_primary_dazzle_proton", "# dazzle proton", 10, 0, 10},
	     //{"num_primary_dazzle_muon", "# dazzle_muon", 10, 0, 10},
	     //{"num_primary_dazzle_pion", "# dazzle_pion", 10, 0, 10},
	     //{"num_primary_dazzle_muon_plus_pion", "# dazzle_muon_plus_pion", 10, 0, 10},
	     //{"X0", "X0 [cm]", 200, -200, 200},
	     //{"Y0", "Y0 [cm]", 200, -200, 200},
	     //{"Z0", "Z0 [cm]", 500, 0, 500}
	    //{"num_primary_dazzle_muon_plus_proton", "# dazzle_muon_plus_proton", 10, 0, 10},
	    //{"num_primary_dazzle_pion_plus_proton", "# dazzle_pion_plus_proton", 10, 0, 10},
	    //{"charge_near_vertex_5", "charge_near_vertex_5 (ADC x TICK)", 100, 0, 400000},
	    //{"charge_near_vertex_10", "charge_near_vertex_10 (ADC x TICK)", 100, 0, 400000},
	    //{"charge_near_vertex_15", "charge_near_vertex_15 (ADC x TICK)", 100, 0, 400000},
	    //{"track_differences", "track_differences [cm]", 100, 0, 250},
        
  	};


    TH1 *h[num_cuts][num_hist][num_final_states];
    for(int i_cut = 0; i_cut < num_cuts; i_cut++) {
        for(int i_hist = 0; i_hist < num_hist; i_hist++) {
            for(int i_fs = 0; i_fs < num_final_states; i_fs++) {
                string title = "hi" + to_string(i_cut) +"j"+ to_string(i_hist)+"k" + to_string(i_fs);
                h[i_cut][i_hist][i_fs] = new TH1D(title.c_str()," ", bin_info[i_hist].nBins, bin_info[i_hist].LowBin, bin_info[i_hist].UpBin);
            }
        }
    }

    map<string, double> Data;

  for (int i_t = 0; i_t < num_trees; i_t++) {
    if (i_t == 0)
      strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
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
    for(int i_fs =0; i_fs < num_final_states; i_fs++) {
        num_events_0[i_fs] = h[3][0][i_fs]->Integral(0,h[3][0][i_fs]->GetNbinsX());
    }

    for(int i_cut = 0; i_cut < num_cuts; i_cut++) {
        cout << endl;
        cout << cut_name[i_cut] << endl;

        for(int i_fs =1; i_fs < num_final_states; i_fs++) {
            cout << final_states[i_fs] << endl;
            
            double num_signal = h[i_cut][0][i_fs]->Integral(0,h[i_cut][0][i_fs]->GetNbinsX());
            double num_background = h[i_cut][0][0]->Integral(0,h[i_cut][0][0]->GetNbinsX()) - num_signal;

            cout << "BackGround: " << num_background << " Signal: " << num_signal << endl;
            cout << "Purity: " << 100*num_signal/(num_background + num_signal) << endl;

            if(i_cut > 0) {
                cout << "Efficiency ->" <<  " BackGround: " << 100*num_background/(num_events_0[0]-num_events_0[i_fs])  << " Signal: " << 100*num_signal/num_events_0[i_fs] << endl;
                cout << "signal/srqt(background) ->" << num_signal/sqrt(num_background)  << endl;
                cout << "EFF * PUR ->" << num_signal/num_events_0[i_fs] * num_signal/(num_background + num_signal)   << endl;

                double num_p_signal = h[i_cut-1][0][i_fs]->Integral(0,h[i_cut-1][0][i_fs]->GetNbinsX());
                double num_p_background = h[i_cut-1][0][0]->Integral(0,h[i_cut-1][0][0]->GetNbinsX()) - num_p_signal; 
            
                cout << "Efficiency (last cut) ->" <<  " BackGround: " << 100*num_background/num_p_background << " Signal: " << 100*num_signal/num_p_signal << endl;
            }
        }
    }

    if(save_graphs) {

    //SAVE THE GRAPHS
    gSystem->Exec("rm -fr /Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/old_selection_graphs/*");

    TH1* h_0[num_final_states];
     for(int i_fs = 0; i_fs < num_final_states; i_fs++) {
        h_0[i_fs] = (TH1D*)h[3][0][i_fs]->Clone();
     }
    
    h_0[0]->Add(h_0[1],-1);

    for(int i_cut = 0; i_cut < num_cuts; i_cut++) {
        gSystem->Exec(("mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/old_selection_graphs/" + cut_name[i_cut]).c_str());
        gSystem->Exec(("mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/old_selection_graphs/" + cut_name[i_cut]+ "/LogScale/").c_str());
       
        //if(i_cut != num_cuts - 1) continue;
        //prepare the background
        for(int i_hist = 0; i_hist < num_hist; i_hist++) {
            TCanvas *ci = new TCanvas();
            h[i_cut][i_hist][0]->SetTitle( (";" + bin_info[i_hist].Title + "; Events").c_str());
            double Max = 0;
            
            h[i_cut][i_hist][0]->Add(h[i_cut][i_hist][1], -1);

            for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
                if(i_fs == 0) h[i_cut][i_hist][i_fs]->Draw("hist");
                if(i_fs != 0) h[i_cut][i_hist][i_fs]->Draw("hist same");
                h[i_cut][i_hist][i_fs]->SetFillColorAlpha(i_fs + 1, 0.1);
                h[i_cut][i_hist][i_fs]->SetLineColor(i_fs + 1);
                
                if(i_fs == 2) {
                    h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kBlue, 0.1);
                    h[i_cut][i_hist][i_fs]->SetLineColor(kBlue);
                }
                    

                h[i_cut][i_hist][i_fs]->SetStats(0);
                 if(Max < h[i_cut][i_hist][i_fs]->GetMaximum()) Max = h[i_cut][i_hist][i_fs]->GetMaximum();
            }
            h[i_cut][i_hist][0]->SetMaximum(1.2*Max);
            h[i_cut][i_hist][0]->SetMinimum(0.1);

            TLegend *leg = new TLegend(0.6, 0.65, .85, .85);
            for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
                if(i_fs == 0) {
                    leg->AddEntry(h[i_cut][i_hist][i_fs], "background CC1Pi","lf");
                } else {
                    leg->AddEntry(h[i_cut][i_hist][i_fs], final_states[i_fs].c_str(),"lf");                
                }
            }
            leg->Draw();
            gPad->SetLogy(0);
            ci->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/old_selection_graphs/" + cut_name[i_cut] + "/" + bin_info[i_hist].FillDataType + ".pdf").c_str());
            gPad->SetLogy();
            ci->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/old_selection_graphs/" + cut_name[i_cut] + "/LogScale/" + bin_info[i_hist].FillDataType + ".pdf").c_str());
            ci->Close();

            if(i_hist == 0) {
                TCanvas *cj = new TCanvas();
                
                h[i_cut][i_hist][0]->SetTitle( (";" + bin_info[i_hist].Title + "; Events/Events_0").c_str());
                for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
                    for(int i_b = 1; i_b < h_0[i_fs]->GetNbinsX();i_b++){
                        if(h_0[i_fs]->GetBinContent(i_b) == 0) {
                            h[i_cut][i_hist][i_fs]->SetBinContent(i_b, 0);
                        } else {    
                            h[i_cut][i_hist][i_fs]->SetBinContent(i_b, h[i_cut][i_hist][i_fs]->GetBinContent(i_b)*1./h_0[i_fs]->GetBinContent(i_b));
                        } 
                    }
                    
                    if(i_fs == 0) h[i_cut][i_hist][i_fs]->Draw("hist");
                    if(i_fs != 0) h[i_cut][i_hist][i_fs]->Draw("hist same");
                    h[i_cut][i_hist][i_fs]->SetFillColorAlpha(i_fs + 1, 0.1);
                    h[i_cut][i_hist][i_fs]->SetLineColor(i_fs + 1);
                    if(i_fs == 2) {
                        h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kBlue, 0.1);
                        h[i_cut][i_hist][i_fs]->SetLineColor(kBlue);
                    }
                    

                    h[i_cut][i_hist][i_fs]->SetStats(0);
                    
                    h[i_cut][i_hist][0]->SetMaximum(1.2);
                    h[i_cut][i_hist][0]->SetMinimum(0);
                }
                TGraph* g1 = new TGraph();
                g1->SetPoint(1,100,1);
                g1->SetPoint(0,-10,1);
                g1->SetLineStyle(2);
                g1->Draw("l");

                leg->Draw();
                gPad->SetLogy(0);
                cj->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/old_selection_graphs/" + cut_name[i_cut] + "/" + bin_info[i_hist].FillDataType + "Comp.pdf").c_str());
                cj->Close();
                
                TCanvas *ck = new TCanvas();
                if(i_cut == 0) continue;
                TH1 *h_last[num_final_states];

                for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
                    h_last[i_fs] = (TH1D*) h[i_cut][i_hist][0]->Clone();

                    for(int i_b = 1; i_b < h_0[i_fs]->GetNbinsX();i_b++){
                        if( h[i_cut-1][i_hist][i_fs]->GetBinContent(i_b) == 0) {
                            h_last[i_fs]->SetBinContent(i_b, 0);
                        } else {    
                            h_last[i_fs]->SetBinContent(i_b, h[i_cut][i_hist][i_fs]->GetBinContent(i_b)*1./h[i_cut-1][i_hist][i_fs]->GetBinContent(i_b));
                        } 
                    }
                    //"EX0" es otra opcion pero no salen bien
                    if(i_fs == 0) h_last[i_fs]->Draw("hist");
                    if(i_fs != 0) h_last[i_fs]->Draw("hist same");
                    h_last[i_fs]->SetFillColorAlpha(i_fs + 1, 0.1);
                    h_last[i_fs]->SetLineColor(i_fs + 1);
                    if(i_fs == 2) {
                        h_last[i_fs]->SetFillColorAlpha(kBlue, 0.1);
                        h_last[i_fs]->SetLineColor(kBlue);
                    }
                    h_last[i_fs]->SetStats(0);
                    
                    h_last[i_fs]->SetMaximum(1.2);
                    h_last[i_fs]->SetMinimum(0.1);
                }
                g1->Draw("l");
                leg->Draw();
                gPad->SetLogy(0);
                ck->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/old_selection_graphs/" + cut_name[i_cut] + "/" + bin_info[i_hist].FillDataType + "CompLast.pdf").c_str());
                ck->Close();

            }      
                    

        }        
    }

    }

}