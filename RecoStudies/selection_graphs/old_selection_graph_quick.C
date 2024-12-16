
#include "../../Includes.h"
#include "../tree_utils.cpp"

bool save_graphs = false;
bool Normalize = true;


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

void graphs_quick()
{
    Cut_Parameters cut_p;

    //Declare the variables
    string strRuta;

    const int num_trees = 5;

    TTree *tree[num_trees];
    TFile *input[num_trees];

    for(int i_t = 0; i_t < num_trees; i_t++) {
        if(i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/MuPiProject/Functions/StandardMuPi/final_data/final_data_nu.root";
        if(i_t == 1) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/MuPiProject/Functions/StandardMuPi/final_data/final_data_nu_coh.root";
        if(i_t == 2) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/MuPiProject/Functions/StandardMuPi/final_data/final_data_nu_ncc1pi.root";
        if(i_t == 3) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/MuPiProject/Functions/StandardMuPi/final_data/final_data_nu_cosmics.root";
        if(i_t == 4) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/MuPiProject/Functions/StandardMuPi/final_data/final_data_in_time_cosmics.root";
        input[i_t] = new TFile(strRuta.c_str());
        tree[i_t] =(TTree*)input[i_t]->Get("tree");
        set_Branch(tree[i_t]);
    }


    const int num_cuts = 7;
    //VERTEX CUT AND HIT CT


    string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "crumbs_cut", "fv_cut", "track_cut", "razzled_muon_like_cut"};
    //string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut", "crumbs_cut", "track_cut", "chi2_muon_like_cut", "razzled_candidates_muon_like_cut"};

    const int num_final_states = 6;
    string final_states[num_final_states] = {"no_cut", "nu_back", "out_av_nu","cosmic", "CC1Pi", "N_CC1Pi"};


  cut_p.min_distance_to_wall_x_y = 20;
  cut_p.min_distance_to_last_z_wall = 20;
  cut_p.min_distance_to_first_z_wall = 30;
  cut_p.min_distance_to_CPA = 5;
  cut_p.min_crumbs_score = 0;
  cut_p.apply_quality_cuts = false;
  cut_p.max_primary_distance_to_vertex = 15;
  cut_p.min_track_lenght = 10;
  cut_p.min_track_score = 0.5;


    const int num_hist = 13;
    struct LocalBinInformation bin_info[num_hist] = {

            {"num_primary_tracks", "# Primary Tracks", 10, 0, 10},
            {"num_primary_showers", "# Primary Showers", 10, 0, 10},


            {"crumbs_score", "Crumbs score", 100, -0.2, 1},

            {"num_razzled_primary_photons", "# Razzled primary photons", 10, 0, 10},
            {"num_razzled_primary_protons", "# Razzled primary protons", 10, 0, 10},
            {"num_razzled_primary_electrons", "# Razzled primary electron", 10, 0, 10},
            {"num_razzled_primary_muon_like", "# Razzled primary muon_like", 10, 0, 10},

            {"num_razzled_primary_pions", "num_razzled_primary_pions",5, 0, 5},
            {"num_razzled_primary_muons", "num_razzled_primary_muons",5, 0, 5},
            {"x", "x",200, -200, 200},
            {"y", "y",200, -200, 200},
            {"z", "z",250, 0, 500},
            {"E_mu_reco", "E_{#mu}",250, 0, 500}
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

    for(int i_t = 0; i_t < num_trees; i_t++) {
        int num_entries = tree[i_t]->GetEntries();
        for (int i_e = 0; i_e < num_entries; ++i_e) {
            tree[i_t]->GetEntry(i_e);
            if (i_e % 1000 == 0) cout << "Entry:" << i_e << endl;
            double w = 1;
            if (Normalize) w = weight;

            for (int i_cut = 0; i_cut < num_cuts; i_cut++) {
                bool pass_cut = true;
                for (int i_pcut = 0; i_pcut <= i_cut; i_pcut++) {
                    if (!pass_cut) continue;

                    if (cut_name[i_pcut].compare("no_cut") == 0) {
                        pass_cut = true;
                    } else if (cut_name[i_pcut].compare("is_clear_cosmic_cut") == 0) {
                        if (is_clear_cosmic) pass_cut = false;

                    } else if (cut_name[i_pcut].compare("reco_cut") == 0) {
                        if (!is_reconstructed) pass_cut = false;

                    } else if (cut_name[i_pcut].compare("fv_cut") == 0) {
                        if (v_x > 200 - cut_p.min_distance_to_wall_x_y) pass_cut = false;
                        if (v_x < -200 + cut_p.min_distance_to_wall_x_y) pass_cut = false;
                        if ((v_x < 0 + cut_p.min_distance_to_CPA) && (v_x > 0 - cut_p.min_distance_to_CPA))
                            pass_cut = false;
                        if (v_y > 200 - cut_p.min_distance_to_wall_x_y) pass_cut = false;
                        if (v_y < -200 + cut_p.min_distance_to_wall_x_y) pass_cut = false;
                        if (v_z > 500 - cut_p.min_distance_to_last_z_wall) pass_cut = false;
                        if (v_z < 0 + cut_p.min_distance_to_first_z_wall) pass_cut = false;


                        } else if (cut_name[i_pcut].compare("crumbs_cut") == 0) {
                        if (crumbs_score < cut_p.min_crumbs_score) pass_cut = false;
                    } else if (cut_name[i_pcut].compare("track_cut") == 0) {
                        if (num_primary_tracks < 2) pass_cut = false;

                    } else if (cut_name[i_pcut].compare("no_razzled_primary_electron_cut") == 0) {
                        if (num_razzled_primary_electrons != 0) pass_cut = false;

                    } else if (cut_name[i_pcut].compare("no_razzled_primary_photon_cut") == 0) {
                        if (num_razzled_primary_photons != 0) pass_cut = false;

                    } else if (cut_name[i_pcut].compare("razzled_muon_like_cut") == 0) {
                        if (num_razzled_primary_muon_like != 2) pass_cut = false;

                    } else if (cut_name[i_pcut].compare("one_razzled_primary_muon_cut") == 0) {
                        if (num_razzled_primary_muons != 1) pass_cut = false;

                    } else if (cut_name[i_pcut].compare("one_razzled_primary_pion_cut") == 0) {
                        if (num_razzled_primary_pions != 1) pass_cut = false;
                    }else if (cut_name[i_pcut].compare("chi2_muon_like_cut") == 0) {
                        if (num_chi2_primary_muon_like_candidates != 2) pass_cut = false;
                    }  else if (cut_name[i_pcut].compare("razzled_candidates_muon_like_cut") == 0) {
                        if (num_razzled_primary_muon_like_candidates != 2) pass_cut = false;
                    } else {
                        cout << "ERROR" << cut_name[i_pcut] << endl;
                        pass_cut = false;
                    }
                }
                if (!pass_cut) continue;

                Data["x"] = v_x;
                Data["y"] = v_y;
                Data["z"] = v_z;

                Data["num_primary_tracks"] = num_primary_tracks;
                Data["num_primary_showers"] = num_primary_showers;

                Data["crumbs_score"] = crumbs_score;

                Data["num_razzled_primary_photons"] = num_razzled_primary_photons;
                Data["num_razzled_primary_protons"] = num_razzled_primary_protons;
                Data["num_razzled_primary_electrons"] = num_razzled_primary_electrons;
                Data["num_razzled_primary_muon_like"] = num_razzled_primary_muon_like;

                Data["num_razzled_primary_pions"] = num_razzled_primary_pions;
                Data["num_razzled_primary_muons"] = num_razzled_primary_muons;

                //cout << i_e << "Pass data" << endl;

                for (int i_hist = 0; i_hist < num_hist; i_hist++) {
                    h[i_cut][i_hist][0]->Fill(Data[bin_info[i_hist].FillDataType], w);
                }

                if (gen_index == 2) {
                    for (int i_hist = 0; i_hist < num_hist; i_hist++)
                        h[i_cut][i_hist][3]->Fill(Data[bin_info[i_hist].FillDataType], w);
                } else if (gen_index == 1) {

                    if (is_inside_AV(v_x_true, v_y_true, v_z_true)) {
                        //cout << "2" << endl;


                        bool is_CC1Pi = false;
                        bool is_NCC1Pi = false;

                        if (nu_type->compare("CC1Pi") == 0) is_CC1Pi = true;
                        if (nu_type->compare("N_CC1Pi") == 0) is_CC1Pi = true;
                        if (nu_type->compare("N_CC1Pi_Coherent") == 0) is_CC1Pi = true;

                        if (nu_type->compare("N_CC1Pi") == 0) is_NCC1Pi = true;
                        if (nu_type->compare("N_CC1Pi_Coherent") == 0) is_NCC1Pi = true;

                        if (is_CC1Pi) {
                            for (int i_hist = 0; i_hist < num_hist; i_hist++)
                                h[i_cut][i_hist][4]->Fill(Data[bin_info[i_hist].FillDataType], w);
                        }
                        if(is_NCC1Pi){
                            for (int i_hist = 0; i_hist < num_hist; i_hist++)
                                h[i_cut][i_hist][5]->Fill(Data[bin_info[i_hist].FillDataType], w);

                        }
                        if(!is_CC1Pi){
                            for (int i_hist = 0; i_hist < num_hist; i_hist++) h[i_cut][i_hist][1]->Fill(Data[bin_info[i_hist].FillDataType], w);
                        }
                    } else {
                        for (int i_hist = 0; i_hist < num_hist; i_hist++) h[i_cut][i_hist][2]->Fill(Data[bin_info[i_hist].FillDataType], w);
                    }
                }

            }
        }
    }

    //Cout the result
    cout << "Efficiencies" << endl;

    //int comparison_index = 3;
    int comparison_index = 3;
    int interesting_sample_index = 4;

    double num_events_0[num_final_states];
    for(int i_fs =0; i_fs < num_final_states; i_fs++) {
        num_events_0[i_fs] = h[comparison_index][0][i_fs]->Integral(0,h[comparison_index][0][i_fs]->GetNbinsX());
    }

    double pur[num_cuts][num_final_states];
    double eff[num_cuts][num_final_states];
    for(int i_cut = 0; i_cut < num_cuts; i_cut++) {
        cout << endl << endl;
        cout << cut_name[i_cut] << endl;
        eff[i_cut][0] = 0;

        for(int i_fs = 1; i_fs < num_final_states; i_fs++) {
            //if(i_fs != 4) continue;

            cout << endl << final_states[i_fs] << endl;

            double num_signal = h[i_cut][0][i_fs]->Integral(0,h[i_cut][0][i_fs]->GetNbinsX());
            double num_background = h[i_cut][0][0]->Integral(0,h[i_cut][0][0]->GetNbinsX()) - num_signal;

            cout << "BackGround: " << num_background << " Signal: " << num_signal << endl;
            cout << "Purity: " << 100*num_signal/(num_background + num_signal) << endl;
            pur[i_cut][i_fs] = 100*num_signal/(num_background + num_signal);

            if(i_cut > 0) {
                cout << "Efficiency ->" <<  " BackGround: " << 100*num_background/(num_events_0[0]-num_events_0[i_fs])  << " Signal: " << 100*num_signal/num_events_0[i_fs] << endl;
                cout << "signal/srqt(background) ->" << num_signal/sqrt(num_background)  << endl;
                cout << "EFF * PUR ->" << num_signal/num_events_0[i_fs] * num_signal/(num_background + num_signal)   << endl;

                eff[i_cut][i_fs] = 100*num_signal/num_events_0[i_fs];

                double num_p_signal = h[i_cut-1][0][i_fs]->Integral(0,h[i_cut-1][0][i_fs]->GetNbinsX());
                double num_p_background = h[i_cut-1][0][0]->Integral(0,h[i_cut-1][0][0]->GetNbinsX()) - num_p_signal;

                cout << "Efficiency (last cut) ->" <<  " BackGround: " << 100*num_background/num_p_background << " Signal: " << 100*num_signal/num_p_signal << endl;
            }
        }
    }

    if(save_graphs) {

        double scale_factor = 1;
        //SAVE THE GRAPHS
        gSystem->Exec("rm -fr  /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/graphs_cut/");
        gSystem->Exec("mkdir  /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/graphs_cut/");

        for(int i_cut = 0; i_cut < num_cuts; i_cut++) {
            //if(i_cut < 6) continue;
            gSystem->Exec(("mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/graphs_cut/" + cut_name[i_cut]).c_str());
            gSystem->Exec(("mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/graphs_cut/" + cut_name[i_cut]+ "/LogScale/").c_str());

            //if(i_cut != num_cuts - 1) continue;
            //prepare the background
            for(int i_hist = 0; i_hist < num_hist; i_hist++) {
                TCanvas *ci = new TCanvas();
                h[i_cut][i_hist][0]->SetTitle( (";" + bin_info[i_hist].Title + "; Events").c_str());
                double Max = 0;


                h[i_cut][i_hist][0]->Add(h[i_cut][i_hist][interesting_sample_index], -1);

                for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
                    if(i_fs == 4) {
                        h[i_cut][i_hist][i_fs]->Scale(scale_factor);
                    }

                }

                for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
                    if(i_fs == 0) h[i_cut][i_hist][i_fs]->Draw("hist");
                    if(i_fs != 0) h[i_cut][i_hist][i_fs]->Draw("hist same");
                    h[i_cut][i_hist][i_fs]->SetFillColorAlpha(i_fs + 1, 0.1);
                    h[i_cut][i_hist][i_fs]->SetLineColor(i_fs + 1);

                    if(i_fs == 4) {
                        h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kGreen+2, 0.1);
                        h[i_cut][i_hist][i_fs]->SetLineColor(kGreen+2);
                    } else if (i_fs == 2) {
                        h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kOrange+2, 0.1);
                        h[i_cut][i_hist][i_fs]->SetLineColor(kOrange+2);
                    } else if (i_fs == 5) {
                        h[i_cut][i_hist][i_fs]->SetFillColorAlpha(kGreen, 0.1);
                        h[i_cut][i_hist][i_fs]->SetLineColor(kGreen);
                    }


                    h[i_cut][i_hist][i_fs]->SetStats(0);
                    if(Max < h[i_cut][i_hist][i_fs]->GetMaximum()) Max = h[i_cut][i_hist][i_fs]->GetMaximum();
                }
                h[i_cut][i_hist][0]->SetMaximum(1.2*Max);
                h[i_cut][i_hist][0]->SetMinimum(0.1);



                //double num_signal = h[i_cut][i_hist][interesting_sample_index]->Integral(0,h[i_cut][i_hist][interesting_sample_index]->GetNbinsX())*1.0/scale_factor;
                //double num_background = h[i_cut][i_hist][0]->Integral(0,h[i_cut][i_hist][0]->GetNbinsX()) - num_signal;

                //double eff = 100num_signal/num_events_0[interesting_sample_index] * num_signal/(num_background + num_signal);
                //double pur = 100*num_signal/(num_background + num_signal);
                TText *t_pur;
                cout << pur[i_cut][interesting_sample_index] << " " <<eff[i_cut][interesting_sample_index] << endl;
                TString ts_pur = TString::Format("%3.2f",pur[i_cut][interesting_sample_index]);
                TString ts_tot_pur = TString("Purity: ") + ts_pur +TString("%");

                TText *t_eff;
                TString ts_eff = TString::Format("%3.2f",eff[i_cut][interesting_sample_index]);
                TString ts_tot_eff = TString("Efficiency: ") + ts_eff +TString("%");

                TLegend *leg;
                if(bin_info[i_hist].FillDataType.compare("crumbs_score") == 0) {
                    leg = new TLegend(0.2, 0.6, .5, .88);
                    t_eff = new TText(0.261,0.55,ts_tot_eff);
                    t_pur = new TText(0.25,0.58,ts_tot_pur);
                } else if(bin_info[i_hist].FillDataType.compare("dQdx_score") == 0) {
                    leg = new TLegend(0.3, 0.6, .6, .88);
                    t_eff = new TText(0.361,0.55,ts_tot_eff);
                    t_pur = new TText(0.35,0.58,ts_tot_pur);
                } else if(bin_info[i_hist].FillDataType.compare("bucket_time") == 0) {
                    if(i_cut > 7) {
                        leg = new TLegend(0.2, 0.6, .5, .88);
                        t_eff = new TText(0.261,0.55,ts_tot_eff);
                        t_pur = new TText(0.25,0.58,ts_tot_pur);
                    } else {
                        leg = new TLegend(0.2, 0.6, .5, .88);
                        t_eff = new TText(0.262,0.55,ts_tot_eff);
                        t_pur = new TText(0.25,0.58,ts_tot_pur);
                    }
                } else if(bin_info[i_hist].FillDataType.compare("corrected_time") == 0) {
                    leg = new TLegend(0.2, 0.6, .5, .88);
                    t_eff = new TText(0.261,0.55,ts_tot_eff);
                    t_pur = new TText(0.25,0.58,ts_tot_pur);
                } else if(bin_info[i_hist].FillDataType.compare("nu_score") == 0) {
                    leg = new TLegend(0.2, 0.6, .5, .88);
                    t_eff = new TText(0.261,0.55,ts_tot_eff);
                    t_pur = new TText(0.25,0.58,ts_tot_pur);
                } else {
                    leg = new TLegend(0.5, 0.6, .8, .88);
                    t_eff = new TText(0.561,0.55,ts_tot_eff);
                    t_pur = new TText(0.55,0.58,ts_tot_pur);
                };
                t_eff->SetTextSize(0.02);
                t_eff->SetNDC(true);
                t_pur->SetTextSize(0.02);
                t_pur->SetNDC(true);


                for (int i_fs = 0; i_fs < num_final_states; i_fs++) {
                    if(i_fs == 0) {
                        TString s_sample = TString::Format("%8.0f", h[i_cut][i_hist][i_fs]->Integral(0,h[i_cut][i_hist][i_fs]->GetNbinsX()+1));
                        TString s = TString("Total background: ") + s_sample + TString(" Events");;
                        leg->AddEntry(h[i_cut][i_hist][i_fs], s,"lf");
                    } else {
                        TString s_sample;
                        s_sample = TString::Format("%8.2f", h[i_cut][i_hist][i_fs]->Integral(0,h[i_cut][i_hist][i_fs]->GetNbinsX()+1));

                        TString s;
                        if(final_states[i_fs].compare("nu_back") == 0) s = TString("#nu background: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("out_av_nu") == 0) s = TString("Non-AV #nu: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("cosmic") == 0) s = TString("Cosmics: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("CC1Pi") == 0) s = TString("CC1Pi: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("N_CC1Pi") == 0) s = TString("N_CC1Pi: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("CC_mu_coh") == 0) s = TString("CC_mu_coh: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("CC_mu_qe") == 0) s = TString("CC_mu_qe: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("N_CC1Pi_Coherent") == 0) s = TString("N_CC1Pi_Coherent: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("CC_mu") == 0) s = TString("CC_mu: ") + s_sample + TString(" Events");
                        if(final_states[i_fs].compare("NC") == 0) s = TString("NC: ") + s_sample + TString(" Events");
                        leg->AddEntry(h[i_cut][i_hist][i_fs], s,"lf");
                    }
                }
                const int num_final_states = 7;

                t_pur->Draw();
                t_eff->Draw();

                leg->Draw();
                gPad->SetLogy(0);
                ci->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/graphs_cut/" + cut_name[i_cut] + "/" + bin_info[i_hist].FillDataType + ".png").c_str());
                gPad->SetLogy();
                ci->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/graphs_cut/" + cut_name[i_cut] + "/LogScale/" + bin_info[i_hist].FillDataType + ".png").c_str());
                ci->Close();
                delete ci;
                delete leg;
                //for (int i_fs = 0; i_fs < num_final_states; i_fs++) delete h[i_cut][i_hist][i_fs];
            }
        }

    }

    TCanvas *c12313 = new TCanvas();

    for(int i_t = 0; i_t < num_trees; i_t++) input[i_t]->Close();
}