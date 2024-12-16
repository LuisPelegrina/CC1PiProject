#include "../../Includes.h"
#include "../tree_utils.cpp"
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

bool Normalize = false;
void background_composition_quick()
{
    
    Cut_Parameters cut_p;

    cut_p.min_distance_to_wall_x_y = 20;
    cut_p.min_distance_to_last_z_wall = 20;
    cut_p.min_distance_to_first_z_wall = 30;
    cut_p.min_distance_to_CPA = 5;

    cut_p.min_crumbs_score = 0;

    cut_p.min_track_lenght = 10;
    cut_p.min_track_score = 0.5;

    Slice *slice = 0;

    const int num_trees = 5;

    TTree *tree[num_trees];
    TFile *input[num_trees];
    string strRuta;

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


    const int num_cuts = 6;
    string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "fv_cut", "crumbs_cut", "track_cut", "razzled_muon_like_cut"};
    //VERTEX CUT AND HIT CT


    vector<string> final_states_particle = {"no_cut", "out_av_nu","cosmic", "CC1Pi", "NC", "CC_e", "CC_mu_0pi_0p", "CC_mu_0pi_1p", "CC_mu_0pi_2p", "CC_mu_2pi"};
    vector<string> final_states_particle_title = { "no_cut", "Out AV #nu", "Cosmic","Signal", "NC", "CC1e", "CC1#mu0#pi0p",  "CC1#mu0#pi1p",  "CC1#mu0#pi2p", "CC1#mu2#pi"};

    //, "NC", "CCe", "CC#mu0#pi0p", "CC#mu0#pi1p", "CC#mu0#pi>2p", "CC#mu>2#pi"};

    const int num_hist = 3;
    struct LocalBinInformation bin_info[num_hist] = {
	    {"E_nu", "E_{#nu} [GeV]", 35, 0, 7},
        {"crumbs_score", "CRUMBS score", 60, 0, 1.2},
        {"E_mu_reco", "E_{#mu} [GeV]", 50, 0, 1.5}
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

            double w = weight;

            bool pass_cut = true;
            for (int i_pcut = 0; i_pcut < num_cuts; i_pcut++) {
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
            Data["E_nu"] = E_0;
            Data["crumbs_score"] = crumbs_score;
            Data["E_mu_reco"] = E_mu_reco/1000;

            for (int i_hist = 0; i_hist < num_hist; i_hist++) {
                h[i_hist][0]->Fill(Data[bin_info[i_hist].FillDataType], w);
            }
            cout << gen_index << endl;

            if (gen_index == 2) {
                for (int i_hist = 0; i_hist < num_hist; i_hist++)
                    h[i_hist][2]->Fill(Data[bin_info[i_hist].FillDataType], w);
            } else if (gen_index == 1) {

                if (is_inside_AV(v_x_true, v_y_true, v_z_true)) {
                    if(final_state->compare("CC1Pi")== 0) {
                        for(int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][3]->Fill(Data[bin_info[i_hist].FillDataType], w);
                    } else if(final_state->compare("NC") == 0) {
                        for(int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][4]->Fill(Data[bin_info[i_hist].FillDataType], w);
                    } else if(final_state->compare("CC_e") == 0) {
                        for(int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][5]->Fill(Data[bin_info[i_hist].FillDataType], w);
                    } else if(final_state->compare("CC_mu_0pi_0p") == 0) {
                        for(int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][6]->Fill(Data[bin_info[i_hist].FillDataType], w);
                    } else if(final_state->compare("CC_mu_0pi_1p") == 0) {
                        for(int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][7]->Fill(Data[bin_info[i_hist].FillDataType], w);
                    } else if(final_state->compare("CC_mu_0pi_2p") == 0) {
                        for(int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][8]->Fill(Data[bin_info[i_hist].FillDataType], w);
                    } else if(final_state->compare("CC_mu_2pi") == 0) {
                        for(int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][9]->Fill(Data[bin_info[i_hist].FillDataType], w);
                    }
                } else {
                    for (int i_hist = 0; i_hist < num_hist; i_hist++) h[i_hist][1]->Fill(Data[bin_info[i_hist].FillDataType], w);
                }
            }

        }
    }


    cout << endl;
    double cont = 0;
    for(int i_fs =0; i_fs < num_final_states; i_fs++) {

        double num_signal = h[1][i_fs]->Integral();
        cout << i_fs << " " << final_states_particle[i_fs] << " events: " << num_signal << " percent: " << num_signal/h[1][0]->Integral() << " percent back: " << num_signal/(h[1][0]->Integral() - h[1][3]->Integral()) << endl;

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
                string current_background =final_states_particle.at(i_fs);

                if(current_background.compare("out_av_nu") == 0) {
                    h[i_hist][i_fs]->SetFillColorAlpha(kBlue, 0.1);
                    h[i_hist][i_fs]->SetLineColor(kBlue+2 );
                } else if(current_background.compare("cosmic") == 0) {
                    h[i_hist][i_fs]->SetFillColorAlpha(kRed, 0.1);
                    h[i_hist][i_fs]->SetLineColor(kRed + 2);
                } else if(current_background.compare("CC1Pi") == 0) {
                    h[i_hist][i_fs]->SetFillColorAlpha(kBlack, 0.1);
                    h[i_hist][i_fs]->SetLineColor(kBlack);
                } else if(current_background.compare("NC") == 0) {
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
            hs->Add(h[i_hist][1]);
            hs->Add(h[i_hist][4]);
            hs->Add(h[i_hist][2]);
            hs->Add(h[i_hist][6]);
            hs->Add(h[i_hist][5]);
            hs->Add(h[i_hist][3]);


            cout << "NICE 2" << endl;
            hs->SetTitle( (";" + bin_info[i_hist].Title + "; Events").c_str());
           
            hs->Draw("hist");

            cout << "NICE 3" << endl;
            TLegend *leg = new TLegend(0.6, 0.4, .88, .85);

            leg->AddEntry(h[i_hist][3], final_states_particle_title[3].c_str(),"lf");
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
            ci->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/background_composition/" + bin_info[i_hist].FillDataType + ".pdf").c_str());
            ci->Close();


        }        



}