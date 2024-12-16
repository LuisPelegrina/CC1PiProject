#include "../../Includes.h"

bool Normalize = false;
void mupi_differentiation()
{
    GenerateDictionaries();
    TTree *tree;
    TFile *input;
    
    Cut_Parameters cut_p;
    Slice *slice = 0;
 
    string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Data/processed_data/processed_data_nu_single.root";
    input = new TFile(strRuta.c_str());
    tree =(TTree*)input->Get("tree");

    tree->SetBranchAddress("slice", &slice);
    int nEntries = tree->GetEntries();

    double cont_fitted = 0;
    double cont_muon_fitted = 0;
    double cont_pion_fitted = 0;
    double cont_total = 0;

    TH2 *h_confusion = new TH2D("h_confusion", "h_confusion",6, 0, 6, 6, 0, 6);
    TH2 *h_confusion_trick = new TH2D("h_confusion_trick", "h_confusion_trick", 6, 0, 6, 6, 0, 6);

    map<int, int> map_pdg_position = {{0, 0}, {11, 1}, {13, 2}, {22, 3}, {211, 4}, {2212, 5}};

    //const int num_contaiment = 11;
    //string s_containment[num_contaiment] = {"no_cut", "both_contained_fv", "none_contained_fv", "only_one_contained_fv", "muon_contained_fv", "pion_contained_fv", "only_muon_contained_fv", "only_pion_contained_fv", "both_candidates_contained_fv", "one_candidate_contained_fv", "no_candidates_contained_fv"};

    const int num_contaiment = 4;
    string s_containment[num_contaiment] = {"no_cut", "both_candidates_contained_fv", "one_candidate_contained_fv", "no_candidates_contained_fv"};

    TH2 *h_confusion_vec[num_contaiment];
    TH2 *h_confusion_trick_vec[num_contaiment];

    const int num_cuts = 8;
    //string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut", "crumbs_cut"};
    string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut", "crumbs_cut", "track_cut", "muon_like_cut", "razzled_muon_like_cut"};


    for(int i_c = 0;i_c < num_contaiment ;i_c++) {
        h_confusion_vec[i_c] = new TH2D(("h_confusion" + to_string(i_c)).c_str(), "h_confusion", 6, 0, 6, 6, 0, 6);
        h_confusion_trick_vec[i_c] = new TH2D(("h_confusion_trick" + to_string(i_c)).c_str(), "h_confusion", 6, 0, 6, 6, 0, 6);
    }


    //CHANGE CUTP FOR ANALYSIS
    cut_p.min_distance_to_CPA = 5;
    cut_p.min_distance_to_first_z_wall = 10;
    cut_p.min_distance_to_last_z_wall = 50;
    cut_p.min_distance_to_wall_x_y = 20;

    for(int i_e = 0; i_e < nEntries; ++i_e) {
        tree->GetEntry(i_e);
        if(i_e%100 == 0) cout << "Entry:" << i_e << endl;
        //cout << endl << endl;

        bool pass_cut = true;

        for(int i_pcut = 0; i_pcut < num_cuts; i_pcut++) {
            if(!pass_cut) continue;
            if(!slice->pass_analysis_cut(cut_p, cut_name[i_pcut])) {
                pass_cut = false;
            }
        }
        if(!pass_cut) continue;

        //TTRABAJO SOLO CON CC1Pi
        if(!slice->true_interaction.is_selected_final_state("CC1Pi")) continue;

        double w = 1;
        //slice->print(cut_p, true, false, false, false, false, false);
        if(Normalize) w = slice->weight;
        //IMPORTANTE PARA NO COMERSE MIERDA
        if(slice->slice_ID != 0) continue;
        cont_total += w;


        //HAS MUON AND PION FITTED TO A MUON-LIKE PARTICLE 
        bool is_muon_fitted = false;
        bool is_pion_fitted = false;

        double max_pion_score = 0;
        double max_muon_score = 0;
        int num_primary_muon_like = 0;
        for(int i_p = 0; i_p < slice->pandora_particle.size();i_p++) {
            Reco_Particle particle = slice->pandora_particle.at(i_p);

            if(!particle.is_muon_like_candidate(slice->primary_vertex_reco.vertex_cm ,cut_p) ) continue;


            double pion_score = particle.razzled_score.pion_score;
            if(pion_score > max_pion_score) max_pion_score = pion_score;

            double muon_score = particle.razzled_score.muon_score;
            if(muon_score > max_muon_score) max_muon_score = muon_score;

            num_primary_muon_like++;
            if(abs(particle.matched_pdg) == 13) is_muon_fitted = true;
            if(abs(particle.matched_pdg) == 211) is_pion_fitted = true;
        }

        if(is_muon_fitted) cont_muon_fitted +=w;
        if(is_pion_fitted) cont_pion_fitted +=w;

        if(is_pion_fitted && is_muon_fitted) cont_fitted += w;

        for(int i_p = 0; i_p < slice->pandora_particle.size();i_p++) {
            Reco_Particle particle = slice->pandora_particle.at(i_p);


            //STUDY WITH PURITY AND COMPLETENESS
            if(particle.track_lenght < 3) continue;
            if(!particle.is_pandora_primary) continue;
            if(particle.purity < 0.5) continue;
            if(particle.completeness < 0.5) continue;

            //cout << particle.matched_pdg << "  " << particle.get_razzled_pdg() << endl;
            h_confusion->Fill(map_pdg_position[abs(particle.matched_pdg)] ,map_pdg_position[particle.get_razzled_pdg()], w);

            //if(map_pdg_position[abs(particle.matched_pdg)] == 0) cout << abs(particle.matched_pdg) << endl;
            //if(map_pdg_position[particle.get_razzled_pdg()] == 0) cout << particle.get_razzled_pdg() << endl;
            for(int i_c = 0;i_c < num_contaiment ;i_c++) {
                if(!slice->pass_containment_cut(s_containment[i_c], 10, cut_p)) continue;
                h_confusion_vec[i_c]->Fill(map_pdg_position[abs(particle.matched_pdg)] ,map_pdg_position[particle.get_razzled_pdg()], w);
            }

            if(num_primary_muon_like < 2) continue;

            if(particle.razzled_score.pion_score == max_pion_score) {
                h_confusion_trick->Fill(map_pdg_position[abs(particle.matched_pdg)], map_pdg_position[211]);
                for(int i_c = 0;i_c < num_contaiment ;i_c++) {
                    if(!slice->pass_containment_cut(s_containment[i_c], 5, cut_p)) continue;
                    h_confusion_trick_vec[i_c]->Fill(map_pdg_position[abs(particle.matched_pdg)], map_pdg_position[211]);
                }
            } else if(particle.razzled_score.muon_score == max_muon_score) {
                h_confusion_trick->Fill(map_pdg_position[abs(particle.matched_pdg)], map_pdg_position[13]);
                for(int i_c = 0;i_c < num_contaiment ;i_c++) {
                    if(!slice->pass_containment_cut(s_containment[i_c], 5, cut_p)) continue;
                    h_confusion_trick_vec[i_c]->Fill(map_pdg_position[abs(particle.matched_pdg)], map_pdg_position[13]);
                }
            }
        }
        
    //double muon_like_dazzle_score = reco_particle.dazzle_score.muon_score + reco_particle.dazzle_score.pion_score;
   

    
        
    }
   
    cout << "total: " << cont_total << endl;
    cout << "fitted: " << cont_fitted << " " << cont_fitted*100.0/cont_total << endl;
    cout << "muon fitted: " << cont_muon_fitted<< " " << cont_muon_fitted*100.0/cont_total <<  endl;
    cout << "pion fitted: " << cont_pion_fitted << " " << cont_pion_fitted*100.0/cont_total << endl;

    h_confusion->SetStats(0);

    h_confusion->GetXaxis()->SetBinLabel(1, "others");
    h_confusion->GetXaxis()->SetBinLabel(2, "electron");
    h_confusion->GetXaxis()->SetBinLabel(3, "muon");
    h_confusion->GetXaxis()->SetBinLabel(4, "gamma");
    h_confusion->GetXaxis()->SetBinLabel(5, "pion");
    h_confusion->GetXaxis()->SetBinLabel(6, "proton");
    
    h_confusion->GetYaxis()->SetBinLabel(1, "others");
    h_confusion->GetYaxis()->SetBinLabel(2, "electron");
    h_confusion->GetYaxis()->SetBinLabel(3, "muon");
    h_confusion->GetYaxis()->SetBinLabel(4, "gamma");
    h_confusion->GetYaxis()->SetBinLabel(5, "pion");
    h_confusion->GetYaxis()->SetBinLabel(6, "proton");

    TCanvas *c_confusion_eff = new TCanvas();

    TH2D *h_confusion_eff = (TH2D*)h_confusion->Clone();
    h_confusion_eff->SetTitle("Efficieny; Matched; Razzled; Events");

    h_confusion_eff->Draw("colz");
    for(int i = 0; i< h_confusion->GetNbinsX(); i++) {
        int TotalP = 0;
        for (int j= 0; j< h_confusion->GetNbinsY(); j++) {
          TotalP +=  h_confusion->GetBinContent(h_confusion->GetBin(i +1, j+1));
        }

      for(int j = 0; j< h_confusion->GetNbinsY(); j++) {
        if(int( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))) == 0) continue;

        TString sPN = TString::Format("%8.2f", h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))*100./TotalP);
        sPN += "%";
        TText *tPN = new TText(i+0.5,j+.5,sPN);
        
        string sPN1 = to_string(int( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))));
        TText *tPN1 = new TText(i+0.5, j+0.6,sPN1.c_str());    
        tPN1->SetTextSize(0.04);
        //tPN1->Draw();
    
        if( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))*100./TotalP > 50) {
          tPN1->SetTextColor(kWhite);
          tPN->SetTextColor(kWhite);
        }
        tPN->SetTextSize(0.04);
        tPN->Draw();

        if(h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1)) == 0) continue;
        h_confusion_eff->SetBinContent(h_confusion->GetBin(i+1, j+1), h_confusion->GetBinContent(h_confusion->GetBin(i+1, j+1))*100./TotalP);
        
      }  
    }

    TCanvas *c_confusion_pur = new TCanvas();

    TH2D *h_confusion_pur = (TH2D*)h_confusion->Clone();
    h_confusion_pur->SetTitle("Purity; Matched; Razzled; Events");

    h_confusion_pur->Draw("colz");
    for(int j = 0; j< h_confusion->GetNbinsY(); j++) {
    
      int TotalP = 0;
      for (int i= 0; i< h_confusion->GetNbinsX(); i++) {
        TotalP +=  h_confusion->GetBinContent(h_confusion->GetBin(i +1, j+1));
      }

      for(int i = 0; i < h_confusion->GetNbinsX(); i++) {
        if(int( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))) == 0) continue;

        TString sPN = TString::Format("%8.2f", h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))*100./TotalP);
        sPN += "%";
        TText *tPN = new TText(i+0.5,j+.5,sPN);
        
        string sPN1 = to_string(int( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))));
        TText *tPN1 = new TText(i+0.5, j+0.6,sPN1.c_str());    
        tPN1->SetTextSize(0.04);
        //tPN1->Draw();

        if( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))*100./TotalP > 50) {
          tPN1->SetTextColor(kWhite);
          tPN->SetTextColor(kWhite);
        }
        tPN->SetTextSize(0.04);
        tPN->Draw();

        if(h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1)) == 0) continue;
        h_confusion_pur->SetBinContent(h_confusion->GetBin(i+1, j+1), h_confusion->GetBinContent(h_confusion->GetBin(i+1, j+1))*100./TotalP);
        

        
      }  
    }
    

    /*
    TCanvas *c_confusion_trick = new TCanvas();

    h_confusion_trick->Draw("colz");
    for(int i = 0; i< 5; i++) {
      for(int j = 0; j< 5; j++) {

        //if(int( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))) == 0) continue;

        string sPN1 = to_string(int( h_confusion_trick->GetBinContent( h_confusion_trick->GetBin(i+1, j+1))));
        TText *tPN1 = new TText(i+0.5, j+0.6,sPN1.c_str());

        int TotalP = 0;
        for (int ji = 0; ji< 5; ji++) {
          TotalP +=  h_confusion_trick->GetBinContent(h_confusion_trick->GetBin(i +1, ji+1));
        }
        TString sPN = TString::Format("%8.2f", h_confusion_trick->GetBinContent( h_confusion_trick->GetBin(i+1, j+1))*100./TotalP);
        sPN += "%";
        TText *tPN = new TText(i+0.5,j+.4,sPN);
      
        if ( h_confusion_trick->GetBinContent( h_confusion_trick->GetBin(i+1, j+1))*100./TotalP > 60) {
          tPN1->SetTextColor(kWhite);
          tPN->SetTextColor(kWhite);
        }
        if(h_confusion_trick->GetBinContent( h_confusion_trick->GetBin(i+1, j+1)) == 0) continue;

        tPN1->SetTextSize(0.04);
        tPN->SetTextSize(0.04);
        tPN1->Draw();
        tPN->Draw();
      }  
    }
    
    h_confusion_trick->SetStats(0);
    h_confusion_trick->SetTitle("confusion forcing; Matched; Reco; Events");

    h_confusion_trick->GetXaxis()->SetBinLabel(1, "others");
    h_confusion_trick->GetXaxis()->SetBinLabel(2, "muon");
    h_confusion_trick->GetXaxis()->SetBinLabel(3, "pion");
    h_confusion_trick->GetXaxis()->SetBinLabel(4, "proton");
    h_confusion_trick->GetXaxis()->SetBinLabel(5, "shower");
    
    h_confusion_trick->GetYaxis()->SetBinLabel(1, "others");
    h_confusion_trick->GetYaxis()->SetBinLabel(2, "muon");
    h_confusion_trick->GetYaxis()->SetBinLabel(3, "pion");
    h_confusion_trick->GetYaxis()->SetBinLabel(4, "proton");
    h_confusion_trick->GetYaxis()->SetBinLabel(5, "shower");

    c_confusion_trick->Update();

  */


    gSystem->Exec("rm -fr /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/mu_pi_diff/trick/*");
    gSystem->Exec("rm -fr /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/mu_pi_diff/normal/*");
    gSystem->Exec("mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/mu_pi_diff/trick");
    gSystem->Exec("mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/mu_pi_diff/normal");

    TH2D *h_confusion_trick_eff_vec[num_contaiment];
    TH2D *h_confusion_trick_pur_vec[num_contaiment];
    for(int i_c = 0;i_c < num_contaiment ;i_c++) {
    
        h_confusion_trick_vec[i_c]->SetStats(0);
   
        h_confusion_trick_vec[i_c]->GetXaxis()->SetBinLabel(2, "electron");
        h_confusion_trick_vec[i_c]->GetXaxis()->SetBinLabel(3, "muon");
        h_confusion_trick_vec[i_c]->GetXaxis()->SetBinLabel(4, "gamma");
        h_confusion_trick_vec[i_c]->GetXaxis()->SetBinLabel(5, "pion");
        h_confusion_trick_vec[i_c]->GetXaxis()->SetBinLabel(6, "proton");
    
        h_confusion_trick_vec[i_c]->GetYaxis()->SetBinLabel(1, "others");
        h_confusion_trick_vec[i_c]->GetYaxis()->SetBinLabel(2, "electron");
        h_confusion_trick_vec[i_c]->GetYaxis()->SetBinLabel(3, "muon");
        h_confusion_trick_vec[i_c]->GetYaxis()->SetBinLabel(4, "gamma");
        h_confusion_trick_vec[i_c]->GetYaxis()->SetBinLabel(5, "pion");
        h_confusion_trick_vec[i_c]->GetYaxis()->SetBinLabel(6, "proton");
   
        h_confusion_trick_eff_vec[i_c] = (TH2D*)h_confusion_trick_vec[i_c]->Clone();
        h_confusion_trick_eff_vec[i_c]->SetTitle("Efficieny; Matched; Razzled; Events");

        h_confusion_trick_pur_vec[i_c] = (TH2D*)h_confusion_trick_vec[i_c]->Clone();
        h_confusion_trick_pur_vec[i_c]->SetTitle("Purity; Matched; Razzled; Events");
        
        double total_pur[int(h_confusion_trick_vec[i_c]->GetNbinsY())];
        double total_eff[int(h_confusion_trick_vec[i_c]->GetNbinsX())];

        TCanvas *c_confusion_trick_eff_vec = new TCanvas();
        h_confusion_trick_eff_vec[i_c]->Draw("col");

        TCanvas *c_confusion_trick_pur_vec = new TCanvas();
        h_confusion_trick_pur_vec[i_c]->Draw("col");

        for (int j= 0; j< h_confusion_trick_vec[i_c]->GetNbinsY(); j++) {
          total_pur[j] = 0;
          for(int i = 0; i< h_confusion_trick_vec[i_c]->GetNbinsX(); i++) {
            total_pur[j] += h_confusion_trick_vec[i_c]->GetBinContent(h_confusion_trick_vec[i_c]->GetBin(i +1, j+1));
          }
        }

        for(int i = 0; i< h_confusion_trick_vec[i_c]->GetNbinsX(); i++) {
          total_eff[i] = 0;
          for (int j= 0; j< h_confusion_trick_vec[i_c]->GetNbinsY(); j++) {
            total_eff[i] += h_confusion_trick_vec[i_c]->GetBinContent(h_confusion_trick_vec[i_c]->GetBin(i +1, j+1));
          }
        }

        for(int i = 0; i< h_confusion_trick_vec[i_c]->GetNbinsX(); i++) {
            for(int j = 0; j< h_confusion_trick_vec[i_c]->GetNbinsY(); j++) {
            if(int( h_confusion_trick_vec[i_c]->GetBinContent( h_confusion_trick_vec[i_c]->GetBin(i+1, j+1))) == 0) continue;
            c_confusion_trick_pur_vec->cd();

            TString sPN_pur = TString::Format("%8.2f", h_confusion_trick_vec[i_c]->GetBinContent( h_confusion_trick_vec[i_c]->GetBin(i+1, j+1))*100./total_pur[j] );
            sPN_pur += "%";
            TText *tPN_pur = new TText(i+0.5,j+.5,sPN_pur);
             
            if(h_confusion_trick_vec[i_c]->GetBinContent( h_confusion_trick_vec[i_c]->GetBin(i+1, j+1))*100./total_pur[j]  > 50) {
              tPN_pur->SetTextColor(kWhite);
            }
            tPN_pur->SetTextSize(0.04);
            tPN_pur->Draw();

            h_confusion_trick_pur_vec[i_c]->SetBinContent(h_confusion_trick_vec[i_c]->GetBin(i+1, j+1), h_confusion_trick_vec[i_c]->GetBinContent(h_confusion_trick_vec[i_c]->GetBin(i+1, j+1))*100./total_pur[j]);



            c_confusion_trick_eff_vec->cd();
            if(int( h_confusion_trick_vec[i_c]->GetBinContent( h_confusion_trick_vec[i_c]->GetBin(i+1, j+1))) == 0) continue;

            TString sPN_eff = TString::Format("%8.2f", h_confusion_trick_vec[i_c]->GetBinContent( h_confusion_trick_vec[i_c]->GetBin(i+1, j+1))*100./total_eff[i]);
            sPN_eff += "%";
            TText *tPN_eff = new TText(i+0.5,j+.5,sPN_eff);
             
            if(h_confusion_trick_vec[i_c]->GetBinContent( h_confusion_trick_vec[i_c]->GetBin(i+1, j+1))*100./total_eff[i]  > 50) {
              tPN_eff->SetTextColor(kWhite);
            }
            tPN_eff->SetTextSize(0.04);
            tPN_eff->Draw();

            h_confusion_trick_eff_vec[i_c]->SetBinContent(h_confusion_trick_eff_vec[i_c]->GetBin(i+1, j+1), h_confusion_trick_vec[i_c]->GetBinContent(h_confusion_trick_vec[i_c]->GetBin(i+1, j+1))*100./total_eff[i] );
          }  
        }
        
        c_confusion_trick_eff_vec->Update();
        c_confusion_trick_pur_vec->Update();

        c_confusion_trick_eff_vec->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/mu_pi_diff/trick/" + s_containment[i_c] + "_eff.pdf").c_str());
        c_confusion_trick_pur_vec->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/mu_pi_diff/trick/" + s_containment[i_c] + "_pur.pdf").c_str());
        
        c_confusion_trick_eff_vec->Close();
        c_confusion_trick_pur_vec->Close();
    }



    TH2D *h_confusion_eff_vec[num_contaiment];
    TH2D *h_confusion_pur_vec[num_contaiment];
    for(int i_c = 0;i_c < num_contaiment ;i_c++) {
    
        h_confusion_vec[i_c]->SetStats(0);
   
        h_confusion_vec[i_c]->GetXaxis()->SetBinLabel(2, "electron");
        h_confusion_vec[i_c]->GetXaxis()->SetBinLabel(3, "muon");
        h_confusion_vec[i_c]->GetXaxis()->SetBinLabel(4, "gamma");
        h_confusion_vec[i_c]->GetXaxis()->SetBinLabel(5, "pion");
        h_confusion_vec[i_c]->GetXaxis()->SetBinLabel(6, "proton");
    
        h_confusion_vec[i_c]->GetYaxis()->SetBinLabel(1, "others");
        h_confusion_vec[i_c]->GetYaxis()->SetBinLabel(2, "electron");
        h_confusion_vec[i_c]->GetYaxis()->SetBinLabel(3, "muon");
        h_confusion_vec[i_c]->GetYaxis()->SetBinLabel(4, "gamma");
        h_confusion_vec[i_c]->GetYaxis()->SetBinLabel(5, "pion");
        h_confusion_vec[i_c]->GetYaxis()->SetBinLabel(6, "proton");
   
        h_confusion_eff_vec[i_c] = (TH2D*)h_confusion_vec[i_c]->Clone();
        h_confusion_eff_vec[i_c]->SetTitle("Efficieny; Matched; Razzled; Events");

        h_confusion_pur_vec[i_c] = (TH2D*)h_confusion_vec[i_c]->Clone();
        h_confusion_pur_vec[i_c]->SetTitle("Purity; Matched; Razzled; Events");
        
        double total_pur[int(h_confusion_vec[i_c]->GetNbinsY())];
        double total_eff[int(h_confusion_vec[i_c]->GetNbinsX())];

        TCanvas *c_confusion_eff_vec = new TCanvas();
        h_confusion_eff_vec[i_c]->Draw("col");

        TCanvas *c_confusion_pur_vec = new TCanvas();
        h_confusion_pur_vec[i_c]->Draw("col");

        for (int j= 0; j< h_confusion_vec[i_c]->GetNbinsY(); j++) {
            total_pur[j] = 0;
          for(int i = 0; i< h_confusion_vec[i_c]->GetNbinsX(); i++) {
            total_pur[j] += h_confusion_vec[i_c]->GetBinContent(h_confusion_vec[i_c]->GetBin(i +1, j+1));
          }
        }

        for(int i = 0; i< h_confusion_vec[i_c]->GetNbinsX(); i++) {
          total_eff[i] = 0;
          for (int j= 0; j< h_confusion_vec[i_c]->GetNbinsY(); j++) {
            total_eff[i] += h_confusion_vec[i_c]->GetBinContent(h_confusion_vec[i_c]->GetBin(i +1, j+1));
          }
        }

        for(int i = 0; i< h_confusion_vec[i_c]->GetNbinsX(); i++) {
            for(int j = 0; j< h_confusion_vec[i_c]->GetNbinsY(); j++) {
            if(int( h_confusion_vec[i_c]->GetBinContent( h_confusion_vec[i_c]->GetBin(i+1, j+1))) == 0) continue;
            c_confusion_pur_vec->cd();

            TString sPN_pur = TString::Format("%8.2f", h_confusion_vec[i_c]->GetBinContent( h_confusion_vec[i_c]->GetBin(i+1, j+1))*100./total_pur[j]);
            sPN_pur += "%";
            TText *tPN_pur = new TText(i+0.5,j+.5,sPN_pur);
             
            if(h_confusion_vec[i_c]->GetBinContent( h_confusion_vec[i_c]->GetBin(i+1, j+1))*100./total_pur[j] > 50) {
              tPN_pur->SetTextColor(kWhite);
            }
            tPN_pur->SetTextSize(0.04);
            tPN_pur->Draw();

            h_confusion_pur_vec[i_c]->SetBinContent(h_confusion_pur_vec[i_c]->GetBin(i+1, j+1), h_confusion_vec[i_c]->GetBinContent( h_confusion_vec[i_c]->GetBin(i+1, j+1))*100./total_pur[j]);


            c_confusion_eff_vec->cd();
            if(int( h_confusion_vec[i_c]->GetBinContent( h_confusion_vec[i_c]->GetBin(i+1, j+1))) == 0) continue;

            TString sPN_eff = TString::Format("%8.2f", h_confusion_vec[i_c]->GetBinContent( h_confusion_vec[i_c]->GetBin(i+1, j+1))*100./total_eff[i]);
            sPN_eff += "%";
            TText *tPN_eff = new TText(i+0.5,j+.5,sPN_eff);
             
            if(h_confusion_vec[i_c]->GetBinContent( h_confusion_vec[i_c]->GetBin(i+1, j+1))*100./total_eff[i] > 50) {
              tPN_eff->SetTextColor(kWhite);
            }
            tPN_eff->SetTextSize(0.04);
            tPN_eff->Draw();

            h_confusion_eff_vec[i_c]->SetBinContent(h_confusion_eff_vec[i_c]->GetBin(i+1, j+1), h_confusion_vec[i_c]->GetBinContent(h_confusion_vec[i_c]->GetBin(i+1, j+1))*100./total_eff[i]);
          }  
        } 
        c_confusion_eff_vec->Update();
        c_confusion_pur_vec->Update();

        c_confusion_eff_vec->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/mu_pi_diff/normal/" + s_containment[i_c] + "_eff.pdf").c_str());
        c_confusion_pur_vec->SaveAs(("/Users/luispelegrinagutierrez/Desktop/Doctorado/Graphs/MuPiReco/mu_pi_diff/normal/" + s_containment[i_c] + "_pur.pdf").c_str());
        
        c_confusion_eff_vec->Close();
        c_confusion_pur_vec->Close();
    }

}
