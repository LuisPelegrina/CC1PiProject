
#include "../../Includes.h"
bool save_graphs = true;
bool Normalize = true;


bool is_inside_AV(double v_x, double v_y, double v_z){
    bool is_inside = true;

    if (v_x > 180) is_inside = false;
    if (v_x < -180) is_inside = false;

    if (v_y > 180) is_inside = false;
    if (v_y < -180) is_inside = false;

    if (v_z > 480) is_inside = false;
    if (v_z < 20) is_inside = false;

    return is_inside;
}

void pur_comp_study()
{
  GenerateDictionaries();

    Cut_Parameters cut_p;

    //Declare the variables
    string strRuta;

    const int num_trees = 1;

    TTree *tree[num_trees];
    TFile *input[num_trees];

    Slice* slice = nullptr;

  int num_hits = 20;
    TH2* h_minimum_hits_in_two_planes= new TH2D("h_minimum_hits_in_two_planes","h_minimum_hits_in_two_planes",20, 0, 1,num_hits, 0,num_hits);
    map<string, TH2*> TH2_map;
    std::vector<string> TH2_cases_string = {"TL_hist", "hitsC_hist", "hitsU_hist", "hitsV_hist", "totalhits_hist", "hitsUV_hist", "hits2lower_hist" };

  int num_bins_x = 20;
  int num_bins_y = 20;
  int upper_bin = 20;
  for(int i_s = 0; i_s < TH2_cases_string.size();i_s++) {
    if(TH2_cases_string.at(i_s) == "hitsC_hist") upper_bin = 20;
    if(TH2_cases_string.at(i_s) == "hitsC_hist") num_bins_x = 20;
    TH2_map[("pur_" + TH2_cases_string.at(i_s))] = new TH2D(("pur_"+TH2_cases_string.at(i_s)).c_str(), ("pur_"+TH2_cases_string.at(i_s)).c_str(), num_bins_x, 0, upper_bin,num_bins_y, 0, 1);
    TH2_map[("comp_" + TH2_cases_string.at(i_s))] = new TH2D(("comp_"+TH2_cases_string.at(i_s)).c_str(), ("comp_"+TH2_cases_string.at(i_s)).c_str(), num_bins_x, 0, upper_bin,num_bins_y, 0, 1);
  }

  for(int i_t = 0; i_t < num_trees; i_t++) {
        if (i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
        input[i_t] = new TFile(strRuta.c_str());
        tree[i_t] = (TTree *) input[i_t]->Get("tree");
        tree[i_t]->SetBranchAddress("slice", &slice);
    }

    for(int i_t = 0; i_t < num_trees; i_t++) {
        int num_entries = tree[i_t]->GetEntries();
        for (int i_e = 0; i_e < num_entries; ++i_e) {
          tree[i_t]->GetEntry(i_e);
          if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;

          //if(slice->slice_ID != 0) continue;
          TrueInteraction true_interaction = slice->true_interaction;

          double w = 1;
          if (Normalize) w = slice->weight;
          //if(!true_interaction.is_selected_final_state("CC1Pi")) continue;
          if(!is_inside_AV(slice->primary_vertex_true.vertex_cm.X(), slice->primary_vertex_true.vertex_cm.Y(), slice->primary_vertex_true.vertex_cm.Z())) continue;

          for(int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
            Reco_Particle p = slice->pandora_particle.at(i_p);
            if((abs(p.matched_pdg) != 13) &&(abs(p.matched_pdg) != 211)) continue;
            if(!p.is_pandora_primary || !p.has_track || p.track_score < 0.5) continue;

            int has_chi2 = 0;
            if(p.chi2_score.muon_score != -1) has_chi2 = 1;

            int num_hits_U = 0;
            int num_hits_V = 0;
            int num_hits_C = 0;
            int total_hits = 0;
            for(int i_h = 0; i_h < slice->hits.size(); i_h++){
              Hit hit = slice->hits.at(i_h);
              if(hit.associated_pfp_ID == p.ID) {
                if(hit.Plane_ID == 2) num_hits_C++;
                if(hit.Plane_ID == 1) num_hits_U++;
                if(hit.Plane_ID == 0) num_hits_V++;
                if(hit.Plane_ID != -1) total_hits++;
              }
            }
            if(total_hits < 15) continue;

            int planes_ok = 0;
            if(num_hits_U >= 5) planes_ok++;
            if(num_hits_V >= 5) planes_ok++;
            if(num_hits_C >= 5) planes_ok++;
            if(planes_ok < 2) continue;

            //if((num_hits_U < 1) || (num_hits_V < 1) || (num_hits_C < 1)) continue;
            if((num_hits_C < 1)) continue;
            TH2_map["pur_TL_hist"]->Fill(p.track_lenght, p.purity, w);
            TH2_map["pur_hitsC_hist"]->Fill(num_hits_C, p.purity, w);
            TH2_map["pur_hitsU_hist"]->Fill(num_hits_U, p.purity, w);
            TH2_map["pur_hitsV_hist"]->Fill(num_hits_V, p.purity, w);
            TH2_map["pur_totalhits_hist"]->Fill(total_hits, p.purity, w);
            TH2_map["pur_hitsUV_hist"]->Fill(num_hits_V + num_hits_U, p.purity, w);

            TH2_map["comp_TL_hist"]->Fill(p.track_lenght, p.completeness, w);
            TH2_map["comp_hitsC_hist"]->Fill(num_hits_C, p.completeness, w);
            TH2_map["comp_hitsU_hist"]->Fill(num_hits_U, p.completeness, w);
            TH2_map["comp_hitsV_hist"]->Fill(num_hits_V, p.completeness, w);
            TH2_map["comp_totalhits_hist"]->Fill(total_hits, p.completeness, w);
            TH2_map["comp_hitsUV_hist"]->Fill(num_hits_V + num_hits_U, p.completeness, w);

            int first_hits = TMath::Max(TMath::Max(num_hits_C,num_hits_U),TMath::Max(num_hits_C,num_hits_V));
            int second_hits = 0;
            int third_hits = 0;
            if(first_hits == num_hits_C) {
             if(num_hits_U > num_hits_V) {
                second_hits = num_hits_U;
                third_hits = num_hits_V;
             } else {
               second_hits = num_hits_V;
               third_hits = num_hits_U;
             }
            } else if(first_hits == num_hits_U){
              if(num_hits_C > num_hits_V) {
                second_hits = num_hits_C;
                third_hits = num_hits_V;
              } else {
                second_hits = num_hits_V;
                third_hits = num_hits_C;
              }
            } else if(first_hits == num_hits_V){
              if(num_hits_U > num_hits_C) {
                second_hits = num_hits_U;
                third_hits = num_hits_C;
              } else {
                second_hits = num_hits_C;
                third_hits = num_hits_U;
              }
            }

            TH2_map["pur_hits2lower_hist"]->Fill(second_hits+third_hits, p.purity, w);
            TH2_map["comp_hits2lower_hist"]->Fill(second_hits+third_hits, p.completeness, w);
            if(p.track_lenght > 20) continue;

            for(int i_c = 0; i_c < num_hits; i_c++) {
              int planes_ok = 0;
              if(num_hits_U >= i_c) planes_ok++;
              if(num_hits_V >= i_c) planes_ok++;
              if(num_hits_C >= i_c) planes_ok++;
              //cout << planes_ok << " " << p.completeness<< endl;
              if(planes_ok >=2) h_minimum_hits_in_two_planes->Fill(p.completeness, i_c,w);
            }

          }
        }
    }


   for(int i_s = 0; i_s < TH2_cases_string.size();i_s++) {
    TCanvas *ci = new TCanvas();
    if(TH2_cases_string.at(i_s) == "TL_hist")  TH2_map[("pur_" + TH2_cases_string.at(i_s))]->SetTitle(";TL [cm];Purity");
     if(TH2_cases_string.at(i_s) == "hitsC_hist") TH2_map[("pur_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in C plane;Purity");
     if(TH2_cases_string.at(i_s) == "hitsU_hist") TH2_map[("pur_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in C plane;Purity");
     if(TH2_cases_string.at(i_s) == "hitsV_hist") TH2_map[("pur_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in C plane;Purity");
     if(TH2_cases_string.at(i_s) == "totalhits_hist") TH2_map[("pur_" + TH2_cases_string.at(i_s))]->SetTitle(";# total of hits in all planes;Purity");
     if(TH2_cases_string.at(i_s) == "hitsUV_hist") TH2_map[("pur_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in U+V plane;Purity");
     if(TH2_cases_string.at(i_s) == "hits2lower_hist") TH2_map[("pur_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in 2 less populated planes;Purity");

     if(TH2_cases_string.at(i_s) == "TL_hist")  TH2_map[("comp_" + TH2_cases_string.at(i_s))]->SetTitle(";TL [cm]; Completeness");
     if(TH2_cases_string.at(i_s) == "hitsC_hist") TH2_map[("comp_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in C plane;Completeness");
     if(TH2_cases_string.at(i_s) == "hitsU_hist") TH2_map[("comp_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in C plane;Completeness");
     if(TH2_cases_string.at(i_s) == "hitsV_hist") TH2_map[("comp_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in C plane;Completeness");
     if(TH2_cases_string.at(i_s) == "totalhits_hist") TH2_map[("comp_" + TH2_cases_string.at(i_s))]->SetTitle(";# total of hits in all planes;Completeness");
     if(TH2_cases_string.at(i_s) == "hitsUV_hist") TH2_map[("comp_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in U+V plane;Completeness");
     if(TH2_cases_string.at(i_s) == "hits2lower_hist") TH2_map[("comp_" + TH2_cases_string.at(i_s))]->SetTitle(";# of hits in 2 less populated planes;Completeness");
     ci->cd();
    ci->Divide(2, 1);
    ci->cd(1);
    TH2_map[("pur_" + TH2_cases_string.at(i_s))]->SetStats(0);
    TH2_map[("pur_" + TH2_cases_string.at(i_s))]->Draw("colz");
    TH2_map[("pur_" + TH2_cases_string.at(i_s))]->ProfileX()->Draw("same");
    ci->cd(2);
    TH2_map[("comp_" + TH2_cases_string.at(i_s))]->SetStats(0);
    TH2_map[("comp_" + TH2_cases_string.at(i_s))]->Draw("colz");
    TH2_map[("comp_" + TH2_cases_string.at(i_s))]->ProfileX()->Draw("same");

    string ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/PurCompTLHits/" + TH2_cases_string.at(i_s) +"_mupi_total_hits_cut.pdf";
    ci->SaveAs(ruta.c_str());


   }

  TH2F *h_minimum_hits_in_two_planes_norm_rows = (TH2F*)h_minimum_hits_in_two_planes->Clone("h_minimum_hits_in_two_planes_norm");
  TCanvas *c1 = new TCanvas();
  for(int i_by =1; i_by <= h_minimum_hits_in_two_planes->GetNbinsY();i_by++) {
    double integral_in_x = 0;
    double mean_pur = 0;

    for(int i_bx =1; i_bx <= h_minimum_hits_in_two_planes->GetNbinsX();i_bx++) {
      integral_in_x +=  h_minimum_hits_in_two_planes->GetBinContent( h_minimum_hits_in_two_planes->GetBin(i_bx, i_by));
      mean_pur += h_minimum_hits_in_two_planes->GetBinContent( h_minimum_hits_in_two_planes->GetBin(i_bx, i_by))*(i_by/20);
    }
    for(int i_bx =1; i_bx <= h_minimum_hits_in_two_planes->GetNbinsX();i_bx++) {
      double bin_content =  h_minimum_hits_in_two_planes->GetBinContent( h_minimum_hits_in_two_planes->GetBin(i_bx, i_by));
      cout <<  integral_in_x << " " << h_minimum_hits_in_two_planes->GetBinContent( h_minimum_hits_in_two_planes->GetBin(i_bx, i_by))<<  endl;
      h_minimum_hits_in_two_planes_norm_rows->SetBinContent(h_minimum_hits_in_two_planes->GetBin(i_bx, i_by), bin_content/integral_in_x );
    }
  }

  TH2F *h_minimum_hits_in_two_planes_norm_columns = (TH2F*)h_minimum_hits_in_two_planes->Clone("h_minimum_hits_in_two_planes_norm_columns");
  for(int i_bx =1; i_bx <= h_minimum_hits_in_two_planes->GetNbinsY();i_bx++) {
    double integral_in_x = 0;

    for(int i_by =1; i_by <= h_minimum_hits_in_two_planes->GetNbinsX();i_by++) {
     if( h_minimum_hits_in_two_planes->GetBinContent( h_minimum_hits_in_two_planes->GetBin(i_bx, i_by)) > integral_in_x)  integral_in_x =  h_minimum_hits_in_two_planes->GetBinContent( h_minimum_hits_in_two_planes->GetBin(i_bx, i_by));
    }
    for(int i_by =1; i_by <= h_minimum_hits_in_two_planes->GetNbinsX();i_by++) {
      double bin_content =  h_minimum_hits_in_two_planes->GetBinContent( h_minimum_hits_in_two_planes->GetBin(i_bx, i_by));
      cout <<  integral_in_x << " " << h_minimum_hits_in_two_planes->GetBinContent( h_minimum_hits_in_two_planes->GetBin(i_bx, i_by))<<  endl;
      h_minimum_hits_in_two_planes_norm_columns->SetBinContent(h_minimum_hits_in_two_planes->GetBin(i_bx, i_by), bin_content/integral_in_x );
    }
  }
  c1->Divide(2, 1);
  c1->cd(1);
  h_minimum_hits_in_two_planes_norm_rows->SetStats(0);
  h_minimum_hits_in_two_planes_norm_rows->SetTitle("Completeness distribution by rows;Completeness; Minimum hits in 2 planes");
  h_minimum_hits_in_two_planes_norm_rows->Draw("colz");

  c1->cd(2);
  h_minimum_hits_in_two_planes_norm_columns->SetStats(0);
  h_minimum_hits_in_two_planes_norm_columns->SetTitle("Efficiency plot;Completeness; Minimum hits in 2 planes");
  h_minimum_hits_in_two_planes_norm_columns->GetZaxis()->SetRangeUser(0,1);
  h_minimum_hits_in_two_planes_norm_columns->Draw("colz");

  string ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/PurCompTLHits/Minimum_hits_mupi_total_hits_cut.pdf";
  c1->SaveAs(ruta.c_str());



}