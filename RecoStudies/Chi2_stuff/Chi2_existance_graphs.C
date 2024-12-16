
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

void Chi2_existance_graphs()
{
  GenerateDictionaries();

    Cut_Parameters cut_p;

    //Declare the variables
    string strRuta;

    const int num_trees = 1;

    TTree *tree[num_trees];
    TFile *input[num_trees];

    Slice* slice = nullptr;

    TH2* chi2_TL_hist = new TH2D("chi2_TL_hist", "chi2_TL_hist", 20, 0, 20,2, 0, 2);
    TH2* chi2_hitsC_hist = new TH2D("chi2_hitsC_hist", "chi2_hitsC_hist", 100, 0, 200,2, 0, 2);
    TH2* chi2_hitsU_hist = new TH2D("chi2_hitsU_hist", "chi2_hitsU_hist", 100, 0, 200,2, 0, 2);
    TH2* chi2_hitsV_hist = new TH2D("chi2_hitsV_hist", "chi2_hitsV_hist", 100, 0, 200,2, 0, 2);

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

          if(slice->slice_ID != 0) continue;
          TrueInteraction true_interaction = slice->true_interaction;
          //if((tm_purity < 0.5)||(tm_completeness < 0.5)) continue;

          double w = 1;
          if (Normalize) w = slice->weight;
          //if(!true_interaction.is_selected_final_state("CC1Pi")) continue;
          if(!is_inside_AV(slice->primary_vertex_true.vertex_cm.X(), slice->primary_vertex_true.vertex_cm.Y(), slice->primary_vertex_true.vertex_cm.Z())) continue;

          for(int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
            Reco_Particle p = slice->pandora_particle.at(i_p);
            if(!p.is_pandora_primary || !p.has_track || p.track_score < 0.5) continue;



            int has_chi2 = 0;
            if(p.chi2_score.muon_score != -1) has_chi2 = 1;

            int num_hits_U = 0;
            int num_hits_V = 0;
            int num_hits_C = 0;
            for(int i_h = 0; i_h < slice->hits.size(); i_h++){
              Hit hit = slice->hits.at(i_h);
              if(hit.associated_pfp_ID == p.ID) {
                if(hit.Plane_ID == 2) num_hits_C++;
                if(hit.Plane_ID == 1) num_hits_U++;
                if(hit.Plane_ID == 0) num_hits_V++;
              }
            }
            chi2_TL_hist->Fill(p.track_lenght, has_chi2, w);
            chi2_hitsC_hist->Fill(num_hits_C, has_chi2, w);
            chi2_hitsU_hist->Fill(num_hits_U, has_chi2, w);
            chi2_hitsV_hist->Fill(num_hits_V, has_chi2, w);

          }
        }
    }

    //Normalize by columns
    for(int i_b = 1; i_b <= chi2_TL_hist->GetNbinsX(); i_b++) {
      double colum_down_value = chi2_TL_hist->GetBinContent(chi2_TL_hist->GetBin(i_b, 1));
      double colum_up_value = chi2_TL_hist->GetBinContent(chi2_TL_hist->GetBin(i_b, 2));
      double colum_total = colum_up_value + colum_down_value;
      chi2_TL_hist->SetBinContent(chi2_TL_hist->GetBin(i_b, 1),colum_down_value/colum_total);
      chi2_TL_hist->SetBinContent(chi2_TL_hist->GetBin(i_b, 2),colum_up_value/colum_total);
    }


  for(int i_b = 1; i_b <= chi2_hitsC_hist->GetNbinsX(); i_b++) {
    double colum_down_value = chi2_hitsC_hist->GetBinContent(chi2_hitsC_hist->GetBin(i_b, 1));
    double colum_up_value = chi2_hitsC_hist->GetBinContent(chi2_hitsC_hist->GetBin(i_b, 2));
    double colum_total = colum_up_value + colum_down_value;
    chi2_hitsC_hist->SetBinContent(chi2_hitsC_hist->GetBin(i_b, 1),colum_down_value/colum_total);
    chi2_hitsC_hist->SetBinContent(chi2_hitsC_hist->GetBin(i_b, 2),colum_up_value/colum_total);
  }

  for(int i_b = 1; i_b <= chi2_hitsU_hist->GetNbinsX(); i_b++) {
    double colum_down_value = chi2_hitsU_hist->GetBinContent(chi2_hitsU_hist->GetBin(i_b, 1));
    double colum_up_value = chi2_hitsU_hist->GetBinContent(chi2_hitsU_hist->GetBin(i_b, 2));
    double colum_total = colum_up_value + colum_down_value;
    chi2_hitsU_hist->SetBinContent(chi2_hitsU_hist->GetBin(i_b, 1),colum_down_value/colum_total);
    chi2_hitsU_hist->SetBinContent(chi2_hitsU_hist->GetBin(i_b, 2),colum_up_value/colum_total);
  }

  for(int i_b = 1; i_b <= chi2_hitsV_hist->GetNbinsX(); i_b++) {
    double colum_down_value = chi2_hitsV_hist->GetBinContent(chi2_hitsV_hist->GetBin(i_b, 1));
    double colum_up_value = chi2_hitsV_hist->GetBinContent(chi2_hitsV_hist->GetBin(i_b, 2));
    double colum_total = colum_up_value + colum_down_value;
    chi2_hitsV_hist->SetBinContent(chi2_hitsV_hist->GetBin(i_b, 1),colum_down_value/colum_total);
    chi2_hitsV_hist->SetBinContent(chi2_hitsV_hist->GetBin(i_b, 2),colum_up_value/colum_total);
  }


  TCanvas *c1 = new TCanvas();
    c1->cd();
    chi2_TL_hist->SetStats(0);
    chi2_TL_hist->Draw("colz");



  TCanvas *c2 = new TCanvas();
  chi2_hitsC_hist->Draw("colz");
  chi2_hitsC_hist->SetStats(0);

  TCanvas *c3 = new TCanvas();
  chi2_hitsU_hist->Draw("colz");
  chi2_hitsU_hist->SetStats(0);

  TCanvas *c4 = new TCanvas();
  chi2_hitsV_hist->Draw("colz");
  chi2_hitsV_hist->SetStats(0);
}