#include "../../Includes.h"
#include "Graphs_utils.cpp"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}

bool Normalize = false;

void Extrap_distance_distribution()
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
  GenerateDictionaries();
  TTree *tree;
  TFile *input;



  //Declare the variables
  Slice *slice = 0;

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
  input = new TFile(strRuta.c_str());
  tree =(TTree*)input->Get("tree");
  tree->SetBranchAddress("slice", &slice);
  int nEntries = tree->GetEntries();


  Cut_Parameters cut_p;
  cut_p.min_track_lenght = 0;
  cut_p.min_track_score = 0;
  cut_p.min_distance_to_CPA = 20;
  cut_p.min_distance_to_first_z_wall = 10;
  cut_p.min_distance_to_last_z_wall = 25;
  cut_p.min_distance_to_wall_x_y = 20;

  const int num_cuts = 4;
  string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut"};
   

  TH1 *extrap_vertex_dist_hist = new TH1D("extrap_vertex_dist_hist", "extrap_vertex_dist_hist", 50, 0, 25);
  TH1 *mu_pi_vertex_dist_hist = new TH1D("mu_pi_vertex_dist_hist", "extrap_vertex_dist_hist", 50, 0, 25);
  TH1 *extrap_vertex_dist_hist_wo_prompt = new TH1D("extrap_vertex_dist_hist_wo_prompt", "extrap_vertex_dist_hist", 50, 0, 25);
  TH1 *extrap_vertex_dist_hist_wo_prompt_w_purcomp = new TH1D("extrap_vertex_dist_hist_wo_prompt_w_purcomp", "extrap_vertex_dist_hist", 50, 0, 25);
  TH1 *extrap_vertex_dist_hist_wo_prompt_w_purcomp_08 = new TH1D("extrap_vertex_dist_hist_wo_prompt_w_purcomp_08", "extrap_vertex_dist_hist_wo_prompt_w_purcomp_08", 50, 0, 25);

  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;

    if (slice->slice_ID != 0) continue;
    if (!slice->true_interaction.is_selected_final_state("N_CC1Pi")) continue;

    double w = 1;
    if (Normalize) w = slice->weight;

    bool pass_cut = true;
    //Plot only if correctly reconstructed
    for (int i_cut = 0; i_cut < num_cuts; i_cut++) {
      if (!pass_cut) continue;
      if (!slice->pass_analysis_cut(cut_p, cut_name[i_cut])) {
        pass_cut = false;
      }
    }
    if (!pass_cut) continue;
    if((slice->primary_vertex_true.vertex_cm - slice->primary_vertex_reco.vertex_cm).Mag() > 2) continue;

    G4Particle primary_muon = slice->true_interaction.get_g4_primary_p(13).at(0);
    G4Particle primary_pion = slice->true_interaction.get_g4_primary_p(211).at(0);


    bool has_prompt_neutron = false;
    for(int i_pt = 0; i_pt < slice->true_interaction.g4_particles.size();i_pt++) {
      G4Particle g4_part = slice->true_interaction.g4_particles.at(i_pt);
      if (!g4_part.is_primary()) continue;
      if ((g4_part.PDG == 2112) && (g4_part.TL < 2)) has_prompt_neutron = true;
    }

    bool show = false;
    for (int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
      Reco_Particle pandora_p = slice->pandora_particle.at(i_p);
      if (slice->pandora_particle.at(i_p).is_pandora_primary) {

        for(int i_pt = 0; i_pt < slice->true_interaction.g4_particles.size();i_pt++) {

          G4Particle g4_part = slice->true_interaction.g4_particles.at(i_pt);
          if(pandora_p.true_track_id == g4_part.ID){
            if(!g4_part.is_primary()) {
              extrap_vertex_dist_hist->Fill((pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag(), w);
              if((primary_muon.TL > 2) && (primary_pion.TL > 2) && !has_prompt_neutron) {
                extrap_vertex_dist_hist_wo_prompt->Fill((pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag(), w);
                if((abs(g4_part.PDG) != 13) && (abs(g4_part.PDG) != 211)) {
                  extrap_vertex_dist_hist_wo_prompt_w_purcomp->Fill((pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag(), w);
                  if(pandora_p.purity > 0.8)extrap_vertex_dist_hist_wo_prompt_w_purcomp_08->Fill((pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag(), w);
                  //if((pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag() < 2) show = true;
                }
              }
            }
            if(g4_part.is_primary() && ((abs(g4_part.PDG) == 13)|| (abs(g4_part.PDG) == 211) ) && ((pandora_p.completeness > 0.5) && (pandora_p.purity > 0.5))) mu_pi_vertex_dist_hist->Fill((pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag(), w);
          }

        }

      }
    }
    if(show)  slice->print(cut_p, true, true, true, true,false, true);
  }

  TCanvas *c1 = new TCanvas();
  c1->cd();
  mu_pi_vertex_dist_hist->Draw("hist");
  mu_pi_vertex_dist_hist->SetTitle(";Distance to vertex [cm]; Events");
  mu_pi_vertex_dist_hist->SetStats(0);
  mu_pi_vertex_dist_hist->SetLineColor(kBlue+2);
  //mu_pi_vertex_dist_hist->SetFillColorAlpha(kBlue+2, 0.1);

  extrap_vertex_dist_hist->Draw("hist same");
  extrap_vertex_dist_hist->SetStats(0);
  extrap_vertex_dist_hist->SetLineColor(kRed+2);
  //extrap_vertex_dist_hist->SetFillColorAlpha(kRed+2, 0.1);

  extrap_vertex_dist_hist_wo_prompt->Draw("hist same");
  extrap_vertex_dist_hist_wo_prompt->SetStats(0);
  extrap_vertex_dist_hist_wo_prompt->SetLineColor(kGreen+2);
  //extrap_vertex_dist_hist_wo_prompt->SetFillColorAlpha(kGreen+2, 0.1);

  extrap_vertex_dist_hist_wo_prompt_w_purcomp->Draw("hist same");
  extrap_vertex_dist_hist_wo_prompt_w_purcomp->SetStats(0);
  extrap_vertex_dist_hist_wo_prompt_w_purcomp->SetLineColor(kOrange+2);
  //extrap_vertex_dist_hist_wo_prompt_w_purcomp->SetFillColorAlpha(kOrange+2, 0.1);

  extrap_vertex_dist_hist_wo_prompt_w_purcomp_08->Draw("hist same");
  extrap_vertex_dist_hist_wo_prompt_w_purcomp_08->SetStats(0);
  extrap_vertex_dist_hist_wo_prompt_w_purcomp_08->SetLineColor(kMagenta+2);

  TLegend *leg  = new TLegend(0.35,0.55,0.88,0.88);
  leg->AddEntry(mu_pi_vertex_dist_hist,"primary #mu/#pi","l");
  leg->AddEntry(extrap_vertex_dist_hist, "non g4 primary", "l");
  leg->AddEntry(extrap_vertex_dist_hist_wo_prompt, "non g4 primary, wo prompt", "l");
  leg->AddEntry(extrap_vertex_dist_hist_wo_prompt_w_purcomp, "non g4 primary, wo prompt, wo #mu/#pi", "l");
  leg->AddEntry(extrap_vertex_dist_hist_wo_prompt_w_purcomp_08, "non g4 primary, wo prompt, wo #mu/#pi, pur > 0.8", "l");


  leg->Draw();
  c1->Update();


  mu_pi_vertex_dist_hist->GetXaxis()->SetRangeUser(0, 25);
  mu_pi_vertex_dist_hist->GetYaxis()->SetRangeUser(0, mu_pi_vertex_dist_hist->GetMaximum()*1.2);
  string ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/DistanceToExtraP/distance_not_normalized.pdf";
  c1->SaveAs(ruta.c_str());


  mu_pi_vertex_dist_hist->GetXaxis()->SetRangeUser(0, 8);
  mu_pi_vertex_dist_hist->GetYaxis()->SetRangeUser(0, extrap_vertex_dist_hist->GetMaximum()*1.2);
  ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/DistanceToExtraP/distance_not_normalized_zoom.pdf";
  c1->SaveAs(ruta.c_str());

  mu_pi_vertex_dist_hist->Scale(1./mu_pi_vertex_dist_hist->Integral());
  extrap_vertex_dist_hist->Scale(1./extrap_vertex_dist_hist->Integral());
  extrap_vertex_dist_hist_wo_prompt->Scale(1./extrap_vertex_dist_hist_wo_prompt->Integral());
  extrap_vertex_dist_hist_wo_prompt_w_purcomp->Scale(1./extrap_vertex_dist_hist_wo_prompt_w_purcomp->Integral());
  extrap_vertex_dist_hist_wo_prompt_w_purcomp_08->Scale(1./extrap_vertex_dist_hist_wo_prompt_w_purcomp_08->Integral());
  c1->Update();

  mu_pi_vertex_dist_hist->SetTitle(";Distance to vertex [cm]; AU");
  mu_pi_vertex_dist_hist->GetXaxis()->SetRangeUser(0, 25);
  mu_pi_vertex_dist_hist->GetYaxis()->SetRangeUser(0, mu_pi_vertex_dist_hist->GetMaximum()*1.2);
  ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/DistanceToExtraP/distance_normalized.pdf";
  c1->SaveAs(ruta.c_str());


  mu_pi_vertex_dist_hist->GetXaxis()->SetRangeUser(0, 8);
  mu_pi_vertex_dist_hist->GetYaxis()->SetRangeUser(0, extrap_vertex_dist_hist->GetMaximum()*1.2);
  ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/DistanceToExtraP/distance_normalized_zoom.pdf";
  c1->SaveAs(ruta.c_str());

}