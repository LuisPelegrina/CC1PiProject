#include "../../Includes.h"

bool Normalize = false;
void chi2_to_tl()
{
  GenerateDictionaries();
    TTree *tree;
    TFile *input;

  Cut_Parameters cut_p;
  cut_p.min_distance_to_wall_x_y = 20;
  cut_p.min_distance_to_last_z_wall = 20;
  cut_p.min_distance_to_first_z_wall = 30;
  cut_p.min_distance_to_CPA = 5;
  cut_p.min_crumbs_score = 0;

  cut_p.apply_quality_cuts = true;
  cut_p.quality_cut_min_track_lenght = 3;
  cut_p.quality_cut_min_total_hits = 15;
  cut_p.quality_cut_min_hits_in_two_planes = 5;

  cut_p.max_primary_distance_to_vertex = 10;
  cut_p.min_track_lenght = 3;
  cut_p.min_track_score = 0.45;

  cut_p.chi2_min_TL = 0;
  cut_p.chi2_max_muon_like_score = 60;
  cut_p.chi2_min_proton_score = 85;

    //Declare the variables
    Slice *slice = 0;

    string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
    input = new TFile(strRuta.c_str());
    tree =(TTree*)input->Get("tree");

    tree->SetBranchAddress("slice", &slice);
    int nEntries = tree->GetEntries();


    double tl_upper_bound = 100;
    double tl_lower_bound = 0;
    double tl_num_bins = 50;

    double chi2_upper_bound = 100;
    double chi2_lower_bound = 0;
    double chi2_num_bins = 50;

    TH2 *h_tl_chi2_p[3] = {new TH2D("h11", "muon", tl_num_bins, tl_lower_bound, tl_upper_bound, chi2_num_bins, chi2_lower_bound, chi2_upper_bound)
                           , new TH2D("h22", "pion", tl_num_bins, tl_lower_bound, tl_upper_bound, chi2_num_bins, chi2_lower_bound, chi2_upper_bound)
                            ,new TH2D("h33", "proton", tl_num_bins, tl_lower_bound, tl_upper_bound, chi2_num_bins, chi2_lower_bound, chi2_upper_bound)};
  TH2 *h_tl_chi2_mu[3] = {new TH2D("h1", "muon", tl_num_bins, tl_lower_bound, tl_upper_bound, chi2_num_bins, chi2_lower_bound, chi2_upper_bound)
    , new TH2D("h2", "pion", tl_num_bins, tl_lower_bound, tl_upper_bound, chi2_num_bins, chi2_lower_bound, chi2_upper_bound)
    ,new TH2D("h3", "proton", tl_num_bins, tl_lower_bound, tl_upper_bound, chi2_num_bins, chi2_lower_bound, chi2_upper_bound)};


  map<int, int> pdg_map = {{0, 0}, {13 , 1}, {211, 2}, {2212, 3}};

    for(int i_e = 0; i_e < nEntries; ++i_e) {
        tree->GetEntry(i_e);

        if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

        double w = 1;
        if(Normalize) w = slice->weight;

        for(Reco_Particle particle: slice->pandora_particle ) {

            if(!particle.is_track(cut_p)) continue;
            if(!particle.is_primary(slice->primary_vertex_reco.vertex_cm, cut_p)) continue;
            if(!particle.pass_quality_cuts(cut_p, slice->hits) )continue;

            if(abs(particle.matched_pdg) == 13) {
              h_tl_chi2_p[0]->Fill(particle.track_lenght, particle.chi2_score.proton_score,w);
              h_tl_chi2_mu[0]->Fill(particle.track_lenght, particle.chi2_score.muon_score,w);
            } else if((abs(particle.matched_pdg) == 211)) {
              h_tl_chi2_p[1]->Fill(particle.track_lenght, particle.chi2_score.proton_score,w);
              h_tl_chi2_mu[1]->Fill(particle.track_lenght, particle.chi2_score.muon_score,w);
            } else if((abs(particle.matched_pdg) == 2212)) {
              h_tl_chi2_p[2]->Fill(particle.track_lenght, particle.chi2_score.proton_score,w);
              h_tl_chi2_mu[2]->Fill(particle.track_lenght, particle.chi2_score.muon_score,w);
            }
        }

    }

  gStyle->SetPadRightMargin(0.175);
    for(int i = 0; i< 3; i++) {
      TCanvas *c_i = new TCanvas();
      string particle;
      if(i == 0) particle =  "Muon";
      if(i == 1) particle =  "Pion";
      if(i == 2) particle =  "Proton";
      string title_mu = particle + "; Track Lenght [cm]; muon #chi^{2} score; events";
      string path_mu = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/Chi2_graphs/";
      path_mu = path_mu + "Chi2_tl_muon" + particle + ".pdf";
      h_tl_chi2_mu[i]->SetStats(0);
      h_tl_chi2_mu[i]->SetTitle(title_mu.c_str());
      h_tl_chi2_mu[i]->Draw("colz");
      c_i->SaveAs(path_mu.c_str());

      TCanvas *c_j = new TCanvas();
      string title_p = particle + "; Track Lenght [cm]; proton #chi^{2} score; events";
      string path_p = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/Chi2_graphs/";
      path_p = path_p + "Chi2_tl_p_" + particle + ".pdf";

      h_tl_chi2_p[i]->SetStats(0);
      h_tl_chi2_p[i]->SetTitle(title_p.c_str());
      h_tl_chi2_p[i]->Draw("colz");
      c_j->SaveAs(path_p.c_str());


    }



}