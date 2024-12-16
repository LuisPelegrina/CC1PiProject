#include "../../Includes.h"

bool Normalize = false;
void chi2_p_mu_tl_cuts()
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


    int num_tl_cuts = 5;
    double tl_step = 2.5;
    double starting_tl = 2.5;

    double tl_upper_bound = 200;
    double tl_lower_bound = 0;
    double tl_num_bins = 50;

    double chi2_upper_bound = 200;
    double chi2_lower_bound = 0;
    double chi2_num_bins = 50;

    TH2 *h[3][num_tl_cuts];
    for(int i=0; i<3;i++) {
      for(int i_tl = 0; i_tl < num_tl_cuts; i_tl++) {
        string name = "h" + to_string(i)+ "_" + to_string(i_tl);
        h[i][i_tl]  = new TH2D(name.c_str(), name.c_str(), tl_num_bins, tl_lower_bound, tl_upper_bound, chi2_num_bins, chi2_lower_bound, chi2_upper_bound);
      }
    }

  map<double,int> cont_TL_good_mu;
  map<double,int> cont_TL_bad_mu;
  map<double,int> cont_TL_good_p;
  map<double,int> cont_TL_bad_p;
  map<double,int> cont_TL_good_pi;
  map<double,int> cont_TL_bad_pi;


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

          for(int i_tl = 0; i_tl < num_tl_cuts; i_tl++) {
            double TL = starting_tl + tl_step*i_tl;

            if (abs(particle.matched_pdg) == 13) {
              if(particle.track_lenght > TL) h[0][i_tl]->Fill(particle.chi2_score.muon_score, particle.chi2_score.proton_score, w);

              if((particle.track_lenght > TL) && (particle.chi2_score.proton_score > cut_p.chi2_min_proton_score) &&  (particle.chi2_score.muon_score < cut_p.chi2_max_muon_like_score)){
                cont_TL_good_mu[TL]++;
              } else {
                cont_TL_bad_mu[TL]++;
              }

            } else if ((abs(particle.matched_pdg) == 211)) {
              if(particle.track_lenght > TL) h[1][i_tl]->Fill(particle.chi2_score.muon_score, particle.chi2_score.proton_score, w);

              if((particle.track_lenght > TL) && (particle.chi2_score.proton_score > cut_p.chi2_min_proton_score) &&  (particle.chi2_score.muon_score < cut_p.chi2_max_muon_like_score)){
                cont_TL_good_pi[TL]++;
              } else {
                cont_TL_bad_pi[TL]++;
              }

            } else if ((abs(particle.matched_pdg) == 2212)) {
              if(particle.track_lenght > TL) h[2][i_tl]->Fill(particle.chi2_score.muon_score, particle.chi2_score.proton_score, w);

              if((particle.track_lenght > TL) && (particle.chi2_score.proton_score > cut_p.chi2_min_proton_score) &&  (particle.chi2_score.muon_score < cut_p.chi2_max_muon_like_score)){
                cont_TL_bad_p[TL]++;
              } else {
                cont_TL_good_p[TL]++;
              }
            }
          }
        }

    }

  for(int i=0; i<3;i++) {
    if(i == 0) cout << "muon:" << endl;
    if(i == 1) cout << "pion:" << endl;
    if(i == 2) cout << "proton:" << endl;
    for(int i_tl = 0; i_tl < num_tl_cuts; i_tl++) {
      double TL = starting_tl + tl_step*i_tl;
      cout <<"TL: " << TL <<  " good to bad ratio: ";
      if(i == 0) cout << cont_TL_good_mu[TL]*1.0/(cont_TL_good_mu[TL] + cont_TL_bad_mu[TL]);
      if(i == 1) cout << cont_TL_good_pi[TL]*1.0/(cont_TL_good_pi[TL] + cont_TL_bad_pi[TL]);
      if(i == 2) cout << cont_TL_good_p[TL]*1.0/(cont_TL_good_p[TL] + cont_TL_bad_p[TL]);
      cout << endl;
    }

  }
  gStyle->SetPadRightMargin(0.175);
  for(int i=0; i<3;i++) {
    for(int i_tl = 0; i_tl < num_tl_cuts; i_tl++) {
      double TL = starting_tl + tl_step*i_tl;
      TCanvas *c_i = new TCanvas();
      string particle;
      if (i == 0) particle = "Muon";
      if (i == 1) particle = "Pion";
      if (i == 2) particle = "Proton";
      string title_mu = particle + "TL > " + to_string(TL) + "; muon #chi^{2} score; proton #chi^{2} score; events";
      string path_mu = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/Chi2_graphs/TL_comp/";
      path_mu = path_mu + "Chi2_" + particle + "_" + to_string(TL) + ".pdf";
      h[i][i_tl]->SetStats(0);
      h[i][i_tl]->SetTitle(title_mu.c_str());
      h[i][i_tl]->Draw("colz");
      c_i->SaveAs(path_mu.c_str());
      c_i->Close();
    }
  }



}