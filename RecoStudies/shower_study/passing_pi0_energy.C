#include "../../Includes.h"

struct OptimizationParameters {
    string op_name = "ThetaMu";
    double delta = 1;
    double inf_limit = 0;
    double sup_limit = 1;
};

void make_index(vector<vector<int>> &index, int num_op, int iter, int n[], int i[]){
    if(iter == 0) {
        for(i[iter] = 0 ;i[iter] < n[iter] ; i[iter]++) {
        vector<int> v; 
            //for(int j = num_op-1; j > -1 ; j--){
            for(int j = 0; j < num_op ; j++){
                v.push_back(i[j]);
            }
        index.push_back(v);
        }
    } else {
        for(i[iter] = 0 ;i[iter] < n[iter] ; i[iter]++) {
            make_index(index, num_op, iter - 1, n, i);
        }
    }
    
}

int is_contained(G4Particle g4_p) {
  bool is_contained = true;

  if(g4_p.X0.X() > 200) is_contained = false;
  if(g4_p.X0.X() < -200) is_contained = false;
  if(g4_p.X0.Y() > 200) is_contained = false;
  if(g4_p.X0.Y() < -200) is_contained = false;
  if(g4_p.X0.Z() > 500) is_contained = false;
  if(g4_p.X0.Z() < 0) is_contained = false;

  if(g4_p.Xf.X() > 200) is_contained = false;
  if(g4_p.Xf.X() < -200) is_contained = false;
  if(g4_p.Xf.Y() > 200) is_contained = false;
  if(g4_p.Xf.Y() < -200) is_contained = false;
  if(g4_p.Xf.Z() > 500) is_contained = false;
  if(g4_p.Xf.Z() < 0) is_contained = false;

  if(is_contained) return 1;
  if(!is_contained) return 0;

}

bool Normalize = false;
void passing_pi0_energy()
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
  cut_p.chi2_min_proton_score = 80;

    //Declare the variables
    Slice *slice = 0;

    string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
    input = new TFile(strRuta.c_str());
    tree =(TTree*)input->Get("tree");

    tree->SetBranchAddress("slice", &slice);
    int nEntries = tree->GetEntries();
    
    TH1 *h_pion_e = new TH1D("h","h",100,0,5);
  TH1 *h_passing_pion_e = new TH1D("h2","h",100,0,5);
  TH1 *h_passing_pion_cont = new TH1D("h3","h",2,0,2);
  TH1 *h_pion_cont = new TH1D("h3","h",2,0,2);

  //VERTEX CUT AND HIT CT
  vector<string> cut_name = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut", "crumbs_cut", "track_cut",
                             "razzled_muon_like_cut", "shower_cut"};

    for(int i_e = 0; i_e < nEntries; ++i_e) {
        tree->GetEntry(i_e);

        if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

        double w = 1;
        if(Normalize) w = slice->weight;
      if(!slice->true_interaction.is_selected_final_state("old_CC1Pi",0.325)) continue;
      if(slice->true_interaction.is_selected_final_state("CC1Pi",0.325)) continue;

      vector<G4Particle> pi0_vec = slice->true_interaction.get_g4_primary_p(111);
      for(G4Particle pi0: pi0_vec) {
        h_pion_e->Fill(pi0.P0.Mag());
        h_pion_cont->Fill(int(is_contained(pi0)));
      }


      bool pass_cut = true;
        for (int i_cut = 0; i_cut < cut_name.size(); i_cut++) {
          if (!pass_cut) continue;
          if (!slice->pass_analysis_cut(cut_p, cut_name[i_cut])) {
            pass_cut = false;
          }
        }
        if (!pass_cut) continue;

      int num_mu = slice->true_interaction.get_generator_num_primary_p(13);
      int num_pi = slice->true_interaction.get_generator_num_primary_p(211);
      bool is_CC1Pi_old = ((num_pi == 1) && (num_mu == 1));
      for(G4Particle pi0: pi0_vec) {
        h_passing_pion_e->Fill(pi0.P0.Mag());
        h_passing_pion_cont->Fill(int(is_contained(pi0)));
      }

    }

  TCanvas *c1 = new TCanvas();
  h_pion_e->Draw("hist");
    h_passing_pion_e->SetLineColor(kRed +2);
  h_passing_pion_e->Draw("same");

TCanvas *c2 = new TCanvas();
  h_pion_cont->Draw("hist");
  h_pion_cont->GetYaxis()->SetRangeUser(0,2000);
  h_passing_pion_cont->SetLineColor(kRed +2);
  h_passing_pion_cont->Draw("same");

}