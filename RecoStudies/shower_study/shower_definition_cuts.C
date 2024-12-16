#include "../../Includes.h"

bool discard_two_showers = false;

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

bool Normalize = false;
void shower_definition_cuts()
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
  cut_p.min_track_score = 0.5;

  cut_p.chi2_min_TL = 0;
  cut_p.chi2_max_muon_score = 60;
  cut_p.chi2_max_muon_like_score = 60;
  cut_p.chi2_min_proton_score = 85;


    //Declare the variables
    Slice *slice = 0;

    string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_82p7k.root";
    input = new TFile(strRuta.c_str());
    tree =(TTree*)input->Get("tree");

    tree->SetBranchAddress("slice", &slice);
    int nEntries = tree->GetEntries();
    
    TH2 *h_muonlike_proton[4];
    TH2 *h_muonlike_other[4];
    TH1 *h_tl[2] = {new TH1D("h1", "h", 100, 0, 200), new TH1D("h2", "h", 100, 0, 100)};

    
    for(int i = 0; i < 4; i++) {
        string title = "hi" + to_string(i);
        h_muonlike_proton[i] = new TH2D((title).c_str()," ", 25, 0,100, 50, 0, 400);
        h_muonlike_other[i] = new TH2D((title).c_str()," ", 25, 0, 100, 50, 0, 400);
    }


    map<int, int> pdg_map = {{0, 0}, {13 , 1}, {211, 2}, {2212, 3}};

    vector<struct OptimizationParameters> op_par = {
        {"max_ts", 0.05, 0.3, 0.5},
        {"min_shower_lenght", 5, 20, 110},
        {"max_ts_no_e", 0.05, 0.45, 0.45},
    };
  const int num_op = op_par.size();
    cout << num_op << endl;

    int n[num_op];
    for(int i = 0 ;i < num_op ; i++) {
        n[i] = (op_par[i].sup_limit - op_par[i].inf_limit)/op_par[i].delta + 1;
    }



  int size = 1;
    for(int i_op = 0; i_op < num_op; i_op++) {
        size *= n[i_op];
    }
    cout << size << endl;


    double cont_cc1pi_good[size];
    double cont_cc1pi_bad[size];
    double cont_other_good[size];
    double cont_other_bad[size];

    vector<vector<int>> index;
    int i[num_op];
    
    make_index(index, num_op, num_op - 1, n, i);
    cout << index.size() << endl; 
    cout << "A" << endl;

    for(int i=0; i < size;i++) {
        for(int j = 0; j < num_op ;j++){
           cout  << index.at(i).at(j) << " " ;
        }        
        cout  << " " << i << endl;

      cont_cc1pi_good[i] = 0;
      cont_cc1pi_bad[i] = 0;
      cont_other_good[i] = 0;
      cont_other_bad[i] = 0;
    }




  //VERTEX CUT AND HIT CT
  vector<string> cut_name = {"no_cut", "is_clear_cosmic_cut", "reco_cut", "fv_cut", "crumbs_cut", "track_cut",
                             "razzled_muon_like_cut"};
  double prev_run_Id;
  double prev_subrun_Id;
  double prev_event_Id;

    for(int i_e = 0; i_e < nEntries; ++i_e) {
        tree->GetEntry(i_e);

        if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

        double w = 1;
        if(Normalize) w = slice->weight;

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
      //if(!is_CC1Pi_old) continue;

      for(int i_c = 0; i_c < size; i_c++) {
        double max_ts = op_par[0].inf_limit + index.at(i_c).at(0) * op_par[0].delta;
        double min_ke = op_par[1].inf_limit + index.at(i_c).at(1) * op_par[1].delta;
        double max_ts_no_e = op_par[2].inf_limit + index.at(i_c).at(2) * op_par[2].delta;
        //double min_dEdx =  op_par[2].inf_limit + index.at(i_c).at(2)*op_par[2].delta;

        bool has_shower = false;

        int num_showers = 0;
        for(int i_p = 0; i_p <slice->pandora_particle.size(); i_p++) {
            Reco_Particle particle = slice->pandora_particle.at(i_p);
            if(particle.track_score > 0.5) continue;
              if ((particle.track_score < max_ts)
                  && (particle.shower_length > min_ke)) has_shower = true;
              if(particle.track_score < max_ts_no_e) num_showers++;
            }

        if((num_showers > 1) && discard_two_showers) has_shower = true;

        if(slice->true_interaction.is_selected_final_state("CC1Pi", 0.325)) {
          if(has_shower) {
            cont_cc1pi_bad[i_c]++;
          } else {
            cont_cc1pi_good[i_c]++;
          }
        } else {
          if(has_shower) {
            cont_other_good[i_c]++;
          } else {
            cont_other_bad[i_c]++;
          }
        }


      }


    }


  int TH2_plot_index_1 = 0;
  int TH2_plot_index_2 = 1;
  TH2 *h_eff = new TH2D("h_eff", "h_eff", n[TH2_plot_index_1], op_par[TH2_plot_index_1].inf_limit, op_par[TH2_plot_index_1].sup_limit + op_par[TH2_plot_index_1].delta, n[TH2_plot_index_2], op_par[TH2_plot_index_2].inf_limit, op_par[TH2_plot_index_2].sup_limit + op_par[TH2_plot_index_2].delta);
  TH2 *h_pur = new TH2D("h_pur", "h_pur", n[TH2_plot_index_1], op_par[TH2_plot_index_1].inf_limit, op_par[TH2_plot_index_1].sup_limit + op_par[TH2_plot_index_1].delta, n[TH2_plot_index_2], op_par[TH2_plot_index_2].inf_limit, op_par[TH2_plot_index_2].sup_limit + op_par[TH2_plot_index_2].delta);
  TH2 *h_pureff = new TH2D("h_pureff", "h_pureff", n[TH2_plot_index_1], op_par[TH2_plot_index_1].inf_limit, op_par[TH2_plot_index_1].sup_limit + op_par[TH2_plot_index_1].delta, n[TH2_plot_index_2], op_par[TH2_plot_index_2].inf_limit, op_par[TH2_plot_index_2].sup_limit + op_par[TH2_plot_index_2].delta);

  //SHOW OPTIMIZATION

    double best_pureff = 0;
    double best_pur = 0;
    double best_eff = 0;

    double best_cut[num_op];
    //SHOW BACKGROUND AND SIGNAL 
    for(int i_c = 0; i_c < size; i_c++) {

        double eff = cont_cc1pi_good[i_c]/(cont_cc1pi_good[i_c] + cont_cc1pi_bad[i_c]);
        double pur = cont_cc1pi_good[i_c]/(cont_cc1pi_good[i_c] + cont_other_bad[i_c]);

        cout << cont_cc1pi_good[i_c] + cont_cc1pi_bad[i_c] << " CUT: " ;

        for(int j = 0; j < num_op ;j++){
           cout << op_par[j].inf_limit + index.at(i_c).at(j) * op_par[j].delta << " ";
        }
        cout << "pur: " << pur << " eff: " << eff  << " purr*eff: " << pur*eff <<  endl;
        double par_1 = op_par[TH2_plot_index_1].inf_limit + index.at(i_c).at(TH2_plot_index_1) * op_par[TH2_plot_index_1].delta;
        double par_2 = op_par[TH2_plot_index_2].inf_limit + index.at(i_c).at(TH2_plot_index_2) * op_par[TH2_plot_index_2].delta;
        h_eff->Fill(par_1 + 0.00001, par_2+ 0.00001, eff);
        h_pur->Fill(par_1+ 0.00001, par_2+ 0.00001, pur);
        h_pureff->Fill(par_1+ 0.00001, par_2+ 0.00001, pur*eff);

        if(pur*eff > best_pureff) {
         //if(pur < 0.47) continue;
            best_pureff = pur*eff;
            best_pur = pur;
            best_eff = eff;

            for(int j = 0; j < num_op ;j++){
                best_cut[j] = op_par[j].inf_limit + index.at(i_c).at(j) *  op_par[j].delta;
            }
        }
    }


  string folder = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/shower_studies/pur_eff_optimization/2_shower_cut";


  gStyle->SetPadRightMargin(0.175);
    TCanvas *c1 = new TCanvas();
  h_eff->SetStats(0);
    h_eff->SetTitle(("Efficiency;" + op_par[TH2_plot_index_1].op_name + ";" +op_par[TH2_plot_index_2].op_name ).c_str());
    h_eff->Draw("colz");
    c1->SaveAs((folder + "eff_graph" +op_par[TH2_plot_index_1].op_name + "_" + op_par[TH2_plot_index_2].op_name + ".pdf").c_str());

  TCanvas *c2 = new TCanvas();
  h_pur->SetStats(0);
  h_pur->SetTitle(("Purity;" + op_par[TH2_plot_index_1].op_name + ";" +op_par[TH2_plot_index_2].op_name ).c_str());
  h_pur->Draw("colz");
  c2->SaveAs((folder + "pur_graph" +op_par[TH2_plot_index_1].op_name + "_" + op_par[TH2_plot_index_2].op_name + ".pdf").c_str());

  TCanvas *c3 = new TCanvas();
  h_pureff->SetStats(0);
  h_pureff->SetTitle(("Purity * Eff;" + op_par[TH2_plot_index_1].op_name + ";" +op_par[TH2_plot_index_2].op_name ).c_str());
  h_pureff->Draw("colz");
  c3->SaveAs((folder + "pur_eff_graph" +op_par[TH2_plot_index_1].op_name + "_" + op_par[TH2_plot_index_2].op_name + ".pdf").c_str());

  cout << "BEST CUT: ";
  for(int j = 0; j < num_op ;j++){
    cout <<  op_par[j].op_name << " "<< best_cut[j] << ", ";
  }
  cout << "pureff: " << best_pureff << endl;
  cout << "pur: " << best_pur << endl;
  cout << "eff: " << best_eff << endl;


}