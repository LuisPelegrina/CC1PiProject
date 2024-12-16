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

bool Normalize = false;
void chi2_confusion_only_mu_score()
{
  GenerateDictionaries();
    TTree *tree;
    TFile *input;

  Cut_Parameters cut_p;
  cut_p.use_default_cut_p();

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

   const int num_op = 3;
    struct OptimizationParameters op_par[num_op] = {
        {"max_muon_like_score", 5, 0, 100},
        {"min_proton_score", 5, 60, 100},
        {"min_TL", 2.5, 0, 15}
    };

    int n[num_op];
    for(int i = 0 ;i < num_op ; i++) {
        n[i] = (op_par[i].sup_limit - op_par[i].inf_limit)/op_par[i].delta + 1;
    }    
    
    int size = 1;
    for(int i_op = 0; i_op < num_op; i_op++) {
        size *= n[i_op];
    }
    cout << size << endl;


    double cont_muon_true[size];
    double cont_muon_miss[size];
    double cont_other_confusion[size];
    double cont_other_true[size];

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

        cont_muon_true[i] = 0;
        cont_muon_miss[i] = 0;
        cont_other_confusion[i] = 0;
        cont_other_true[i] = 0;
    }

    for(int i_e = 0; i_e < nEntries; ++i_e) {
        tree->GetEntry(i_e);

        if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

        double w = 1;
        if(Normalize) w = slice->weight;

        for(int i_p = 0; i_p <slice->pandora_particle.size(); i_p++) {
            Reco_Particle particle = slice->pandora_particle.at(i_p);

            //if(particle.purity < 0.8) continue;
            //if(particle.completeness < 0.8) continue;
    
            if(!particle.is_track(cut_p)) continue;
            if(!particle.is_primary(slice->primary_vertex_reco.vertex_cm, cut_p)) continue;
            if(!particle.pass_quality_cuts(cut_p, slice->hits)) continue;

            double muon_like_chi2_score = particle.chi2_score.muon_score;
            //if(particle.chi2_score.pion_score < muon_like_chi2_score) muon_like_chi2_score = particle.chi2_score.pion_score;
            //if(particle.track_lenght > 100) continue;

            bool is_candidate = false;

            if((particle.track_lenght > cut_p.chi2_min_TL) && (particle.chi2_score.proton_score > cut_p.chi2_max_proton_score) && ( muon_like_chi2_score < cut_p.chi2_min_muon_score)) is_candidate = true;

            if(is_candidate) {
                h_muonlike_proton[pdg_map[abs(particle.matched_pdg)]]->Fill(muon_like_chi2_score, particle.chi2_score.proton_score);
                h_muonlike_other[pdg_map[abs(particle.matched_pdg)]]->Fill(muon_like_chi2_score, particle.chi2_score.kaon_score);
            }



            if((abs(particle.matched_pdg) == 13) || (abs(particle.matched_pdg) == 211)) {
               // cout << particle.chi2_score.muon_score << endl;
                h_tl[0]->Fill(particle.track_lenght);
            } else {
                h_tl[1]->Fill(particle.track_lenght);
            }

            for(int i_c = 0; i_c < size; i_c++) {
                double max_muon_like_score = op_par[0].inf_limit + index.at(i_c).at(0)*op_par[0].delta;
                double min_proton_score = op_par[1].inf_limit + index.at(i_c).at(1)*op_par[1].delta; 
                double min_track_lenght = op_par[2].inf_limit + index.at(i_c).at(2)*op_par[2].delta;   
                //double max_proton_score = op_par[3].inf_limit + index.at(i_c).at(3)*op_par[3].delta;    
                //double min_only_tl =  op_par[3].inf_limit + index.at(i_c).at(3)*op_par[3].delta;            

                bool pass_cut = false;
                
                //if((muon_like_chi2_score < particle.chi2_score.proton_score) && (particle.track_lenght > min_track_lenght )) pass_cut = true;

                //if(particle.track_lenght > min_only_tl) pass_cut = true;
                if((muon_like_chi2_score < max_muon_like_score) && (particle.chi2_score.proton_score > min_proton_score) && (particle.track_lenght > min_track_lenght )) pass_cut = true;
                if((abs(particle.matched_pdg) == 13) || (abs(particle.matched_pdg) == 211)) {
                    if(pass_cut) cont_muon_true[i_c]++;
                    if(!pass_cut) cont_muon_miss[i_c]++;
                } else {
                    if(pass_cut) cont_other_confusion[i_c]++;   
                    if(!pass_cut) cont_other_true[i_c]++;      
                }
    
            }
        }

    }

    //SHOW OPTIMIZATION

    double best_pureff = 0;
    double best_cut[num_op];
    double best_cont_muon_true;
    double best_cont_muon_miss;
    double best_cont_other_confusion;
    double best_cont_other_true;
    //SHOW BACKGROUND AND SIGNAL 
    for(int i_c = 0; i_c < size; i_c++) {

        double eff = cont_muon_true[i_c]/(cont_muon_miss[i_c] + cont_muon_true[i_c]);
        double pur = cont_muon_true[i_c]/(cont_other_confusion[i_c] + cont_muon_true[i_c]);

        cout << cont_muon_miss[i_c] + cont_muon_true[i_c] << " CUT: " ;

        for(int j = 0; j < num_op ;j++){
           cout << op_par[j].inf_limit + index.at(i_c).at(j) * op_par[j].delta << " ";
        }
        cout << "pur: " << pur << " eff: " << eff  << " purr*eff: " << pur*eff <<  endl;

        if(pur*eff > best_pureff) {
            best_pureff = pur*eff;

            for(int j = 0; j < num_op ;j++){
                best_cut[j] = op_par[j].inf_limit + index.at(i_c).at(j) *  op_par[j].delta;
                best_cont_muon_true = cont_muon_true[i_c];
                best_cont_muon_miss = cont_muon_miss[i_c];
                best_cont_other_confusion = cont_other_confusion[i_c];
                best_cont_other_true = cont_other_true[i_c];
            }
        }
    }

    cout << "BEST CUT: ";
    for(int j = 0; j < num_op ;j++){
        cout <<  op_par[j].op_name << " "<< best_cut[j] << ", ";
    }
    cout << "pureff: " << best_pureff << endl;
    
    TH2 *h_confusion = new TH2D("h", "h",2,0,2,2,0,2);
    TCanvas *CASVA = new TCanvas();
    h_confusion->SetStats(0);
    
    //h_confusion->Fill(0., 0., best_cont_muon_true*1.0/(best_cont_muon_true+best_cont_muon_miss));
    //h_confusion->Fill(0., 1., best_cont_muon_miss*1.0/(best_cont_muon_true+best_cont_muon_miss));
    //h_confusion->Fill(1., 0., best_cont_other_confusion*1.0/(best_cont_other_confusion + best_cont_other_true));
    //h_confusion->Fill(1., 1., best_cont_other_true*1.0/(best_cont_other_confusion + best_cont_other_true));
    
    h_confusion->Fill(0., 0., best_cont_muon_true);
    h_confusion->Fill(0., 1., best_cont_muon_miss);
    h_confusion->Fill(1., 0., best_cont_other_confusion);
    h_confusion->Fill(1., 1., best_cont_other_true);

    h_confusion->GetXaxis()->SetBinLabel(1, "#mu like true");
    h_confusion->GetXaxis()->SetBinLabel(2, "others true");
    
    h_confusion->GetYaxis()->SetBinLabel(1, "#mu like reco");
    h_confusion->GetYaxis()->SetBinLabel(2, "others reco");

    h_confusion->Draw("colz");

    for(int i = 0; i< 2; i++) {
      for(int j = 0; j< 2; j++) {

        //if(int( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))) == 0) continue;

        string sPN1 = to_string(int( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))));
        TText *tPN1 = new TText(i+0.5, j+0.6,sPN1.c_str());

        int TotalP = 0;
        for (int ji = 0; ji< 2; ji++) {
          TotalP +=  h_confusion->GetBinContent(h_confusion->GetBin(i +1, ji+1));
        }
        TString sPN = TString::Format("%8.2f", h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))*100./TotalP);
        sPN += "%";
        TText *tPN = new TText(i+0.5,j+.4,sPN);
      
        cout << TotalP<<  endl;
        if ( h_confusion->GetBinContent( h_confusion->GetBin(i+1, j+1))*100./TotalP > 90) {
          tPN1->SetTextColor(kWhite);
          tPN->SetTextColor(kWhite);
        }

        tPN1->SetTextSize(0.05);
        tPN->SetTextSize(0.05);
        tPN1->Draw();
        tPN->Draw();
      }  
    }




    gStyle->SetPadRightMargin(0.175);
    //Print Tracks and shower plots
    for(int i = 0; i < 4; i++) {
        TCanvas *c_i = new TCanvas();
        string particle;
        if(i == 0) particle =  "Others";
        if(i == 1) particle =  "Muon";
        if(i == 2) particle =  "Pion";
        if(i == 3) particle =  "Proton";
        string title = particle + "; muon #chi^{2} score ; proton #chi^{2} score; events";
        string path = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/Chi2_graphs/";
        path = path + "Chi2_TH2_" + particle + ".pdf";

        h_muonlike_proton[i]->SetStats(0);
        h_muonlike_proton[i]->SetTitle(title.c_str());
        h_muonlike_proton[i]->Draw("colz");
        c_i->SaveAs(path.c_str());
    }
    
    
    for(int i = 0; i < 4; i++) {
        TCanvas *c_i = new TCanvas();
        string title;
        if(i == 0) title =  "Others"; 
        if(i == 1) title =  "Muon";
        if(i == 2) title =  "Pion";
        if(i == 3) title =  "Proton";
        title += "; muon_like_score ; undef_score; events";
        h_muonlike_other[i]->SetStats(0);

        h_muonlike_other[i]->SetTitle(title.c_str());
        //h_muonlike_other[i]->Draw("colz");
    }
    
     
    TCanvas *c2 = new TCanvas();
    h_tl[0]->SetLineColor(kBlue);
    h_tl[0]->SetFillColorAlpha(kBlue, 0.1);
    h_tl[0]->Draw("hist");

    h_tl[1]->SetLineColor(kRed);
    h_tl[1]->SetFillColorAlpha(kRed, 0.1);
    h_tl[1]->Draw("hist same");

    //input->Close();
}