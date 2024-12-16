
#include "../Includes.h"
#include "Resonance_utils.cpp"
bool save_graphs = true;
bool Normalize = true;


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

void CC1PiOldComp()
{
  GenerateDictionaries();

    Cut_Parameters cut_p;

    //Declare the variables
    string strRuta;

    const int num_trees = 1;

    TTree *tree[num_trees];
    TFile *input[num_trees];


    Slice* slice = nullptr;

    for(int i_t = 0; i_t < num_trees; i_t++) {
        if (i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_82p7k.root";
        input[i_t] = new TFile(strRuta.c_str());
        tree[i_t] = (TTree *) input[i_t]->Get("tree");
        tree[i_t]->SetBranchAddress("slice", &slice);
    }

    TH1* h_E_CC1Pi = new TH1D("h_E_CC1Pi","h_E_CC1Pi", 80, 0, 6);
    TH1* h_E_CC1PiPi0 = new TH1D("h_E_CC1PiPi0","h_E_Res", 80, 0, 6);

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

            if(!true_interaction.is_selected_final_state("old_CC1Pi", 0.325)) continue;

            if(true_interaction.is_selected_final_state("CC1Pi", 0.325)) {
              h_E_CC1Pi->Fill(slice->true_interaction.E0_nu);
            } else {
              h_E_CC1PiPi0->Fill(slice->true_interaction.E0_nu);
            }

        }
    }



  TLegend *leg = new TLegend(0.4, 0.5, .88, .85);



  double CC1Pi_cont =  h_E_CC1Pi->Integral()*100./(h_E_CC1Pi->Integral() + h_E_CC1PiPi0->Integral());
  double CC1_piothers_cont =  h_E_CC1PiPi0->Integral()*100./(h_E_CC1Pi->Integral() + h_E_CC1PiPi0->Integral());

  std::ostringstream CC1pi_rounded;
  CC1pi_rounded  << std::fixed << std::setprecision(2) << CC1Pi_cont;

  std::ostringstream CC1pi_others_rounded;
  CC1pi_others_rounded  << std::fixed << std::setprecision(2) << CC1_piothers_cont;

  leg->AddEntry(h_E_CC1Pi, ("CC1#pi: " + CC1pi_rounded.str() + " %").c_str(),"lf");
  leg->AddEntry(h_E_CC1PiPi0, ("Excluded 1#mu1#pi events: " + CC1pi_others_rounded.str() + " %").c_str(),"lf");
  TCanvas *c1 = new TCanvas();
    c1->cd();

    h_E_CC1Pi->SetStats(0);
    h_E_CC1PiPi0->SetStats(0);
    //h_E_CC1Pi->Draw("hist");
    h_E_CC1Pi->SetLineColor(kBlue +2 );
    h_E_CC1Pi->SetFillColorAlpha(kBlue +2 ,0.1);

    //h_E_Res->Draw("hist same");

    h_E_CC1PiPi0->SetLineColor(kRed + 2);
    h_E_CC1PiPi0->SetFillColorAlpha(kRed + 2,0.1);
    //h_E_PN_Res->Draw("hist same");


    auto hs = new THStack("hs","");
    hs->Add(h_E_CC1Pi);
    hs->Add(h_E_CC1PiPi0);
    hs->SetTitle(";E_{#nu} [GeV]; Events");
    hs->Draw("hist");
  leg->Draw();
}