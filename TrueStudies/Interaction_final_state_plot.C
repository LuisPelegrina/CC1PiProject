
#include "../Includes.h"
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

void Interaction_final_state_plot()
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
        if (i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
        input[i_t] = new TFile(strRuta.c_str());
        tree[i_t] = (TTree *) input[i_t]->Get("tree");
        tree[i_t]->SetBranchAddress("slice", &slice);
    }

    TH1* h_E_CC1Pi = new TH1D("h_E_CC1Pi","h_E_CC1Pi", 80, 0, 6);
    TH1* h_E_NCC1Pi_Cohless = new TH1D("h_E_NCC1Pi_Cohless","h_E_NCC1Pi_Cohless", 80, 0, 6);
    TH1* h_E_Coh = new TH1D("h_E_Coh","h_E_Coh", 80, 0, 6);
    TH1* h_E_1p_CC1pi = new TH1D("h_E_1p_CC1pi","h_E_1p_CC1pi", 80, 0, 6);
    TH1* h_E_2p_CC1pi = new TH1D("h_E_2p_CC1pi","h_E_2p_CC1pi", 80, 0, 6);
    TH1* h_CC1Pi_rest = new TH1D("h_CC1Pi_rest","h_CC1Pi_rest", 80, 0, 6);

    map<string, int> cont_resonance_freq;
    map<int, int> cont_resonance_freq_charge_dif;

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

          if(!true_interaction.is_selected_final_state("CC1Pi")) continue;

          h_E_CC1Pi->Fill(true_interaction.E0_nu, w);
          if(true_interaction.is_selected_final_state("Coherent")) {
            h_E_Coh->Fill(true_interaction.E0_nu, w);
          } else if (true_interaction.is_selected_final_state("N_CC1Pi")) {
            h_E_NCC1Pi_Cohless->Fill(true_interaction.E0_nu, w);
          } else if (true_interaction.is_selected_final_state("1P_CC1Pi")) {
            h_E_1p_CC1pi->Fill(true_interaction.E0_nu, w);
          } else if (true_interaction.is_selected_final_state("plus2P_CC1Pi")) {
            h_E_2p_CC1pi->Fill(true_interaction.E0_nu, w);
          } else {
            h_CC1Pi_rest->Fill(true_interaction.E0_nu, w);
          }
        }
    }

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h_E_CC1Pi->SetStats(0);
    h_E_NCC1Pi_Cohless->SetStats(0);
    h_E_Coh->SetStats(0);
    h_E_1p_CC1pi->SetStats(0);
    h_CC1Pi_rest->SetStats(0);
    h_E_2p_CC1pi->SetStats(0);

    h_E_CC1Pi->SetLineColor(kBlue+2);
    h_E_CC1Pi->SetFillColorAlpha(kBlue+2,0.1);

  h_E_NCC1Pi_Cohless->SetLineColor(kRed +2);
  h_E_NCC1Pi_Cohless->SetFillColorAlpha(kRed +2,0.1);

  h_E_Coh->SetLineColor(kBlack);
  h_E_Coh->SetFillColorAlpha(kBlack,0.1);

  h_E_1p_CC1pi->SetLineColor(kBlue +2);
  h_E_1p_CC1pi->SetFillColorAlpha(kBlue,0.1);

  h_CC1Pi_rest->SetLineColor(kOrange+2);
  h_CC1Pi_rest->SetFillColorAlpha(kOrange+2,0.1);

  h_E_2p_CC1pi->SetLineColor(kGreen+2);
  h_E_2p_CC1pi->SetFillColorAlpha(kGreen,0.1);

    TLegend *legRes = new TLegend(0.50,0.65,0.88,0.88);
    legRes->AddEntry(h_E_Coh,"Coherent #pi Production", "f");
    legRes->AddEntry(h_E_NCC1Pi_Cohless,"NCC1#pi wo Coherent", "f");
    legRes->AddEntry(h_E_1p_CC1pi,"CC1#pi + 1p", "f");
    legRes->AddEntry(h_E_2p_CC1pi,"CC1#pi + 2p", "f");
    legRes->AddEntry(h_CC1Pi_rest,"Other CC1#pi", "f");

    auto hs = new THStack("hs","");
    hs->Add(h_E_Coh);
    hs->Add(h_E_NCC1Pi_Cohless);
    hs->Add(h_E_1p_CC1pi);
    hs->Add(h_E_2p_CC1pi);
    hs->Add(h_CC1Pi_rest);
    hs->SetTitle(";E_{#nu} [GeV]; Events");
    hs->Draw("hist");
    legRes->Draw();

    TCanvas *c2 = new TCanvas();
    c2->cd();
    h_E_CC1Pi->Draw("hist");

    cout << "Total CC1Pi events: " << h_E_CC1Pi->Integral() << endl;
    cout << "Total CC1Pi events (Sum of subsections): " << h_CC1Pi_rest->Integral()  + h_E_Coh->Integral() + h_E_NCC1Pi_Cohless->Integral() + h_E_1p_CC1pi->Integral()  +  h_E_2p_CC1pi->Integral() << endl;
    cout << "Coherent pions: " << h_E_Coh->Integral() << endl;
    cout << "NCC1Pi wo Coherent: " << h_E_NCC1Pi_Cohless->Integral() << endl;
    cout << "CC1#pi + 1p: " << h_E_1p_CC1pi->Integral() << endl;
    cout << "CC1#pi + 2p: " << h_E_2p_CC1pi->Integral() << endl;
    cout << "CC1#pi + Rest: " << h_CC1Pi_rest->Integral() << endl;



    //for(int i_t = 0; i_t < num_trees; i_t++) input[i_t]->Close();
}