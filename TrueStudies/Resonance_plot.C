
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

void Resonance_plot()
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
    TH1* h_E_Res = new TH1D("h_E_Res","h_E_Res", 80, 0, 6);
    TH1* h_E_Delta_Res = new TH1D("h_E_Delta_Res","h_E_Delta_Res", 80, 0, 6);
    TH1* h_E_PN_Res = new TH1D("h_E_PN_Res","h_E_PN_Res", 80, 0, 6);
    TH1* h_E_unclass_Res = new TH1D("h_E_unclass_Res","h_E_unclass_Res", 80, 0, 6);


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

            if(!true_interaction.is_selected_final_state("CC1Pi", 0.325)) continue;
            h_E_CC1Pi->Fill(true_interaction.E0_nu, w);

            if(true_interaction.interaction_mode != 1) continue;
            h_E_Res->Fill(true_interaction.E0_nu, w);

            vector<GeneratorParticle> Particles = true_interaction.gen_particles;
            int resonant_baryon_pdg = 0;
            for(int i_p = 0; i_p< Particles.size(); i_p++) {
                GeneratorParticle Particle = Particles.at(i_p);
                if(resonant_baryon_pdg != 0) continue;
                if(resonance_index_map[Particle.PDG] != 0) resonant_baryon_pdg = Particle.PDG;
                //cout << Particle.ID << " " << Particle.PDG << " " << resonance_index_map[Particle.PDG] << endl;
            }
            if(resonant_baryon_pdg == 0) cout << "Baryon not classified" << endl;
            //cout << "resonance PDG: " << resonant_baryon_pdg << endl;

            int resonance_index = resonance_index_map[resonant_baryon_pdg];

            string resonance_type = resonance_type_name_map[resonance_index];
            cont_resonance_freq[resonance_type]++;
            cont_resonance_freq_charge_dif[resonant_baryon_pdg]++;

            if(resonance_type_delta_or_N_map[resonance_index].compare("np") == 0) h_E_PN_Res->Fill(true_interaction.E0_nu, w);
            if(resonance_type_delta_or_N_map[resonance_index].compare("delta") == 0) h_E_Delta_Res->Fill(true_interaction.E0_nu, w);
            if(resonance_type_delta_or_N_map[resonance_index].compare("unclassified") == 0) h_E_unclass_Res->Fill(true_interaction.E0_nu, w);

        }
    }

    double n_res = h_E_Delta_Res->Integral() + h_E_PN_Res->Integral();
    cout << "Eventos totales: " << h_E_CC1Pi->Integral() << endl;
    cout << "Eventos resonantes: " << n_res  << " Porcentaje (total CC1Pi): " << n_res*100/h_E_CC1Pi->Integral() << " Porcentaje (Res): "<< n_res*100/n_res <<  endl;
  cout << "Eventos Delta resonance: " << h_E_Delta_Res->Integral()  << " Porcentaje (total CC1Pi): " << h_E_Delta_Res->Integral()*100/h_E_CC1Pi->Integral() << " Porcentaje (Res): "<< h_E_Delta_Res->Integral()*100/n_res<< endl;
  cout << "Eventos NP resonance: " <<  h_E_PN_Res->Integral() << " Porcentaje (total CC1Pi): " << h_E_PN_Res->Integral()*100/h_E_CC1Pi->Integral() << " Porcentaje (Res): "<< h_E_PN_Res->Integral()*100/n_res<< endl;


  cout << "Eventos según masa de barion resonante:" << endl;
    for (std::map<string,int>::iterator it=cont_resonance_freq.begin(); it!=cont_resonance_freq.end(); ++it) {
        cout << "Número de eventos con " << it->first << " : " << it->second << endl;
    }
    cout << endl;

    cout << "Eventos según PDG de baryon resonante:" << endl;
    for (std::map<int,int>::iterator it=cont_resonance_freq_charge_dif.begin(); it!=cont_resonance_freq_charge_dif.end(); ++it) {
        cout << "Número de eventos con " << it->first << " : " << it->second << endl;
    }


  TLegend *legRes = new TLegend(0.60,0.55,0.85,0.85);
  //legRes->AddEntry(hPnuResonB[0],"Resonances", "f");
  double no_res_perc =  100 - n_res*100./h_E_CC1Pi->Integral();
  double delta_res_perc =  h_E_Delta_Res->Integral()*100./h_E_CC1Pi->Integral();
  double np_res_perc =  h_E_PN_Res->Integral()*100./h_E_CC1Pi->Integral();

  std::ostringstream no_res_perc_rounded;
  no_res_perc_rounded  << std::fixed << std::setprecision(2) << no_res_perc;

  std::ostringstream delta_res_perc_rounded;
  delta_res_perc_rounded  << std::fixed << std::setprecision(2) << delta_res_perc;

  std::ostringstream np_res_perc_rounded;
  np_res_perc_rounded  << std::fixed << std::setprecision(2) << np_res_perc;

  legRes->AddEntry(h_E_Delta_Res,("#Delta Resonances: " + delta_res_perc_rounded.str() + " %").c_str(), "f");
  legRes->AddEntry(h_E_PN_Res,("N Resonances: " + np_res_perc_rounded.str() + " %").c_str(), "f");
  legRes->AddEntry(h_E_CC1Pi,("Not resonant CC1#pi events: " + no_res_perc_rounded.str()+ " %").c_str(), "f");


  TCanvas *c1 = new TCanvas();
    c1->cd();

    h_E_CC1Pi->SetStats(0);
    h_E_Res->SetStats(0);
    h_E_Delta_Res->SetStats(0);
    h_E_PN_Res->SetStats(0);
    h_E_unclass_Res->SetStats(0);

    //h_E_CC1Pi->Draw("hist");

    h_E_CC1Pi->Add(h_E_Delta_Res, -1);
    h_E_CC1Pi->Add(h_E_PN_Res, -1);
    h_E_CC1Pi->SetLineColor(kBlue +2 );
    h_E_CC1Pi->SetFillColorAlpha(kBlue +2 ,0.1);

    //h_E_Res->Draw("hist same");

    h_E_PN_Res->SetLineColor(kRed + 2);
    h_E_PN_Res->SetFillColorAlpha(kRed + 2,0.1);
    //h_E_PN_Res->Draw("hist same");

    h_E_Delta_Res->SetLineColor(kBlack);
    h_E_Delta_Res->SetFillColorAlpha(kBlack,0.1);
    //h_E_Delta_Res->Draw("hist same");



    auto hs = new THStack("hs","");
    hs->Add(h_E_Delta_Res);
    hs->Add(h_E_PN_Res);
    hs->Add(h_E_CC1Pi);
    hs->SetTitle(";E_{#nu} [GeV]; Events");
    hs->Draw("hist");
    legRes->Draw();

    TCanvas *c2 = new TCanvas();
    c2->cd();
    h_E_unclass_Res->Draw("hist");


    //for(int i_t = 0; i_t < num_trees; i_t++) input[i_t]->Close();
}