
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

void TL_impact_on_energy_spectrum()
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

    int num_TL_steps = 2;
    double TL_steps = 3;
    TH1* h_E[num_TL_steps];
    for(int i_s = 0; i_s < num_TL_steps; i_s++) {
      h_E[i_s]= new TH1D(("h_E" + to_string(i_s)).c_str(),"h_TL_mu_pi_event_cum", 50, 0, 5);
    }

    for(int i_t = 0; i_t < num_trees; i_t++) {
        int num_entries = tree[i_t]->GetEntries();
        //num_entries = 100;
        for (int i_e = 0; i_e < num_entries; ++i_e) {
          tree[i_t]->GetEntry(i_e);
          if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;

          if(slice->slice_ID != 0) continue;
          TrueInteraction true_interaction = slice->true_interaction;
          //if((tm_purity < 0.5)||(tm_completeness < 0.5)) continue;

          double w = 1;
          if (Normalize) w = slice->weight;
          if(!true_interaction.is_selected_final_state("CC1Pi", 0.375)) continue;
          if(!is_inside_AV(slice->primary_vertex_true.vertex_cm.X(), slice->primary_vertex_true.vertex_cm.Y(), slice->primary_vertex_true.vertex_cm.Z())) continue;

          G4Particle prim_mu = true_interaction.get_g4_primary_p(13).at(0);
          G4Particle prim_pi = true_interaction.get_g4_primary_p(211).at(0);

          for(int i_s = 0; i_s < num_TL_steps; i_s++) {
            double TL_min = i_s*TL_steps;
            if((prim_pi.TL > TL_min) && (prim_mu.TL > TL_min)) h_E[i_s]->Fill(slice->true_interaction.E0_nu);
          }

        }
    }

  TCanvas *c1 = new TCanvas("c1", "Split Canvas", 800, 600);

  // Define two pads with unequal space
  TPad *pad1 = new TPad("pad1", "Pad 1", 0.0, 0.3, 1.0, 1.0); // Larger pad (top)
  TPad *pad2 = new TPad("pad2", "Pad 2", 0.0, 0.0, 1.0, 0.3); // Smaller pad (bottom)

  // Set margins for each pad if needed
  pad2->SetTopMargin(0.02);    // Reduces gap between pads
  pad2->SetBottomMargin(0.3);    // Reduces gap between pads
  pad1->SetBottomMargin(0.05);    // Reduces gap between pads

  // Draw pads on the canvas
  pad1->Draw();
  pad2->Draw();


  pad1->cd();
    TLegend *leg  = new TLegend(0.50,0.65,0.88,0.88);
    for(int i_s = 0; i_s < num_TL_steps; i_s++) {
      h_E[i_s]->SetTitle(";#nu Energy [GeV]; Events");
      h_E[i_s]->SetStats(0);
      //h_E[i_s]->SetLineWidth(3);
      if(i_s == 0) {
        h_E[i_s]->SetLineColor(kBlack);
      } else if(i_s == 1) {
        h_E[i_s]->SetLineColor(kBlue +2);
      } else if(i_s == 2) {
        h_E[i_s]->SetLineColor(kRed +2);
      } else if(i_s == 3) {
        h_E[i_s]->SetLineColor(kGreen +2);
      } else if(i_s == 4) {
        h_E[i_s]->SetLineColor(kOrange +2);
      } else if(i_s == 5) {
        h_E[i_s]->SetLineColor(kMagenta +2);
      }
      if(i_s == 0) h_E[i_s]->Draw("hist");
      h_E[i_s]->Draw("hist same");
      double TL_min = i_s*TL_steps;
      if(i_s == 0) {
        leg->AddEntry(h_E[i_s],"No track length requirement", "f");
      } else {
        int integerPart = static_cast<int>(std::round(TL_min)); // Use std::floor or std::ceil if you want specific rounding

        leg->AddEntry(h_E[i_s],("Track length > " + to_string(integerPart) + " cm").c_str(), "f");
      }
    }
    leg->Draw();
  string ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/TLCutImpactOnEnergy/energy_spectrum_comparison.pdf";
  c1->SaveAs(ruta.c_str());


  TLegend *leg2  = new TLegend(0.50,0.25,0.88,0.48);
  pad2->cd();
  TH1* h_E_clone[num_TL_steps];
  for(int i_s = 0; i_s < num_TL_steps; i_s++) {
    h_E_clone[i_s] = (TH1F*)h_E[i_s]->Clone("h_TL_mu_pi_event_cum_percentage");
  }
  for(int i_s = 0; i_s < num_TL_steps; i_s++) {
    h_E_clone[i_s]->SetTitle(";#nu Energy [GeV]; % Total events");
    h_E_clone[i_s]->Divide(h_E[0]);
    h_E_clone[i_s]->SetStats(0);

    h_E_clone[i_s]->GetXaxis()->SetTitleSize(0.13); // Adjust these values as needed
    h_E_clone[i_s]->GetXaxis()->SetLabelSize(0.1); // Adjust these values as needed
    h_E_clone[i_s]->GetYaxis()->SetTitleSize(0.13);
    h_E_clone[i_s]->GetYaxis()->SetTitleOffset(0.45);
    h_E_clone[i_s]->GetYaxis()->SetLabelSize(0.1); // Adjust these values as needed
    if(i_s == 0) {
      h_E_clone[i_s]->SetLineColor(kBlack);
    } else if(i_s == 1) {
      h_E_clone[i_s]->SetLineColor(kBlue +2);
    } else if(i_s == 2) {
      h_E_clone[i_s]->SetLineColor(kRed +2);
    } else if(i_s == 3) {
      h_E_clone[i_s]->SetLineColor(kGreen +2);
    } else if(i_s == 4) {
      h_E_clone[i_s]->SetLineColor(kOrange +2);
    } else if(i_s == 5) {
      h_E_clone[i_s]->SetLineColor(kMagenta +2);
    }
    if(i_s == 0) h_E_clone[i_s]->Draw("hist");
    h_E_clone[i_s]->Draw("hist same");
    double TL_min = i_s*TL_steps;
    leg2->AddEntry(h_E_clone[i_s],("Track length > " + to_string(int(TL_min)) + "cm").c_str(), "f");
  }
  //leg2->Draw();
  ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/TLCutImpactOnEnergy/uncut_comparison.pdf";
  c1->SaveAs(ruta.c_str());
}