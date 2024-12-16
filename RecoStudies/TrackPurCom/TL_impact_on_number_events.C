
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

void TL_impact_on_number_events()
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

    TH1* h_TL_mu = new TH1D("h_TL_mu","h_TL_mu", 100, 0, 200);
    TH1* h_TL_pi = new TH1D("h_TL_pi","h_TL_pi", 100, 0, 200);
    TH2* h_TL_pi_mu = new TH2D("h_TL_pi_mu","h_TL_pi", 25, 0, 100, 25, 0, 100);
    int bin_max_cum = 100;
    int bin_min_cum = 0;
    int n_bins_cum = 1000;
  double pas_cum = (bin_max_cum-bin_min_cum)/n_bins_cum;
   TH2* h_TL_pi_mu_cum = new TH2D("h_TL_pi_mu_cum","h_TL_pi", n_bins_cum, bin_min_cum, bin_max_cum, n_bins_cum, bin_min_cum, bin_max_cum);

    double num_CC1pi_events = 0;
  int bin_max_cum_2D = 100;
  int bin_min_cum_2D = 0;
  int n_bins_cum_2D = 100;
  double pas_cum_2D = (bin_max_cum_2D-bin_min_cum_2D)/n_bins_cum_2D;
    TH1* h_TL_mu_pi_event_cum = new TH1D("h_TL_mu_pi_event_cum","h_TL_mu_pi_event_cum", n_bins_cum_2D, bin_min_cum_2D, bin_max_cum_2D);

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
          if(!true_interaction.is_selected_final_state("CC1Pi",0.325)) continue;
          if(!is_inside_AV(slice->primary_vertex_true.vertex_cm.X(), slice->primary_vertex_true.vertex_cm.Y(), slice->primary_vertex_true.vertex_cm.Z())) continue;
          num_CC1pi_events += w;
          h_TL_mu->Fill(true_interaction.get_g4_primary_p(13).at(0).TL, w);
          h_TL_pi->Fill(true_interaction.get_g4_primary_p(211).at(0).TL, w);
          h_TL_pi_mu->Fill(true_interaction.get_g4_primary_p(211).at(0).TL, true_interaction.get_g4_primary_p(13).at(0).TL, w);

          G4Particle prim_mu = true_interaction.get_g4_primary_p(13).at(0);
          G4Particle prim_pi = true_interaction.get_g4_primary_p(211).at(0);

          for(int i_b =1; i_b <= h_TL_mu_pi_event_cum->GetNbinsX(); i_b++) {
            double TL_min =  h_TL_mu_pi_event_cum->GetBinLowEdge(i_b);
            if((prim_pi.TL > TL_min) && (prim_mu.TL > TL_min)) {
              h_TL_mu_pi_event_cum->SetBinContent(i_b, h_TL_mu_pi_event_cum->GetBinContent(i_b) + w);
            }
          }


          for(int i_bx =0; i_bx < h_TL_pi_mu_cum->GetNbinsX(); i_bx++) {
            double TL_x = bin_min_cum_2D + i_bx*pas_cum_2D;
            for(int i_by =0; i_by < h_TL_pi_mu_cum->GetNbinsY(); i_by++) {
              double TL_y = bin_min_cum_2D + i_by*pas_cum_2D;
              if((prim_pi.TL > TL_x) && (prim_mu.TL > TL_y)) {
                h_TL_pi_mu_cum->Fill(TL_x+0.001, TL_y+0.001, w);
              }
            }
          }

        }
    }

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h_TL_mu->SetStats(0);


    h_TL_pi->SetStats(0);

    h_TL_mu->SetLineColor(kBlue+2);
  h_TL_mu->SetFillColorAlpha(kBlue +2,0.1);
  h_TL_mu->Draw("hist");
  h_TL_mu->SetTitle(";#mu track length [cm]; Events");
  string ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/TLCutImpactOnNEvents/muon_TL_distribution.pdf";
  c1->SaveAs(ruta.c_str());

  TCanvas *c2 = new TCanvas();
  c2->cd();
  h_TL_pi->SetLineColor(kBlue+2);
  h_TL_pi->SetFillColorAlpha(kBlue +2,0.1);
  h_TL_pi->Draw("hist");
  h_TL_pi->SetTitle(";#mu track length [cm]; Events");
  ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/TLCutImpactOnNEvents/pion_TL_distribution.pdf";
  c2->SaveAs(ruta.c_str());

  TCanvas *c3 = new TCanvas();
  c3->cd();
  h_TL_pi_mu->SetStats(0);
  h_TL_pi_mu->Draw("colz");
  h_TL_pi_mu->SetTitle(";#pi track length [cm]; #mu track length [cm]; Events");

  ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/TLCutImpactOnNEvents/muon_pion_TL_distribution.pdf";
  c3->SaveAs(ruta.c_str());

  TCanvas *c4 = new TCanvas();
  c4->cd();
  h_TL_mu_pi_event_cum->Scale(1./num_CC1pi_events);
  h_TL_mu_pi_event_cum->SetStats(0);
  h_TL_mu_pi_event_cum->SetTitle("Cumulative of Events;#mu-like Track length (cm); Event fraction");
  c4->SetGridy();
  c4->SetGridx();
  h_TL_mu_pi_event_cum->GetYaxis()->SetNdivisions(20);   // 8 divisions on the Y-axis
  h_TL_mu_pi_event_cum->Draw("hist");
  h_TL_mu_pi_event_cum->SetLineColor(kBlue+2);



  TH1F *h_TL_mu_pi_event_cum_percentage = (TH1F*)h_TL_mu_pi_event_cum->Clone("h_TL_mu_pi_event_cum_percentage");
  for(int i_b = 1; i_b <= h_TL_mu_pi_event_cum_percentage->GetNbinsX(); i_b++) {
    h_TL_mu_pi_event_cum_percentage->SetBinContent(i_b, 1 - h_TL_mu_pi_event_cum_percentage->GetBinContent(i_b));
  }
  h_TL_mu_pi_event_cum_percentage->SetLineColor(kRed+2);
  h_TL_mu_pi_event_cum_percentage->Draw("hist same");


  TLegend *leg  = new TLegend(0.50,0.35,0.88,0.58);
  leg->AddEntry(h_TL_mu_pi_event_cum,"% events remaining ","l");
  leg->AddEntry(h_TL_mu_pi_event_cum_percentage, "% events lost ", "l") ;
  leg->Draw();
  leg->SetFillColorAlpha(kRed, 0);
  c4->Update();

  ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/TLCutImpactOnNEvents/muon_pion_event_cum.pdf";
  c4->SaveAs(ruta.c_str());

  TCanvas *c5 = new TCanvas();
  c5->cd();
  h_TL_pi_mu_cum->Draw("colz");
  h_TL_pi_mu_cum->Scale(1./num_CC1pi_events);
  h_TL_pi_mu_cum->SetTitle("Cumulative of Events;#pi Track length (cm);#mu Track length (cm); Events");
  h_TL_pi_mu_cum->SetStats(0);
  h_TL_pi_mu_cum->SetContour(40);  // 20 contour levels

  TH2F *h_TL_pi_mu_cum_contour = (TH2F*)h_TL_pi_mu_cum->Clone("h_TL_mu_pi_event_cum_percentage");
  h_TL_pi_mu_cum_contour->SetContour(10);  // 20 contour levels
  h_TL_pi_mu_cum_contour->SetLineColor(kBlack);
  h_TL_pi_mu_cum_contour->SetLineStyle(8);
  h_TL_pi_mu_cum_contour->Draw("CONT3 same");
  c5->Update();


  ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MissRecoGraphs/TLCutImpactOnNEvents/muon_pion_kinematic_restriction.pdf";
  c5->SaveAs(ruta.c_str());
    //for(int i_t = 0; i_t < num_trees; i_t++) input[i_t]->Close();
}