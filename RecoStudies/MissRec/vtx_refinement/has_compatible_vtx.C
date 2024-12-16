#include "../../../Includes.h"
#include "../Graphs_utils.cpp"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}

void has_compatible_vtx()
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
  GenerateDictionaries();
  TTree *tree;
  TFile *input;
    
  Cut_Parameters cut_p;
  cut_p.min_track_lenght = 0;
  cut_p.min_track_score = 0;


  //Declare the variables
  Slice *slice = 0;

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_83k.root";
  input = new TFile(strRuta.c_str());
  tree =(TTree*)input->Get("tree");
  tree->SetBranchAddress("slice", &slice);
  int nEntries = tree->GetEntries();

  const int num_cuts = 3;
  string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut"};

  TH1 *h_has_comp = new TH1D("h_has_comp", "h_has_comp",20,0,180);
  TH1 *h_miss_comp = new TH1D("h_miss_comp", "h_miss_comp",20,0,180);
  TH1 *h_random_slice_comp = new TH1D("h_random_slice_comp", "h_miss_comp",20,0,180);


  TH1 *h_has_comp_TS = new TH1D("h_has_comp_TS", "h_has_comp",20,0,180);
  TH1 *h_miss_comp_TS = new TH1D("h_miss_comp_TS", "h_miss_comp",20,0,180);
  TH1 *h_random_slice_comp_TS = new TH1D("h_random_slice_comp_TS", "h_miss_comp",20,0,180);

  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

    if(!slice->true_interaction.is_selected_final_state("N_CC1Pi")) continue;
    if(!slice->true_interaction.is_vertex_contained(10)) continue;

    bool pass_cut = true;
    //Plot only if correctly reconstructed
    for(int i_cut = 0; i_cut < num_cuts; i_cut++) {
      if(!pass_cut) continue;
      if(!slice->pass_analysis_cut(cut_p, cut_name[i_cut])) {
        pass_cut = false;
      }
    }
    if(!pass_cut) continue;

    double distance_to_true_vertex = (slice->primary_vertex_true.vertex_cm - slice->primary_vertex_reco.vertex_cm).Mag();

    G4Particle primary_muon = slice->true_interaction.get_g4_primary_p(13).at(0);
    G4Particle primary_pion = slice->true_interaction.get_g4_primary_p(211).at(0);

    if ((primary_muon.TL < 3) || (primary_pion.TL < 3)) continue;
    //Really bad vertex
    if(distance_to_true_vertex < 2) continue;

    //Check if there is an extra slice or not
    bool extra_slice = false;

    double run_Id = slice->run_ID;
    double subrun_Id = slice->subrun_ID;
    double event_Id = slice->event_ID;

    tree->GetEntry(i_e-1);
    if((run_Id = slice->run_ID) && (subrun_Id = slice->subrun_ID) && (event_Id = slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e+1);
    if((run_Id = slice->run_ID) && (subrun_Id = slice->subrun_ID) && (event_Id = slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e-2);
    if((run_Id = slice->run_ID) && (subrun_Id = slice->subrun_ID) && (event_Id = slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e+2);
    if((run_Id = slice->run_ID) && (subrun_Id = slice->subrun_ID) && (event_Id = slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e);

    int mother_PDG = -1;
    bool first = false;
    if(extra_slice) {
      for (int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
        if((slice->pandora_particle.at(i_p).track_start -slice->primary_vertex_reco.vertex_cm).Mag() < 2) {
          if(first) continue;
          first = true;
          for (int i_pg4 = 0; i_pg4 < slice->true_interaction.g4_particles.size(); i_pg4++) {
            G4Particle g4_p = slice->true_interaction.g4_particles.at(i_pg4);
            if (slice->pandora_particle.at(i_p).true_track_id == g4_p.ID) {
              for (int i_pg4_2 = 0; i_pg4_2 < slice->true_interaction.g4_particles.size(); i_pg4_2++) {
                if (slice->true_interaction.g4_particles.at(i_pg4_2).ID == g4_p.mother) {
                  mother_PDG = slice->true_interaction.g4_particles.at(i_pg4_2).PDG;
                }
              }
            }
          }
        }
      }
    }

    bool has_compatible_vtx = false;
    for(int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      if((slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag() < 4)  has_compatible_vtx = true;
    }

    if(has_compatible_vtx)  {
      h_has_comp->Fill(primary_muon.P0.Angle(primary_pion.P0)*360/(2*TMath::Pi()));
    } else if((!has_compatible_vtx) && extra_slice && ((mother_PDG == 2112) || (mother_PDG == 111))) {
      h_random_slice_comp->Fill(primary_muon.P0.Angle(primary_pion.P0)*360/(2*TMath::Pi()));
    } else {
      h_miss_comp->Fill(primary_muon.P0.Angle(primary_pion.P0)*360/(2*TMath::Pi()));
    }


    for(int i_p = 0; i_p < slice->pandora_particle.size() ;i_p++) {
      if(((slice->pandora_particle.at(i_p).track_start - slice->primary_vertex_true.vertex_cm).Mag() < 4) ||((slice->pandora_particle.at(i_p).track_end - slice->primary_vertex_true.vertex_cm).Mag() < 4))  has_compatible_vtx = true;
    }

    if(has_compatible_vtx) {
      h_has_comp_TS->Fill(primary_muon.P0.Angle(primary_pion.P0)*360/(2*TMath::Pi())) ;
    } else if((!has_compatible_vtx) && extra_slice && ((mother_PDG == 2112) || (mother_PDG == 111))) {
      h_random_slice_comp_TS->Fill(primary_muon.P0.Angle(primary_pion.P0)*360/(2*TMath::Pi()));
    } else {
      h_miss_comp_TS->Fill(primary_muon.P0.Angle(primary_pion.P0)*360/(2*TMath::Pi()));
    }

  }

  h_miss_comp->SetLineColor(kRed+2);
  h_random_slice_comp->SetLineColor(kGreen+2);
  THStack *hst1 =new  THStack();
  hst1->Add(h_has_comp);
  hst1->Add(h_miss_comp);
  hst1->Add(h_random_slice_comp);
  hst1->SetTitle("Checking all vtx;#mu#pi true angle [#circ];Events");
  hst1->Draw();
  hst1->SetMaximum(hst1->GetMaximum()*1.6);
  int lost = h_miss_comp->Integral();
  int saved =  h_has_comp->Integral();
  int extra =  h_random_slice_comp->Integral();

  TLegend *leg = new TLegend(0.55, 0.7, .88, .88);
  leg->AddEntry(h_has_comp, ("Saved: " + to_string(saved)).c_str(),"lf");
  leg->AddEntry(h_miss_comp, ("Still lost: "+ to_string(lost)).c_str(),"lf");
  leg->AddEntry(h_random_slice_comp,("Extra slice: "+ to_string(extra )).c_str(),"lf");
  leg->Draw();
  TCanvas *c2 = new TCanvas();

  h_miss_comp_TS->SetLineColor(kRed+2);
  h_random_slice_comp_TS->SetLineColor(kGreen+2);
  THStack *hst2 =new  THStack();
  hst2->Add(h_has_comp_TS);
  hst2->Add(h_miss_comp_TS);
  hst2->Add(h_random_slice_comp_TS);
  hst2->SetTitle("Checking vtx, track start & track end;#mu#pi true angle [#circ];Events");
  hst2->SetMaximum(hst2->GetMaximum()*1.6);
  hst2->Draw();

  lost = h_miss_comp_TS->Integral();
  saved =  h_has_comp_TS->Integral();
  extra =  h_random_slice_comp_TS->Integral();
  TLegend *leg2 = new TLegend(0.55, 0.7, .88, .88);
  leg2->AddEntry(h_has_comp_TS, ("Saved: " + to_string(saved)).c_str(),"lf");
  leg2->AddEntry(h_miss_comp_TS, ("Still lost: "+ to_string(lost)).c_str(),"lf");
  leg2->AddEntry(h_random_slice_comp_TS, ("Extra slice: "+ to_string(extra )).c_str(),"lf");
  leg2->Draw();
}