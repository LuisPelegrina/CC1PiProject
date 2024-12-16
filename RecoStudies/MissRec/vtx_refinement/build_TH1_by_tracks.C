#include "../../../Includes.h"
#include "../Graphs_utils.cpp"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}








void build_TH1_by_tracks()
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


  std::vector<TGraph*> gr_general_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_muon_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_pion_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_proton_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_other_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_muon_primary_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_pion_primary_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_proton_primary_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_other_primary_track_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_primary_vertex = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_primary_vertex_truth = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_secondary_vertex = {new TGraph(), new TGraph(), new TGraph()};


  TGraph2D* gr_general_sp_2D = new TGraph2D();
  TGraph2D* gr_muon_track_sp_2D = new TGraph2D();
  TGraph2D* gr_pion_track_sp_2D = new TGraph2D();
  TGraph2D* gr_proton_track_sp_2D = new TGraph2D();
  TGraph2D* gr_other_track_sp_2D = new TGraph2D();
  TGraph2D* gr_muon_primary_track_sp_2D = new TGraph2D();
  TGraph2D* gr_pion_primary_track_sp_2D = new TGraph2D();
  TGraph2D* gr_proton_primary_track_sp_2D = new TGraph2D();
  TGraph2D* gr_other_primary_track_sp_2D = new TGraph2D();
  TGraph2D* gr_primary_vertex_2D = new TGraph2D();
  TGraph2D* gr_primary_vertex_truth_2D = new TGraph2D();
  TGraph2D* gr_secondary_vertex_2D = new TGraph2D();

  std::vector<TCanvas*> c = {new TCanvas("c1", "My Canvas", 500, 400), new TCanvas("c2", "My Canvas", 500, 400), new TCanvas("c3", "My Canvas", 500, 400)};
  c[0]->SetWindowPosition(10, 50);
  c[1]->SetWindowPosition(700, 50);
  c[2]->SetWindowPosition(10, 500);

  TCanvas* c2D = new TCanvas("c4", "My Canvas", 500, 400);
  c2D->SetWindowPosition(700, 500);

  TCanvas* button = new TCanvas("c12", "My Canvas", 200, 200);
  button->SetWindowPosition(1400, 50);
  button->SetFillColor(kRed);
  button->Draw();

  const int num_cuts = 3;
  string cut_name[num_cuts] = {"no_cut", "is_clear_cosmic_cut", "reco_cut"};
   
  //pp = Primary Particle
  map<string,double> cases_cont;
  for(int i_e = 9000; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

    //if(slice->slice_ID != 0) continue;
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

    G4Particle primary_muon = slice->true_interaction.get_g4_primary_p(13).at(0);
    G4Particle primary_pion = slice->true_interaction.get_g4_primary_p(211).at(0);

    int cont_primary_p = 0;
    int cont_p = slice->pandora_particle.size();
    map<int, int> cont_by_pdg_pur_80;
    map<int, int> cont_by_pdg_pur_50;
    map<int, int> cont_by_pdg;

    int close_muon = 0;
    int close_pion = 0;
    for(int i_p = 0; i_p < slice->pandora_particle.size() ;i_p++) {
      Reco_Particle pandora_p = slice->pandora_particle.at(i_p);
        if(slice->pandora_particle.at(i_p).is_pandora_primary){
          cont_by_pdg[11111]++;
          cont_by_pdg[abs(pandora_p.matched_pdg)]++;
          if((pandora_p.track_start - slice->primary_vertex_reco.vertex_cm).Mag() < 8) {
            cont_by_pdg[22222]++;
            if(pandora_p.matched_pdg == 211) {
              close_pion++;
            } else if(pandora_p.matched_pdg == 13) {
              close_muon++;
            }
          }
          if((pandora_p.completeness > 0.8) && (pandora_p.purity > 0.8)) cont_by_pdg_pur_80[abs(pandora_p.matched_pdg)]++;
          if((pandora_p.completeness > 0.5) && (pandora_p.purity > 0.5)) cont_by_pdg_pur_50[abs(pandora_p.matched_pdg)]++;
          if((abs(pandora_p.matched_pdg) == 211) || (abs(pandora_p.matched_pdg) == 13)) continue;
          cont_by_pdg[-1]++;
        }
    }

    double distance_to_true_vertex = (slice->primary_vertex_true.vertex_cm - slice->primary_vertex_reco.vertex_cm).Mag();
    bool show = true;


    if ((primary_muon.TL < 3) || (primary_pion.TL < 3)) continue;
    //Really bad vertex
    if(distance_to_true_vertex < 2) continue;

    bool skip = false;
    bool has_vtx_comp = false;
    for(int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      if((slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag() < 4) has_vtx_comp = true;
    }

    double run_Id = slice->run_ID;
    double subrun_Id = slice->subrun_ID;
    double event_Id = slice->event_ID;

    //Check if there is an extra slice or not
    bool extra_slice = false;

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
          cout << "ID: " << slice->pandora_particle.at(i_p).true_track_id << endl;
          for (int i_pg4 = 0; i_pg4 < slice->true_interaction.g4_particles.size(); i_pg4++) {
            G4Particle g4_p = slice->true_interaction.g4_particles.at(i_pg4);
            if (slice->pandora_particle.at(i_p).true_track_id == g4_p.ID) {
              cout << "Mother: " << g4_p.mother << endl;
              for (int i_pg4_2 = 0; i_pg4_2 < slice->true_interaction.g4_particles.size(); i_pg4_2++) {
                if (slice->true_interaction.g4_particles.at(i_pg4_2).ID == g4_p.mother) {
                  cout << slice->true_interaction.g4_particles.at(i_pg4_2).ID << " Mother PDG: "
                       << slice->true_interaction.g4_particles.at(i_pg4_2).PDG << endl;
                  mother_PDG = slice->true_interaction.g4_particles.at(i_pg4_2).PDG;
                }
              }
            }
          }
        }
      }
    }
    if(extra_slice) cout << "EXTRA SLICE" << endl;
    if(extra_slice && ((mother_PDG == 2112) || (mother_PDG == 111))) skip = true;

    if(skip || !has_vtx_comp) continue;

    double mu_pi_angle =primary_muon.P0.Angle(primary_pion.P0)*360/(2*TMath::Pi());

    int num_canvas = slice->pandora_particle.size();
    cout << "VERTEXES" << endl;

    if(!show) continue;
    slice->print(cut_p, true, true, true, true,false, true);
    cout << "DISTANCE TO TRUE VERTEX" << distance_to_true_vertex << endl;
    cout << "PRINT ENDED" << endl;

    int TPC = 0;
    if(slice->primary_vertex_true.vertex_cm.X() > 0) TPC = 1;
    cout << "TPC Chosen" << endl;

    //bucle para los hits
    for(int i_h = 0; i_h < slice->hits.size(); i_h++) {
      if (slice->hits.at(i_h).TPC_ID != TPC) continue;
      Hit hit = slice->hits.at(i_h);
     // if(slice->get_distance_to_vertex(hit) > 200 ) continue;
      gr_general_hits[hit.Plane_ID]->SetPoint(gr_general_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
      if(hit.associated_pfp_ID == -1) continue;

      //Search for the pfp that the hit is associated to:
      int particle_index = 0;
      for(int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
        if(slice->pandora_particle.at(i_p).ID == hit.associated_pfp_ID) particle_index = i_p;
      }
      if(slice->pandora_particle.at( particle_index).is_pandora_primary) {
        if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 211) {
          gr_pion_primary_track_hits[hit.Plane_ID]->SetPoint(gr_pion_primary_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 13) {
          gr_muon_primary_track_hits[hit.Plane_ID]->SetPoint(gr_muon_primary_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 2212) {
          gr_proton_primary_track_hits[hit.Plane_ID]->SetPoint(gr_proton_primary_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else {
          gr_other_primary_track_hits[hit.Plane_ID]->SetPoint(gr_other_primary_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        }
      } else {
        if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 211) {
          gr_pion_track_hits[hit.Plane_ID]->SetPoint(gr_pion_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 13) {
          gr_muon_track_hits[hit.Plane_ID]->SetPoint(gr_muon_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 2212) {
          gr_proton_track_hits[hit.Plane_ID]->SetPoint(gr_proton_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        } else {
          gr_other_track_hits[hit.Plane_ID]->SetPoint(gr_other_track_hits[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
        }
      }

    }
    cout << "hits done" << endl;

    //bucle en las 3 vistas para el vértice primario
    for(int i_p = 0; i_p < 3; i_p++) {
      double vertex_drift_t = slice->primary_vertex_reco.vertex_wire.at(i_p).drift_t;
      int vertex_wire_ID = slice->primary_vertex_reco.vertex_wire.at(i_p).Wire_ID;

      gr_primary_vertex[i_p]->SetPoint(gr_primary_vertex[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
      gr_primary_vertex[i_p]->SetPoint(gr_primary_vertex[i_p]->GetN(), vertex_wire_ID, vertex_drift_t);
      gr_primary_vertex[i_p]->SetPoint(gr_primary_vertex[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
    }
    cout << "Primary vertex done" << endl;

    //bucle en las 3 vistas para el vértice true
    for(int i_p = 0; i_p < 3; i_p++) {
      double vertex_drift_t = slice->primary_vertex_true.vertex_wire.at(i_p).drift_t;
      int vertex_wire_ID = slice->primary_vertex_true.vertex_wire.at(i_p).Wire_ID;
      gr_primary_vertex_truth[i_p]->SetPoint(gr_primary_vertex_truth[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
      gr_primary_vertex_truth[i_p]->SetPoint(gr_primary_vertex_truth[i_p]->GetN(), vertex_wire_ID, vertex_drift_t);
      gr_primary_vertex_truth[i_p]->SetPoint(gr_primary_vertex_truth[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
    }
    cout << "true vertex done" << endl;

    //bucle en las 3 vistas para vertices secundarios true
    for(int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      for(int i_p = 0; i_p < 3; i_p++) {
        double vertex_drift_t = slice->vertex_vec.at(i_v).vertex_wire.at(i_p).drift_t;
        int vertex_wire_ID = slice->vertex_vec.at(i_v).vertex_wire.at(i_p).Wire_ID;
        gr_secondary_vertex[i_p]->SetPoint(gr_secondary_vertex[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
        gr_secondary_vertex[i_p]->SetPoint(gr_secondary_vertex[i_p]->GetN(), vertex_wire_ID, vertex_drift_t);
        gr_secondary_vertex[i_p]->SetPoint(gr_secondary_vertex[i_p]->GetN(), vertex_wire_ID + 0.0001, vertex_drift_t);
      }
    }
    cout << "secondary vertex done" << endl;

    //Dibuja las 3 vistas
    for(int i_p = 0; i_p < 3; i_p++) {
      c[i_p]->cd();

      gr_general_hits[i_p]->SetTitle(";Wire ID; Time [#mu s]");
      gr_general_hits[i_p]->SetMarkerSize(0.5);
      gr_general_hits[i_p]->SetMarkerColor(kGray+2);

      gr_muon_track_hits[i_p]->SetMarkerSize(0.5);
      gr_muon_track_hits[i_p]->SetMarkerColor(kViolet+2);
      gr_muon_track_hits[i_p]->SetMarkerStyle(8);

      gr_pion_track_hits[i_p]->SetMarkerSize(0.5);
      gr_pion_track_hits[i_p]->SetMarkerColor(kViolet+2);
      gr_pion_track_hits[i_p]->SetMarkerStyle(8);

      gr_proton_track_hits[i_p]->SetMarkerSize(0.5);
      gr_proton_track_hits[i_p]->SetMarkerColor(kViolet+2);
      gr_proton_track_hits[i_p]->SetMarkerStyle(8);

      gr_other_track_hits[i_p]->SetMarkerSize(0.5);
      gr_other_track_hits[i_p]->SetMarkerColor(kViolet+2);
      gr_other_track_hits[i_p]->SetMarkerStyle(8);


      gr_muon_primary_track_hits[i_p]->SetMarkerSize(0.5);
      gr_muon_primary_track_hits[i_p]->SetMarkerColor(kGreen+2);
      gr_muon_primary_track_hits[i_p]->SetMarkerStyle(8);

      gr_pion_primary_track_hits[i_p]->SetMarkerSize(0.5);
      gr_pion_primary_track_hits[i_p]->SetMarkerColor(kRed+2);
      gr_pion_primary_track_hits[i_p]->SetMarkerStyle(8);

      gr_proton_primary_track_hits[i_p]->SetMarkerSize(0.5);
      gr_proton_primary_track_hits[i_p]->SetMarkerColor(kBlue+2);
      gr_proton_primary_track_hits[i_p]->SetMarkerStyle(8);

      gr_other_primary_track_hits[i_p]->SetMarkerSize(0.5);
      gr_other_primary_track_hits[i_p]->SetMarkerColor(kOrange+2);
      gr_other_primary_track_hits[i_p]->SetMarkerStyle(8);

      gr_primary_vertex[i_p]->SetMarkerSize(1);
      gr_primary_vertex[i_p]->SetMarkerColor(kBlack);
      gr_primary_vertex[i_p]->SetMarkerStyle(8);

      gr_primary_vertex_truth[i_p]->SetMarkerSize(1);
      gr_primary_vertex_truth[i_p]->SetMarkerColor(kBlue +2 );
      gr_primary_vertex_truth[i_p]->SetMarkerStyle(22);

      gr_secondary_vertex[i_p]->SetMarkerSize(1);
      gr_secondary_vertex[i_p]->SetMarkerColor(kBlack);
      gr_secondary_vertex[i_p]->SetMarkerStyle(34);



      TMultiGraph *mg = new TMultiGraph();
      if(gr_general_hits[i_p]->GetN() > 0)mg->Add(gr_general_hits[i_p],"AP");
      if(gr_muon_track_hits[i_p]->GetN() > 0)mg->Add(gr_muon_track_hits[i_p],"AP");
      if(gr_pion_track_hits[i_p]->GetN() > 0)mg->Add(gr_pion_track_hits[i_p],"AP");
      if(gr_proton_track_hits[i_p]->GetN() > 0)mg->Add(gr_proton_track_hits[i_p],"AP");
      if(gr_other_track_hits[i_p]->GetN() > 0)mg->Add(gr_other_track_hits[i_p],"AP");
      if(gr_muon_primary_track_hits[i_p]->GetN() > 0)mg->Add(gr_muon_primary_track_hits[i_p],"AP");
      if(gr_pion_primary_track_hits[i_p]->GetN() > 0)mg->Add(gr_pion_primary_track_hits[i_p],"AP");
      if(gr_proton_primary_track_hits[i_p]->GetN() > 0)mg->Add(gr_proton_primary_track_hits[i_p],"AP");
      if(gr_other_primary_track_hits[i_p]->GetN() > 0)mg->Add(gr_other_primary_track_hits[i_p],"AP");
      if(gr_primary_vertex_truth[i_p]->GetN() > 0)mg->Add(gr_primary_vertex_truth[i_p],"AP");
      if(gr_primary_vertex[i_p]->GetN() > 0)mg->Add(gr_primary_vertex[i_p],"AP");
      if(gr_secondary_vertex[i_p]->GetN() > 0)mg->Add(gr_secondary_vertex[i_p],"AP");

      mg->SetTitle(";Wire ID; Time [#mu s]");
      mg->Draw("A");

      c[i_p]->Update();
    }

    //Lo mimso pero en 3D

    double max_sp_x=0;
    double max_sp_y=0;
    double max_sp_z=0;
    double min_sp_x=2000;
    double min_sp_y=2000;
    double min_sp_z=2000;


    for(int i_sp = 0; i_sp < slice->space_point_vec.size(); i_sp++) {
      SpacePoint sp = slice->space_point_vec.at(i_sp);
      gr_general_sp_2D->SetPoint(gr_general_sp_2D->GetN(), sp.z, sp.x, sp.y);
      if(sp.associated_pfp_ID == -1) continue;

      if(max_sp_x > sp.x) max_sp_x = sp.x;
      if(max_sp_y > sp.y) max_sp_y = sp.y;
      if(max_sp_z > sp.z) max_sp_z = sp.z;
      if(min_sp_x < sp.x) min_sp_x = sp.x;
      if(min_sp_y < sp.y) min_sp_y = sp.y;
      if(min_sp_z < sp.z) min_sp_z = sp.z;

      //Search for the pfp that the hit is associated to:
      int particle_index = 0;
      for(int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
        if(slice->pandora_particle.at(i_p).ID == sp.associated_pfp_ID) particle_index = i_p;
      }
      if(slice->pandora_particle.at( particle_index).is_pandora_primary) {
        if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 211) {
          gr_pion_primary_track_sp_2D->SetPoint(gr_pion_primary_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 13) {
          gr_muon_primary_track_sp_2D->SetPoint(gr_muon_primary_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 2212) {
          gr_proton_primary_track_sp_2D->SetPoint(gr_proton_primary_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else {
          gr_other_primary_track_sp_2D->SetPoint(gr_other_primary_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        }
      } else {
        if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 211) {
          gr_pion_track_sp_2D->SetPoint(gr_pion_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 13) {
          gr_muon_track_sp_2D->SetPoint(gr_muon_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else if(abs(slice->pandora_particle.at( particle_index).matched_pdg) == 2212) {
          gr_proton_track_sp_2D->SetPoint(gr_proton_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        } else {
          gr_other_track_sp_2D->SetPoint(gr_other_track_sp_2D->GetN(), sp.z, sp.x, sp.y);
        }
      }

    }

    cout << "sp done" << endl;
    //bucle en las 3 vistas para el vértice primario
    gr_primary_vertex_2D->SetPoint(gr_primary_vertex_2D->GetN(),slice->primary_vertex_reco.vertex_cm.Z(),slice->primary_vertex_reco.vertex_cm.X(),slice->primary_vertex_reco.vertex_cm.Y());
    cout << "Primary vertex done" << endl;

    gr_primary_vertex_truth_2D->SetPoint(gr_primary_vertex_2D->GetN(),slice->primary_vertex_true.vertex_cm.Z(),slice->primary_vertex_true.vertex_cm.X(),slice->primary_vertex_true.vertex_cm.Y());
    cout << "True vertex done" << endl;

    //bucle en las 3 vistas para vertices secundarios true
    for(int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      gr_secondary_vertex_2D->SetPoint(gr_primary_vertex_2D->GetN(),slice->vertex_vec.at(i_v).vertex_cm.Z(),slice->vertex_vec.at(i_v).vertex_cm.X(),slice->vertex_vec.at(i_v).vertex_cm.Y());
    }
    cout << "secondary vertex done" << endl;

    c2D->cd();
    gr_general_sp_2D->SetTitle(";z [cm]; y [cm];x [cm]");
    gr_general_sp_2D->SetMarkerSize(0.5);
    gr_general_sp_2D->SetMarkerColor(kViolet+2);

    gr_muon_track_sp_2D->SetMarkerSize(0.5);
    gr_muon_track_sp_2D->SetMarkerColor(kViolet+2);
    gr_muon_track_sp_2D->SetMarkerStyle(8);

    gr_pion_track_sp_2D->SetMarkerSize(0.5);
    gr_pion_track_sp_2D->SetMarkerColor(kViolet+2);
    gr_pion_track_sp_2D->SetMarkerStyle(8);

    gr_proton_track_sp_2D->SetMarkerSize(0.5);
    gr_proton_track_sp_2D->SetMarkerColor(kViolet+2);
    gr_proton_track_sp_2D->SetMarkerStyle(8);

    gr_other_track_sp_2D->SetMarkerSize(0.5);
    gr_other_track_sp_2D->SetMarkerColor(kViolet+2);
    gr_other_track_sp_2D->SetMarkerStyle(8);

    gr_muon_primary_track_sp_2D->SetMarkerSize(0.5);
    gr_muon_primary_track_sp_2D->SetMarkerColor(kGreen+2);
    gr_muon_primary_track_sp_2D->SetMarkerStyle(8);

    gr_pion_primary_track_sp_2D->SetMarkerSize(0.5);
    gr_pion_primary_track_sp_2D->SetMarkerColor(kRed+2);
    gr_pion_primary_track_sp_2D->SetMarkerStyle(8);

    gr_proton_primary_track_sp_2D->SetMarkerSize(0.5);
    gr_proton_primary_track_sp_2D->SetMarkerColor(kBlue+2);
    gr_proton_primary_track_sp_2D->SetMarkerStyle(8);

    gr_other_primary_track_sp_2D->SetMarkerSize(0.5);
    gr_other_primary_track_sp_2D->SetMarkerColor(kOrange+2);
    gr_other_primary_track_sp_2D->SetMarkerStyle(8);

    gr_primary_vertex_2D->SetMarkerSize(1);
    gr_primary_vertex_2D->SetMarkerColor(kBlack);
    gr_primary_vertex_2D->SetMarkerStyle(8);

    gr_primary_vertex_truth_2D->SetMarkerSize(1);
    gr_primary_vertex_truth_2D->SetMarkerColor(kBlue +2 );
    gr_primary_vertex_truth_2D->SetMarkerStyle(22);

    gr_secondary_vertex_2D->SetMarkerSize(1);
    gr_secondary_vertex_2D->SetMarkerColor(kBlack);
    gr_secondary_vertex_2D->SetMarkerStyle(34);

    //SetPoints2D(gr_general_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_muon_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_pion_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_proton_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_other_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_muon_primary_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_pion_primary_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_proton_primary_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_other_primary_track_sp_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_primary_vertex_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_primary_vertex_truth_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);
    //SetPoints2D(gr_secondary_vertex_2D,max_sp_x, max_sp_y, max_sp_z, min_sp_x, min_sp_y, min_sp_z);

    if(gr_general_sp_2D->GetN() > 0) gr_general_sp_2D->Draw("P");
    if(gr_muon_track_sp_2D->GetN() > 0) gr_muon_track_sp_2D->Draw("SAME P");
    if(gr_pion_track_sp_2D->GetN() > 0) gr_pion_track_sp_2D->Draw("SAME P");
    if(gr_proton_track_sp_2D->GetN() > 0) gr_proton_track_sp_2D->Draw("SAME P");
    if(gr_other_track_sp_2D->GetN() > 0) gr_other_track_sp_2D->Draw("SAME P");
    if(gr_muon_primary_track_sp_2D->GetN() > 0) gr_muon_primary_track_sp_2D->Draw("SAME P");
    if(gr_pion_primary_track_sp_2D->GetN() > 0) gr_pion_primary_track_sp_2D->Draw("SAME P");
    if(gr_proton_primary_track_sp_2D->GetN() > 0) gr_proton_primary_track_sp_2D->Draw("SAME P");
    if(gr_other_primary_track_sp_2D->GetN() > 0) gr_other_primary_track_sp_2D->Draw("SAME P");
    if(gr_primary_vertex_2D->GetN() > 0) gr_primary_vertex_2D->Draw("SAME P");
    if(gr_primary_vertex_truth_2D->GetN() > 0) gr_primary_vertex_truth_2D->Draw("SAME P");
    if(gr_secondary_vertex_2D->GetN() > 0) gr_secondary_vertex_2D->Draw("SAME P");
    c2D->Update();




    double nearest_distance = 1000;
    int nearest_vertex_index = 1231;
    //Look for the nearest vertex compatible vertex
    for(int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      if((slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag() < nearest_distance) {
        nearest_distance = (slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag();
        nearest_vertex_index = i_v;

      }
    }

    for(int i_v = 0; i_v < slice->vertex_vec.size(); i_v++) {
      if(i_v != nearest_vertex_index) continue;
      cout << slice->vertex_vec.at(i_v).vertex_cm.X() << " " << slice->vertex_vec.at(i_v).vertex_cm.Y() << " " << slice->vertex_vec.at(i_v).vertex_cm.Z() << endl;

      if((slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag() > 4) continue;

      TH1* h_FRANS_Like;
      TCanvas* c_display_FRANS =  new TCanvas(to_string(i_v).c_str(), "My Canvas", 400, 400);
      TText *t = new TText(.8,.8,to_string((slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag()).c_str());

      vector<TH1*> h_Track_FRANS_Like;
      vector<TCanvas*> c_Frans_Tracks;

      int target_plane_ID = 2;
      int histogram_low_edge = 0;
      int histogram_high_edge = 100;
      int num_signal_points = 100;
      TH1* h_total = new TH1D("h_total","h_total", num_signal_points, histogram_low_edge, histogram_high_edge);

      for(int i_p = 0; i_p < slice->pandora_particle.size();i_p++) {
        Reco_Particle p =  slice->pandora_particle.at(i_p);

        int num_bins = int(histogram_high_edge - histogram_low_edge);
        c_Frans_Tracks.push_back(new TCanvas(to_string(i_p + 10*i_v).c_str(), "My Canvas", 400, 400));
        cout << ( slice->vertex_vec.at(i_v).vertex_cm - p.track_start).Mag() << endl;
        cout << ( slice->vertex_vec.at(i_v).vertex_cm - p.track_end).Mag() << endl;

        if(p.pass_quality_cuts(cut_p, slice->hits)) {
        cout << p.track_end_wires.at(2).Wire_ID << " " << p.track_end_wires.at(2).drift_t << endl;
          if(( slice->vertex_vec.at(i_v).vertex_cm - p.track_start).Mag() < 4) {
            cout << "NICE" << endl;
            h_Track_FRANS_Like.push_back(slice->build_track_dQdx_hist(p.track_start_wires, 2, histogram_low_edge, histogram_high_edge, i_p, num_signal_points, p.ID));
            cout << "Integral: "<<  h_Track_FRANS_Like.at(i_p)->Integral() << endl;
          } else if(( slice->vertex_vec.at(i_v).vertex_cm - p.track_end).Mag() < 4) {
            cout << "NICE 2" << endl;
            h_Track_FRANS_Like.push_back(slice->build_track_dQdx_hist(p.track_end_wires, 2, histogram_low_edge, histogram_high_edge, i_p, num_signal_points, p.ID));
            cout << "Integral: "<<  h_Track_FRANS_Like.at(i_p)->Integral() << endl;
          } else {
            h_Track_FRANS_Like.push_back(new TH1D(("dQ" +to_string(i_p) ).c_str(), "dQ", num_bins, histogram_low_edge, histogram_high_edge));
          }
        } else {
          h_Track_FRANS_Like.push_back(new TH1D(("dQ" +to_string(i_p) ).c_str(), "dQ", num_bins, histogram_low_edge, histogram_high_edge));
        }
      }
      cout << "END VERTEXES" << endl;

      for(int i_p = 0; i_p < slice->pandora_particle.size();i_p++) {
        //if(!slice->pandora_particle.at(i_p).pass_quality_cuts(cut_p, slice->hits)) continue;
        c_Frans_Tracks.at(i_p)->cd();
        h_Track_FRANS_Like.at(i_p)->Draw("hist");
        h_Track_FRANS_Like.at(i_p)->SetLineColor(kGreen);
        c_Frans_Tracks.at(i_p)->Update();
      }
      cout << "END Tracks" << endl;

      h_FRANS_Like = slice->build_dQdx_hist(slice->vertex_vec.at(i_v).vertex_wire, 2 , 0, 100,i_v, 100);
      c_display_FRANS->cd();
      h_FRANS_Like->SetLineColor(kBlue +2 );
      if((slice->vertex_vec.at(i_v).vertex_cm - slice->primary_vertex_true.vertex_cm).Mag() > 4) h_FRANS_Like->SetLineColor(kRed+2);
      h_FRANS_Like->Draw("hist");
      t->Draw();
      c_display_FRANS->Update();

      cout << "END Display" << endl;

      button->WaitPrimitive();
      bool clicked = false;
      while(!clicked) {
        if(button->WaitPrimitive() == nullptr) {
          std::cout << "Canvas 1 clicked! Proceeding..." << std::endl;
          clicked = true;  // Exit the loop if the correct canvas is clicked
        }
      }

      for(int i_p = 0; i_p < slice->pandora_particle.size();i_p++) {
        c_Frans_Tracks.at(i_p)->Close();
        delete c_Frans_Tracks.at(i_p);
        delete h_Track_FRANS_Like.at(i_p);
      }
      c_Frans_Tracks.clear();
      h_Track_FRANS_Like.clear();

      delete h_FRANS_Like;
      c_display_FRANS->Close();
       delete c_display_FRANS;
    }


    bool clicked = false;
    while(!clicked) {
      if(button->WaitPrimitive() == nullptr) {
        std::cout << "Canvas 1 clicked! Proceeding..." << std::endl;
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }
    for(int i_p = 0; i_p < 3; i_p++) {
      int N_points = gr_general_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_general_hits[i_p]->RemovePoint(0);
      N_points = gr_pion_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_pion_track_hits[i_p]->RemovePoint(0);
      N_points = gr_muon_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_muon_track_hits[i_p]->RemovePoint(0);
      N_points = gr_proton_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_proton_track_hits[i_p]->RemovePoint(0);
      N_points = gr_other_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_other_track_hits[i_p]->RemovePoint(0);

      N_points = gr_pion_primary_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_pion_primary_track_hits[i_p]->RemovePoint(0);
      N_points = gr_muon_primary_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_muon_primary_track_hits[i_p]->RemovePoint(0);
      N_points = gr_proton_primary_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_proton_primary_track_hits[i_p]->RemovePoint(0);
      N_points = gr_other_primary_track_hits[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_other_primary_track_hits[i_p]->RemovePoint(0);


      N_points = gr_primary_vertex[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_primary_vertex[i_p]->RemovePoint(0);
      N_points = gr_primary_vertex_truth[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_primary_vertex_truth[i_p]->RemovePoint(0);
      N_points = gr_secondary_vertex[i_p]->GetN();
      for (int i_point = 0; i_point < N_points; i_point++) gr_secondary_vertex[i_p]->RemovePoint(0);
    }


    int N_points = gr_general_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_general_sp_2D->RemovePoint(0);
    N_points = gr_pion_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_pion_track_sp_2D->RemovePoint(0);
    N_points = gr_muon_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_muon_track_sp_2D->RemovePoint(0);
    N_points = gr_proton_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_proton_track_sp_2D->RemovePoint(0);
    N_points = gr_other_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_other_track_sp_2D->RemovePoint(0);

    N_points = gr_pion_primary_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_pion_primary_track_sp_2D->RemovePoint(0);
    N_points = gr_muon_primary_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_muon_primary_track_sp_2D->RemovePoint(0);
    N_points =  gr_proton_primary_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_proton_primary_track_sp_2D->RemovePoint(0);
    N_points = gr_other_primary_track_sp_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_other_primary_track_sp_2D->RemovePoint(0);

    N_points = gr_primary_vertex_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_primary_vertex_2D->RemovePoint(0);
    N_points = gr_primary_vertex_truth_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_primary_vertex_truth_2D->RemovePoint(0);
    N_points = gr_secondary_vertex_2D->GetN();
    for (int i_point = 0; i_point < N_points; i_point++) gr_secondary_vertex_2D->RemovePoint(0);

  }


  for (std::map<string,double>::iterator it=cases_cont.begin(); it!=cases_cont.end(); ++it)
    std::cout << it->first << ": " << it->second << '\n';


}