#include "../Includes.h"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}

void event_display()
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
  GenerateDictionaries();
  TTree *tree;
  TFile *input;


  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0);

  //Declare the variables
  Slice *slice = 0;

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_82p7k.root";
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

  TCanvas* c = new TCanvas("c1", "My Canvas", 1280, 1280);
  c->Divide(2,2);

  //std::vector<TCanvas*> c = {, new TCanvas("c2", "My Canvas", 500, 400), new TCanvas("c3", "My Canvas", 500, 400)};
  //c[0]->SetWindowPosition(10, 50);
  //c[1]->SetWindowPosition(700, 50);
  //c[2]->SetWindowPosition(10, 500);

  //TCanvas* c2D = new TCanvas("c4", "My Canvas", 500, 400);
  //c2D->SetWindowPosition(700, 500);

  TCanvas* button = new TCanvas("c12", "My Canvas", 200, 200);
  button->SetWindowPosition(1400, 50);
  button->SetFillColor(kRed);
  button->Draw();

  vector<string> cut_name = {"no_cut", "is_clear_cosmic_cut", "reco_cut" , "fv_cut",  "crumbs_cut", "track_cut",
                             "razzled_muon_like_cut", "shower_cut", "at_most_one_pfp_scapping_cut", "proton_rejection_razzled"};

  const int num_cuts = cut_name.size();
  Cut_Parameters cut_p;
  cut_p.use_default_cut_p();


  //pp = Primary Particle
  map<string,double> cases_cont;

  for(int i_e = 0; i_e < nEntries; ++i_e) {
    tree->GetEntry(i_e);
    if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

    //if(slice->slice_ID != 0) continue;

    if(slice->true_interaction.is_selected_final_state("CC1Pi", 0.325)) continue;
    if(!slice->true_interaction.is_selected_background_final_state("CC_mu_0pi_1p", 0.325)) continue;

    bool pass_cut = true;
    //Plot only if correctly reconstructed
    for(int i_cut = 0; i_cut < num_cuts; i_cut++) {
      if(!pass_cut) continue;
      if(!slice->pass_analysis_cut(cut_p, cut_name[i_cut])) {
        pass_cut = false;
      }
    }
    if(!pass_cut) continue;

    //Check if there is an extra slice or not
    bool extra_slice = false;

    double run_Id = slice->run_ID;
    double subrun_Id = slice->subrun_ID;
    double event_Id = slice->event_ID;

    tree->GetEntry(i_e-1);
    if((run_Id == slice->run_ID) && (subrun_Id == slice->subrun_ID) && (event_Id == slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e+1);
    if((run_Id == slice->run_ID) && (subrun_Id == slice->subrun_ID) && (event_Id == slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e-2);
    if((run_Id == slice->run_ID) && (subrun_Id == slice->subrun_ID) && (event_Id == slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e+2);
    if((run_Id == slice->run_ID) && (subrun_Id == slice->subrun_ID) && (event_Id == slice->event_ID))extra_slice = true;
    tree->GetEntry(i_e);

    //Search for pfps in slice
    vector<int> pfp_ID_vector;
    for(Reco_Particle slice_pfp: slice->pandora_particle) {
      pfp_ID_vector.push_back(slice_pfp.ID);
    }

    bool has_activity_near_true_vertex = false;
    for(int i_sp = 0; i_sp < slice->space_point_vec.size(); i_sp++) {
      SpacePoint sp = slice->space_point_vec.at(i_sp);
      if(find(pfp_ID_vector.begin(), pfp_ID_vector.end(), sp.associated_pfp_ID) == pfp_ID_vector.end()) continue;
      double delta_x = sp.x - slice->primary_vertex_true.vertex_cm.x();
      double delta_y = sp.y - slice->primary_vertex_true.vertex_cm.y();
      double delta_z = sp.z - slice->primary_vertex_true.vertex_cm.z();

      if(sqrt(pow(delta_x,2) + sqrt(pow(delta_y,2)) +sqrt(pow(delta_z,2))) < 10) has_activity_near_true_vertex = true;
    }

    bool skip = false;
    if(extra_slice && !has_activity_near_true_vertex) {
      skip = true;
    }




    bool show = true;
    if(slice->get_num_showers(cut_p, true, false) > 0) show = true;



    if(skip) continue;
    if(!show) continue;
    slice->print(cut_p, true, true, true, true,false, true);
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
      c->cd(i_p + 1);
      string plane = "Collection plane";
      if(i_p == 0) plane = "Induction plane 1";
      if(i_p == 1) plane = "Induction plane 2";

      gr_general_hits[i_p]->SetTitle((plane + ";Wire ID; Time [#mu s]").c_str());
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

      mg->SetTitle((plane + ";Wire ID; Time [#mu s]").c_str());
      mg->Draw("A");

      c->Update();
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

    c->cd(4);
    gr_general_sp_2D->SetTitle("3D view;z [cm]; y [cm];x [cm]");
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
    c->Update();

    button->WaitPrimitive();

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