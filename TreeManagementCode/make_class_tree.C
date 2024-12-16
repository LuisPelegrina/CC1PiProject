#include "tree_utils.cpp"
#include "../Includes.h"
bool debug = false;

void make_class_tree()
{
  GenerateDictionaries();
  //Definition of the file from which the tree is readed and the folder inside the ".root" file where the tree is stored
  TFile *input_file;
  TDirectoryFile *tree_dir;

  //Definition of the tree that stores information about each LArSoft event
  TTree *event_tree;

  //Definition of the tree that stores information about each LArSoft subrun (number of events generated in the subrun, POT...)
  TTree *subrun_tree;

  //Definition of the path where the file containing the truth tree is stored, along with the code that opens the file and both the subrun and the event by event tree
  string file_name = "82p7k";
  string path_to_tree = "../Data/analysis_output_" + file_name + ".root";
  input_file = new TFile(path_to_tree.c_str());
  tree_dir = (TDirectoryFile*)input_file->Get("ana");
  event_tree =(TTree*)tree_dir->Get("tree");
  subrun_tree =(TTree*)tree_dir->Get("subrun_tree");

  //Code that sets the Branches fot the subrun and event trees in order to use the information they store, please visit "tree_utils.cpp" to see how the functions are defined
  set_branch(event_tree);
  set_branch_subtree(subrun_tree);

  bool has_gen = true;
  if(event_tree->GetBranch("gen_part_start_pos_X") == NULL) has_gen = false;
  bool has_g4 = true;
  if(event_tree->GetBranch("g4_part_trackID") == NULL) has_g4 = false;
  bool has_reco = true;
  if(event_tree->GetBranch("slc_ID") == NULL) has_g4 = false;
  bool has_sp = true;
  if(event_tree->GetBranch("slc_sp_x") == NULL) has_sp = false;
  bool has_hits = true;
  if(event_tree->GetBranch("slc_hits_wire_ID") == NULL) has_hits = false;
  bool has_calo = true;
  if(event_tree->GetBranch("slc_pfp_calo_dEdx") == NULL) has_calo = false;

  //POT NORMALIZATION
  double spill_POT = 5*pow(10,12);
  double expected_POT = pow(10,21);
  double num_months = 2;

  expected_POT = expected_POT*num_months/(9*3);
  double sim_POT = 0;

  int num_subrun_entries = subrun_tree->GetEntries();
  cout << "Total number of subruns: " << num_subrun_entries << endl;
  for(int i_e = 0; i_e < num_subrun_entries; i_e++) {
    subrun_tree->GetEntry(i_e);
    sim_POT += POT;
  }
  cout << "Simulated POT: " << sim_POT << endl;


  Slice tree_slice;

  TFile* output_file_nu;
  output_file_nu = TFile::Open(("../Data/processed_data/processed_data_" + file_name + ".root").c_str(),"RECREATE");
  TTree *outputTree_nu = new TTree("tree", "Collection of Events" );
  outputTree_nu->Branch("slice", "Slice", &tree_slice, 3200, 2);

  //Example of how the event_tree is used. In this case the number of entries in the tree is readed. After that a loop prints on screen all the information in an human-readable way. Please refer to "tree_utils.cpp" for a definition of each variable printed.
  //Calculate the number of entries of the event_tree
  int num_entries = event_tree->GetEntries();
  //Print the number of events on the truth_tree
  cout << "Total number of events: " << num_entries << endl;
  //Start the loop on each event (each event_tree entry)
  for(int i_e = 0; i_e < num_entries; i_e++) {
    tree_slice = Slice();
    //Start reading the entry with index i_e
    if(debug) cout << "Entry: " << i_e << "started" << endl;
    event_tree->GetEntry(i_e);
    if(i_e%100 == 0)cout << i_e << endl;
    if(debug) cout << "Entry: " << i_e << "passed get" << endl;
    //cout << i_e << endl;
    tree_slice.target_months = num_months;
    tree_slice.event_ID = event_ID;
    tree_slice.run_ID = run_ID;
    tree_slice.subrun_ID = subrun_ID;
    tree_slice.file_name = *reco2_file_name;
    tree_slice.gen_index = gen_index;
    if(debug) cout << "event info done" << endl;

    if(nu_weight == -1) nu_weight = 0;
    tree_slice.set_weight(nu_weight, expected_POT, sim_POT);

    //Fill true interaction
    tree_slice.true_interaction.E0_nu = nu_E0;
    tree_slice.true_interaction.PDG_nu = nu_PDG;
    tree_slice.true_interaction.interaction_mode = nu_interaction_mode;
    tree_slice.true_interaction.interaction_type = nu_interaction_type;
    tree_slice.true_interaction.CCNC = nu_CC_NC;
    tree_slice.true_interaction.target = nu_target;
    tree_slice.true_interaction.HitNuc = nu_HitNuc;
    tree_slice.true_interaction.HitQuark= nu_HitQuark;
    tree_slice.true_interaction.W = nu_W;
    tree_slice.true_interaction.X = nu_X;
    tree_slice.true_interaction.Y = nu_Y;
    tree_slice.true_interaction.QSqr = nu_QSqr;
    if(debug) cout << "true interaction done" << endl;

     if(has_gen) {
      tree_slice.true_interaction.primary_vertex = TVector3(gen_part_start_pos_X->at(0), gen_part_start_pos_Y->at(0), gen_part_start_pos_Z->at(0));
      //Fill gen particles
      vector<GeneratorParticle> gen_p_vec;
      vector<GeneratorParticle> gen_primary_p_vec;
      for(int i_p = 0; i_p < gen_part_trackID->size(); i_p++) {
        GeneratorParticle gen_p;
        gen_p.set_gen_particle(gen_part_trackID->at(i_p), gen_part_mother->at(i_p), gen_part_PDGcode->at(i_p), gen_part_mass->at(i_p), gen_part_E0->at(i_p),
                               TVector3(gen_part_start_pos_X->at(i_p), gen_part_start_pos_Y->at(i_p), gen_part_start_pos_Z->at(i_p)),
                               TVector3(gen_part_P0_X->at(i_p), gen_part_P0_Y->at(i_p), gen_part_P0_Z->at(i_p)), gen_part_statusCode->at(i_p));

        gen_p_vec.push_back(gen_p);
        if(gen_p.is_primary()) gen_primary_p_vec.push_back(gen_p);
      }
      tree_slice.true_interaction.gen_particles = gen_p_vec;
      tree_slice.true_interaction.gen_primary_particles = gen_primary_p_vec;
    }
    if(debug) cout << "gen part done" << endl;

    if(has_g4) {

      //Fill g4 particles
      vector<G4Particle> g4_p_vec;
      for(int i_p = 0; i_p < g4_part_trackID->size(); i_p++) {
        G4Particle g4_p;
        g4_p.set_g4_particle(g4_part_trackID->at(i_p), g4_part_mother->at(i_p), g4_part_PDGcode->at(i_p),
                             g4_part_mass->at(i_p), g4_part_E0->at(i_p),
                             TVector3(g4_part_start_pos_X->at(i_p), g4_part_start_pos_Y->at(i_p), g4_part_start_pos_Z->at(i_p)),
                             TVector3(g4_part_P0_X->at(i_p), g4_part_P0_Y->at(i_p), g4_part_P0_Z->at(i_p)),
                             g4_part_Ef->at(i_p),
                             TVector3(g4_part_end_pos_X->at(i_p), g4_part_end_pos_Y->at(i_p), g4_part_end_pos_Z->at(i_p)),
                             TVector3(g4_part_Pf_X->at(i_p), g4_part_Pf_Y->at(i_p), g4_part_Pf_Z->at(i_p)),
                             g4_part_start_T->at(i_p), g4_part_end_T->at(i_p), g4_part_TL->at(i_p),
                             g4_part_process->at(i_p),
                             g4_part_end_process->at(i_p));

        g4_p_vec.push_back(g4_p);
      }
      tree_slice.true_interaction.g4_particles = g4_p_vec;
    }
    if(debug) cout << "g4 part done" << endl;

    if(has_reco) {
      tree_slice.slice_ID = slc_ID;
      tree_slice.is_reconstructed = slc_is_reconstructed;
      tree_slice.is_clear_cosmic = slc_is_obvious_cosmic;
      if(debug) cout << "withou reco done" << endl;

      //if(slc_is_reconstructed && !slc_is_obvious_cosmic) {
        tree_slice.truth_matching_purity = slc_truth_matching_purity;
        tree_slice.truth_matching_completeness = slc_truth_matching_completeness;
        tree_slice.crumbs_score = slc_crumbs_score;
        tree_slice.nu_score = slc_nu_score;
        if(debug) cout << "Slice basic done" << endl;

        //Get the true vertex
        tree_slice.primary_vertex_true.vertex_cm = TVector3(slc_true_prim_vtx_x, slc_true_prim_vtx_y, slc_true_prim_vtx_z);
        if(has_hits) {
          tree_slice.primary_vertex_true.vertex_wire.at(0).Wire_ID = slc_true_prim_vtx_wires_U;
          tree_slice.primary_vertex_true.vertex_wire.at(0).Channel_ID = slc_true_prim_vtx_channel_U;
          tree_slice.primary_vertex_true.vertex_wire.at(0).drift_t = slc_true_prim_vtx_time_tick;
          tree_slice.primary_vertex_true.vertex_wire.at(1).Wire_ID = slc_true_prim_vtx_wires_V;
          tree_slice.primary_vertex_true.vertex_wire.at(1).Channel_ID = slc_true_prim_vtx_channel_V;
          tree_slice.primary_vertex_true.vertex_wire.at(1).drift_t = slc_true_prim_vtx_time_tick;
          tree_slice.primary_vertex_true.vertex_wire.at(2).Wire_ID = slc_true_prim_vtx_wires_C;
          tree_slice.primary_vertex_true.vertex_wire.at(2).Channel_ID = slc_true_prim_vtx_channel_C;
          tree_slice.primary_vertex_true.vertex_wire.at(2).drift_t = slc_true_prim_vtx_time_tick;
        }
        if(debug) cout << "true vertex done" << endl;

        //Get the primary vertex
        tree_slice.primary_vertex_reco.vertex_cm = TVector3(slc_prim_vtx_x, slc_prim_vtx_y, slc_prim_vtx_z);
        if(has_hits) {
          tree_slice.primary_vertex_reco.vertex_wire.at(0).Wire_ID = slc_prim_vtx_wires_U;
          tree_slice.primary_vertex_reco.vertex_wire.at(0).Channel_ID = slc_prim_vtx_channel_U;
          tree_slice.primary_vertex_reco.vertex_wire.at(0).drift_t = slc_prim_vtx_time_tick;
          tree_slice.primary_vertex_reco.vertex_wire.at(1).Wire_ID = slc_prim_vtx_wires_V;
          tree_slice.primary_vertex_reco.vertex_wire.at(1).Channel_ID = slc_prim_vtx_channel_V;
          tree_slice.primary_vertex_reco.vertex_wire.at(1).drift_t = slc_prim_vtx_time_tick;
          tree_slice.primary_vertex_reco.vertex_wire.at(2).Wire_ID = slc_prim_vtx_wires_C;
          tree_slice.primary_vertex_reco.vertex_wire.at(2).Channel_ID = slc_prim_vtx_channel_C;
          tree_slice.primary_vertex_reco.vertex_wire.at(2).drift_t = slc_prim_vtx_time_tick;
        }
        if(debug) cout << "primary vertex done" << endl;


        //Set T0 parameters
        tree_slice.t0 = slc_t0;
        tree_slice.t0_score = slc_t0_score;
        tree_slice.OpT0_fraction = slc_OpT0Fraction;
        if(debug) cout << "t0 done" << endl;

        //Set Space Points
        if(has_sp) {
          vector<SpacePoint> sp_vec;
          for(int i_sp = 0; i_sp < slc_sp_x->size();i_sp++) {
            SpacePoint sp;
            sp.set_space_point(slc_sp_x->at(i_sp), slc_sp_y->at(i_sp), slc_sp_z->at(i_sp), slc_sp_integral->at(i_sp), slc_sp_sigma_integral->at(i_sp), slc_sp_associated_pfp_ID->at(i_sp));
            sp_vec.push_back(sp);
          }
          tree_slice.space_point_vec = sp_vec;
        }
        if(debug) cout << "Space points done" << endl;

        //Set Hits
        if(has_hits) {
          vector<Hit> hit_vec;
          for(int i_h = 0; i_h < slc_hits_wire_ID->size();i_h++) {
            Hit hit;
            hit.ID = i_h;
            hit.set_hit(slc_hits_integral->at(i_h), slc_hits_sigma_integral->at(i_h), slc_hits_peaktime->at(i_h), slc_hits_RMS->at(i_h),  slc_hits_TPC_ID->at(i_h),  slc_hits_plane_ID->at(i_h),  slc_hits_wire_ID->at(i_h),  slc_hits_channel_ID->at(i_h), slc_hits_associated_pfp_ID->at(i_h));
            hit_vec.push_back(hit);
          }
          tree_slice.hits = hit_vec;
        }
        if(debug) cout << "Hits done" << endl;

        //Set the primary corrected vertex
        if(has_hits) tree_slice.set_corrected_vertex_wire();
        if(debug) cout << "vertex wire done" << endl;

        //Set the vertex vector
        vector<Vertex> vertex_vec;
        for(int i_p = 0; i_p < slc_pfp_ID->size(); i_p++) {
          Vertex v;
          v.vertex_cm = TVector3(slc_pfp_vtx_x->at(i_p), slc_pfp_vtx_y->at(i_p), slc_pfp_vtx_z->at(i_p));
          if(has_hits) {
            v.vertex_wire.at(0).Wire_ID = slc_pfp_vtx_wires_U->at(i_p);
            v.vertex_wire.at(0).Channel_ID = slc_pfp_vtx_channel_U->at(i_p);
            v.vertex_wire.at(0).drift_t = slc_pfp_vtx_time_tick->at(i_p);
            v.vertex_wire.at(1).Wire_ID = slc_pfp_vtx_wires_V->at(i_p);
            v.vertex_wire.at(1).Channel_ID = slc_pfp_vtx_channel_V->at(i_p);
            v.vertex_wire.at(1).drift_t = slc_pfp_vtx_time_tick->at(i_p);
            v.vertex_wire.at(2).Wire_ID = slc_pfp_vtx_wires_C->at(i_p);
            v.vertex_wire.at(2).Channel_ID = slc_pfp_vtx_channel_C->at(i_p);
            v.vertex_wire.at(2).drift_t = slc_pfp_vtx_time_tick->at(i_p);
          }
          vertex_vec.push_back(v);
        }
        tree_slice.vertex_vec = vertex_vec;
        if(debug) cout << "vertex vec done" << endl;

        //Set the pfps
        vector<Reco_Particle> pandora_p_vec;
        for(int i_p = 0; i_p < slc_pfp_ID->size(); i_p++) {
          Reco_Particle pandora_p;
          pandora_p.matched_pdg = slc_pfp_true_pdg->at(i_p);
          pandora_p.true_track_id = slc_pfp_true_trackid->at(i_p);
          pandora_p.true_energy = slc_pfp_true_energy->at(i_p);
          pandora_p.purity = slc_pfp_purity->at(i_p);
          pandora_p.completeness = slc_pfp_completeness->at(i_p);
          pandora_p.true_P = TVector3(slc_pfp_true_p_x->at(i_p) ,slc_pfp_true_p_y->at(i_p), slc_pfp_true_p_z->at(i_p));
          pandora_p.ID = slc_pfp_ID->at(i_p);
          pandora_p.parent_ID = slc_pfp_parent_ID->at(i_p);
          pandora_p.is_pandora_primary = slc_pfp_mother_is_primary->at(i_p);
          pandora_p.track_score = slc_pfp_track_score->at(i_p);
          pandora_p.has_track = slc_pfp_has_track->at(i_p);
          pandora_p.pandora_PDG = slc_pfp_pandora_PDG->at(i_p);

          Razzled razzled;
          razzled.muon_score = slc_pfp_razzled_score_mu->at(i_p);
          razzled.pion_score = slc_pfp_razzled_score_pi->at(i_p);
          razzled.proton_score =slc_pfp_razzled_score_p->at(i_p);
          razzled.photon_score =slc_pfp_razzled_score_pho->at(i_p);
          razzled.electron_score = slc_pfp_razzled_score_e->at(i_p);
          razzled.razzled_pdg = slc_pfp_razzled_pdg->at(i_p);
          pandora_p.razzled_score = razzled;

          pandora_p.track_lenght = slc_pfp_track_length->at(i_p);
          pandora_p.track_kinetic_energy = slc_pfp_track_kinetic_energy->at(i_p);
          pandora_p.track_visible_energy = slc_pfp_track_visible_energy->at(i_p);
          pandora_p.track_charge = slc_pfp_track_charge->at(i_p);
          pandora_p.start_direction = TVector3(slc_pfp_track_dir_x->at(i_p), slc_pfp_track_dir_y->at(i_p), slc_pfp_track_dir_z->at(i_p));
          pandora_p.track_charge = slc_pfp_track_charge->at(i_p);
          pandora_p.track_theta = slc_pfp_track_theta->at(i_p);
          pandora_p.track_phi = slc_pfp_track_phi->at(i_p);

          pandora_p.shower_direction = TVector3(slc_pfp_shower_dir_x->at(i_p), slc_pfp_shower_dir_y->at(i_p), slc_pfp_shower_dir_z->at(i_p));
          pandora_p.shower_length = slc_pfp_shower_length->at(i_p);
          pandora_p.shower_opening_angle = slc_pfp_shower_opening_angle->at(i_p);
          pandora_p.shower_energy = slc_pfp_shower_energy->at(i_p);
          pandora_p.shower_dEdx = slc_pfp_shower_dEdx->at(i_p);

          if(has_calo) {
            //Set calo
            vector<CaloPoint> calo_vec;
            for(int i_c = 0; i_c < slc_pfp_calo_dEdx->size();i_c++) {
              CaloPoint calo;
              calo.dEdx = slc_pfp_calo_dEdx->at(i_p).at(i_c);
              calo.dQdx =slc_pfp_calo_dEdx->at(i_p).at(i_c);
              calo.pitch = slc_pfp_calo_dEdx->at(i_p).at(i_c);
              calo.residual_range = slc_pfp_calo_dEdx->at(i_p).at(i_c);
              calo.x_pos = slc_pfp_calo_dEdx->at(i_p).at(i_c);
              calo.y_pos = slc_pfp_calo_dEdx->at(i_p).at(i_c);
              calo.z_pos = slc_pfp_calo_dEdx->at(i_p).at(i_c);
              calo_vec.push_back(calo);
            }
            pandora_p.calo_point_vec = calo_vec;
          }

          //Add track start/end
          pandora_p.track_start = TVector3(slc_pfp_trk_start_x->at(i_p), slc_pfp_trk_start_y->at(i_p), slc_pfp_trk_start_z->at(i_p));
          pandora_p.track_start_wires.at(0).Wire_ID = slc_pfp_trk_start_wires_U->at(i_p);
          pandora_p.track_start_wires.at(0).Channel_ID = slc_pfp_trk_start_channel_U->at(i_p);
          pandora_p.track_start_wires.at(0).drift_t = slc_pfp_trk_start_time_tick->at(i_p);
          pandora_p.track_start_wires.at(1).Wire_ID = slc_pfp_trk_start_wires_V->at(i_p);
          pandora_p.track_start_wires.at(1).Channel_ID = slc_pfp_trk_start_channel_V->at(i_p);
          pandora_p.track_start_wires.at(1).drift_t = slc_pfp_trk_start_time_tick->at(i_p);
          pandora_p.track_start_wires.at(2).Wire_ID = slc_pfp_trk_start_wires_C->at(i_p);
          pandora_p.track_start_wires.at(2).Channel_ID = slc_pfp_trk_start_channel_C->at(i_p);
          pandora_p.track_start_wires.at(2).drift_t = slc_pfp_trk_start_time_tick->at(i_p);

          pandora_p.track_end = TVector3(slc_pfp_trk_end_x->at(i_p), slc_pfp_trk_end_y->at(i_p), slc_pfp_trk_end_z->at(i_p));
          pandora_p.track_end_wires.at(0).Wire_ID = slc_pfp_trk_end_wires_U->at(i_p);
          pandora_p.track_end_wires.at(0).Channel_ID = slc_pfp_trk_end_channel_U->at(i_p);
          pandora_p.track_end_wires.at(0).drift_t = slc_pfp_trk_end_time_tick->at(i_p);
          pandora_p.track_end_wires.at(1).Wire_ID = slc_pfp_trk_end_wires_V->at(i_p);
          pandora_p.track_end_wires.at(1).Channel_ID = slc_pfp_trk_end_channel_V->at(i_p);
          pandora_p.track_end_wires.at(1).drift_t = slc_pfp_trk_end_time_tick->at(i_p);
          pandora_p.track_end_wires.at(2).Wire_ID = slc_pfp_trk_end_wires_C->at(i_p);
          pandora_p.track_end_wires.at(2).Channel_ID = slc_pfp_trk_end_channel_C->at(i_p);
          pandora_p.track_end_wires.at(2).drift_t = slc_pfp_trk_end_time_tick->at(i_p);
          //Track chi2
          pandora_p.chi2_score.muon_score = slc_pfp_track_chi2_muon->at(i_p);
          pandora_p.chi2_score.pion_score = slc_pfp_track_chi2_pion->at(i_p);
          pandora_p.chi2_score.proton_score = slc_pfp_track_chi2_proton->at(i_p);
          pandora_p.chi2_score.kaon_score = slc_pfp_track_chi2_kaon->at(i_p);

          pandora_p_vec.push_back(pandora_p);
        }
        tree_slice.pandora_particle = pandora_p_vec;
        if(debug) cout << "everything done" << endl;
      //}
    }

/*
    cout << tree_slice.event_ID << " " << tree_slice.subrun_ID << " " << tree_slice.run_ID << " " << i_e << " " << tree_slice.slice_ID <<  endl;
    cout << "Is cosmic: "<< tree_slice.is_clear_cosmic << " is_reconstructed: " << tree_slice.is_reconstructed<< endl;
    cout << "crumbs score: " << tree_slice.crumbs_score << endl;
    cout << "Pfparicles: " << endl;
    for (Reco_Particle reco_p: tree_slice.pandora_particle) {
      cout << reco_p.ID << " " << reco_p.matched_pdg <<  endl;
    }
    */


    outputTree_nu->Fill();
    if(debug) cout << "everything written" << endl;

  }


  output_file_nu->cd();
  outputTree_nu->Write();

  output_file_nu->Close();



}