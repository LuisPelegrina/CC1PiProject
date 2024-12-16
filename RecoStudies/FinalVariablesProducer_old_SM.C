#include "../Includes.h"
#include "tree_utils.cpp"

void FinalVariablesProducer_old_SM()
{ 
    GenerateDictionaries();

    double muon_m = 105.7;
    int num_trees = 5;
    string strRuta;

    TTree *tree[num_trees];
    TFile *input[num_trees];

    Cut_Parameters cut_p;

    cut_p.min_distance_to_wall_x_y = 20;
    cut_p.min_distance_to_last_z_wall = 20;
    cut_p.min_distance_to_first_z_wall = 30;
    cut_p.min_distance_to_CPA = 5;
    cut_p.min_crumbs_score = 0;
    cut_p.apply_quality_cuts = false;
    cut_p.max_primary_distance_to_vertex = 15;
    cut_p.min_track_lenght = 10;
    cut_p.min_track_score = 0.5;

  //Declare the variables
    Slice *slice = 0;

    double num_months = 9.*3;
    double norm_months = 9.*3;

    double  re_w = num_months/norm_months;

    TTree *output_tree[num_trees];
    TFile *output_file[num_trees];

    cout << "Files loading" << endl;

    num_trees = 5;
    for(int i_t = 0; i_t < num_trees; i_t++) {
        if(i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Data/processed_data/processed_data_nu_cosmics.root";
        if(i_t == 1) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Data/processed_data/processed_data_in_time_cosmics.root";
        if(i_t == 2) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Data/processed_data/processed_data_nu.root";
        if(i_t == 3) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Data/processed_data/processed_data_nu_ncc1pi.root";
        if(i_t == 4) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Data/processed_data/processed_data_nu_coh.root";
        input[i_t] = new TFile(strRuta.c_str());

        tree[i_t] =(TTree*)input[i_t]->Get("tree");
        tree[i_t]->SetBranchAddress("slice", &slice);
        string output_name;

        if(i_t == 0) output_name = "final_data/final_data_nu_cosmics.root";
        if(i_t == 1) output_name = "final_data/final_data_in_time_cosmics.root";
        if(i_t == 2) output_name = "final_data/final_data_nu.root";
        if(i_t == 3) output_name = "final_data/final_data_nu_ncc1pi.root";
        if(i_t == 4) output_name = "final_data/final_data_nu_coh.root";
        output_file[i_t] = TFile::Open(output_name.c_str(),"RECREATE");

        output_tree[i_t] = new TTree("tree", "Collection of Events" );
        Branch(output_tree[i_t]);

    }
    cout << "Files loading done" << endl;

    bool has_nan = false;
    bool has_not_defined = true;

    vector<string> final_states_particle = {"CC1Pi", "NC", "CC_e", "CC_mu_0pi_0p", "CC_mu_0pi_1p", "CC_mu_0pi_2p", "CC_mu_2pi"};


    for(int i_t = 0; i_t < num_trees; i_t++) {
        int num_entries = tree[i_t]->GetEntries();

        for(int i_e = 0; i_e < num_entries; ++i_e) {
            tree[i_t]->GetEntry(i_e);
            if(i_e%100 == 0) cout << "Entry:" << i_e << endl;

            gen_index = slice->gen_index;
            if(gen_index == 1) {
                input_nu_type = slice->true_interaction.get_nu_final_state();
                if(input_nu_type.compare("CC_mu") == 0) {
                    if(slice->true_interaction.is_selected_background_production_state("CC_mu_res")) input_nu_type = "CC_mu_res";
                    if(slice->true_interaction.is_selected_background_production_state("CC_mu_qe")) input_nu_type = "CC_mu_qe";
                    if(slice->true_interaction.is_selected_final_state("CC1Pi")) input_nu_type = "CC1Pi";
                    if(slice->true_interaction.is_selected_final_state("N_CC1Pi")) input_nu_type = "N_CC1Pi";
                    if(slice->true_interaction.is_selected_final_state("Coherent")) input_nu_type = "N_CC1Pi_Coherent";
                }
            } else {
                input_nu_type = "not_nu";
            }

            bool checked = false;
            if(gen_index == 1) {
                for (int i_fs = 0; i_fs < final_states_particle.size(); i_fs++) {
                    if(checked) continue;
                    string current_final_state = final_states_particle.at(i_fs);
                    if (!slice->true_interaction.is_selected_background_final_state(current_final_state)) continue;
                    input_final_state = current_final_state;
                    checked = true;
                }
            }


            weight = slice->weight*re_w;

            E_0 = slice->true_interaction.E0_nu;

            crumbs_score = slice->crumbs_score;
            is_clear_cosmic = slice->is_clear_cosmic;
            is_reconstructed = slice->is_reconstructed;

            v_x = slice->primary_vertex_reco.vertex_cm.X();
            v_y = slice->primary_vertex_reco.vertex_cm.Y();
            v_z = slice->primary_vertex_reco.vertex_cm.Z();

            v_x_true = slice->primary_vertex_true.vertex_cm.X();
            v_y_true = slice->primary_vertex_true.vertex_cm.Y();
            v_z_true = slice->primary_vertex_true.vertex_cm.Z();

            num_primary_tracks = slice->get_num_tracks(cut_p , true, false);
            num_primary_showers = slice->get_num_showers(cut_p , true, false);

            num_razzled_primary_photons = slice->get_num_razzled_particle(22,cut_p, true, false);
            num_razzled_primary_electrons = slice->get_num_razzled_particle(11,cut_p, true, false);
            num_razzled_primary_protons = slice->get_num_razzled_particle(2212,cut_p, true, false);

            num_razzled_primary_muon_like =  slice->get_num_razzled_particle(211,cut_p, true, true) + slice->get_num_razzled_particle(13,cut_p, true, true);
            num_razzled_primary_muons = slice->get_num_razzled_particle(13,cut_p, true, true);
            num_razzled_primary_pions = slice->get_num_razzled_particle(211,cut_p, true, true);

            num_chi2_primary_muon_like_candidates = slice->get_num_primary_razzled_muon_like_candidates(cut_p);
            num_razzled_primary_muon_like_candidates = slice->get_num_primary_muon_like_candidates(cut_p);


            E_mu_reco = 0;
            double best_muon_score = 0;
            for(int i_p = 0; i_p < slice->pandora_particle.size(); i_p++) {
                Reco_Particle reco_particle = slice->pandora_particle.at(i_p);
                if(!reco_particle.is_track(cut_p)) continue;
                if(!reco_particle.is_primary(slice->primary_vertex_reco.vertex_cm, cut_p)) continue;
                if((reco_particle.get_razzled_pdg() != 211) && (reco_particle.get_razzled_pdg() != 13)) continue;

                if(reco_particle.razzled_score.muon_score > best_muon_score) {
                    best_muon_score = reco_particle.razzled_score.muon_score;
                    if(reco_particle.track_kinetic_energy != -1) {
                        E_mu_reco = muon_m + reco_particle.track_kinetic_energy;
                    } else {
                        E_mu_reco = 0;
                    }
                }
            }

            output_tree[i_t]->Fill();
        }
    }

    for(int i_t = 0; i_t < num_trees; i_t++) {
        output_file[i_t]->cd();
        output_tree[i_t]->Write();
    }
    for(int i_t = 0; i_t < num_trees; i_t++) {
        output_file[i_t]->Close();
        input[i_t]->Close();
    }
  
} 