//
// Created by Luis Pelegrina Gutiérrez on 8/5/24.
//
#include "../Includes.h"

double gen_index;
string *nu_type = nullptr;
string input_nu_type;
string *final_state= nullptr;
string input_final_state;
double weight;

double crumbs_score;
bool is_clear_cosmic;
bool is_reconstructed;

double v_x;
double v_y;
double v_z;

double v_x_true;
double v_y_true;
double v_z_true;

double num_primary_tracks;
double num_primary_showers;

double num_razzled_primary_photons;
double num_razzled_primary_electrons;
double num_razzled_primary_protons;
double num_razzled_primary_muons;
double num_razzled_primary_pions;
double num_razzled_primary_muon_like;

double E_0;

double num_razzled_primary_muon_like_candidates;
double num_chi2_primary_muon_like_candidates;

double E_mu_reco;

void Branch(TTree* tree) {
    tree->Branch("gen_index", &gen_index);
    tree->Branch("weight", &weight);
    tree->Branch("E_0", &E_0);

    tree->Branch("nu_type", &input_nu_type);
    tree->Branch("final_state", &input_final_state);
    tree->Branch("crumbs_score", &crumbs_score);

    tree->Branch("is_clear_cosmic", &is_clear_cosmic);
    tree->Branch("is_reconstructed", &is_reconstructed);
    tree->Branch("v_x", &v_x);
    tree->Branch("v_y", &v_y);
    tree->Branch("v_z", &v_z);
    tree->Branch("v_x_true", &v_x_true);
    tree->Branch("v_y_true", &v_y_true);
    tree->Branch("v_z_true", &v_z_true);

    tree->Branch("num_tracks", &num_primary_tracks);
    tree->Branch("num_showers", &num_primary_showers);
    tree->Branch("num_razzled_primary_photons", &num_razzled_primary_photons);
    tree->Branch("num_razzled_primary_electrons", &num_razzled_primary_electrons);
    tree->Branch("num_razzled_primary_protons", &num_razzled_primary_protons);
    tree->Branch("num_razzled_primary_muons", &num_razzled_primary_muons);
    tree->Branch("num_razzled_primary_pions", &num_razzled_primary_pions);
    tree->Branch("num_razzled_primary_muon_like", &num_razzled_primary_muon_like);

    tree->Branch("num_razzled_primary_muon_like_candidates", &num_razzled_primary_muon_like_candidates);
    tree->Branch("num_chi2_primary_muon_like_candidates", &num_chi2_primary_muon_like_candidates);

    tree->Branch("E_mu_reco", &E_mu_reco);
}



void set_Branch(TTree* tree) {
    tree->SetBranchAddress("gen_index", &gen_index);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("E_0", &E_0);

    tree->SetBranchAddress("final_state", &final_state);
    tree->SetBranchAddress("nu_type", &nu_type);
    tree->SetBranchAddress("crumbs_score", &crumbs_score);

    tree->SetBranchAddress("is_clear_cosmic", &is_clear_cosmic);
    tree->SetBranchAddress("is_reconstructed", &is_reconstructed);
    tree->SetBranchAddress("v_x", &v_x);
    tree->SetBranchAddress("v_y", &v_y);
    tree->SetBranchAddress("v_z", &v_z);
    tree->SetBranchAddress("v_x_true", &v_x_true);
    tree->SetBranchAddress("v_y_true", &v_y_true);
    tree->SetBranchAddress("v_z_true", &v_z_true);

    tree->SetBranchAddress("num_tracks", &num_primary_tracks);
    tree->SetBranchAddress("num_showers", &num_primary_showers);
    tree->SetBranchAddress("num_razzled_primary_photons", &num_razzled_primary_photons);
    tree->SetBranchAddress("num_razzled_primary_electrons", &num_razzled_primary_electrons);
    tree->SetBranchAddress("num_razzled_primary_protons", &num_razzled_primary_protons);;
    tree->SetBranchAddress("num_razzled_primary_muons", &num_razzled_primary_muons);
    tree->SetBranchAddress("num_razzled_primary_pions", &num_razzled_primary_pions);
    tree->SetBranchAddress("num_razzled_primary_muon_like", &num_razzled_primary_muon_like);

    tree->SetBranchAddress("num_razzled_primary_muon_like_candidates", &num_razzled_primary_muon_like_candidates);
    tree->SetBranchAddress("num_chi2_primary_muon_like_candidates", &num_chi2_primary_muon_like_candidates);

    tree->SetBranchAddress("E_mu_reco", &E_mu_reco);
}



