//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_SLICE_H
#define MUPIPROJECT_SLICE_H


#include "TrueInteraction.h"
#include "Cut_Parameters.h"
#include "SpacePoint.h"
#include "RecoVertexWire.h"
#include "Hit.h"
#include "TH1.h"
#include "dQdxSegmentInfo.h"
#include "Vertex.h"
#include "Math/DistFunc.h"
#include "RecoAlgorithms/SPoint.h"

#include <iostream>
#include "RecoParticle.h"

using namespace std;


class Slice {
public:

  Slice()
    : true_interaction(TrueInteraction())
    , file_name("")
    , event_ID(-1)
    , target_months(-1)
    , slice_ID(0)
    , run_ID(-1)
    , subrun_ID(-1)
    , gen_index(-1)
    , weight(-1)
    , truth_matching_purity(-1)
    , truth_matching_completeness(-1)
    , is_reconstructed(false)
    , is_clear_cosmic(false)
    , crumbs_score(-1)
    , nu_score(-1)
    , t0(-1)
    , t0_score(-1)
    , OpT0_fraction(-1)
    , primary_vertex_reco(Vertex())
    , primary_vertex_true(Vertex())
    , vertex_vec({Vertex()})
    , pandora_particle({Reco_Particle()})
    , hits({Hit()})
    , space_point_vec({SpacePoint()})
    {}


    TrueInteraction true_interaction;
    string file_name;

    int target_months;
    int event_ID;
    unsigned int run_ID;// Number assigned to each run
    unsigned int subrun_ID;// Number assigned to each subrun¡¡
    int gen_index; // Generator index, it should be 1 for genie, 2 for corsika and 0 for unknown
    double weight;

    unsigned int slice_ID;
    bool is_reconstructed;
    bool is_clear_cosmic;
    double truth_matching_purity;
    double truth_matching_completeness;
    double crumbs_score;
    double nu_score;

    Vertex primary_vertex_reco;
    Vertex primary_vertex_true;

    double t0;
    double t0_score;
    double OpT0_fraction;

    std::vector<SpacePoint> space_point_vec;
    std::vector<Hit> hits;
    std::vector<Vertex> vertex_vec;
    std::vector<Reco_Particle> pandora_particle;

  int correct_Wire_ID(int input_wire_ID, int Plane_ID, int TPC_ID, int Channel_ID);
  void set_weight(double weight, double expected_POT, double sim_POT);
    void print(Cut_Parameters cut_p, bool b_print_ori_file_name, bool b_print_gen_particles, bool b_print_g4_particles , bool print_vector, bool b_print_hits, bool b_print_reco_vertex);

    void set_corrected_vertex_wire();
    void set_corrected_hits();


    //void set_primary_reco_particle(std::vector<double> primary_track_score_vec, std::vector<double> primary_track_lenght_vec,
    //std::vector<std::vector<double>> primary_razzle_scores, std::vector<std::vector<double>> primary_dazzle_scores, std::vector<std::vector<double>> primary_chi2_scores);

    double get_distance_to_vertex(Hit hit);

    std::vector<Reco_Particle> get_primary_razzled_muon_like_particles(Cut_Parameters cut_p);
    std::vector<Reco_Particle> get_primary_track_particles(Cut_Parameters cut_p);
    Reco_Particle get_longest_track_particle(Cut_Parameters cut_p);
    Reco_Particle get_second_longest_track_particle(Cut_Parameters cut_p);

    double get_num_showers(Cut_Parameters cut_p ,bool check_if_primary, bool check_energy);
    double get_num_tracks(Cut_Parameters cut_p, bool check_if_primary, bool check_quality);


    int get_num_primary_muon_like_chi2_candidates(Cut_Parameters cut_p, bool check_primary);
    int get_num_primary_razzled_muon_like(Cut_Parameters cut_p);

    double get_corrected_time();
    double get_bucket_time();

    int get_num_collection_hits();
    double get_charge_near_vertex(double max_distance);

    bool is_reco_vertex_inside_FV(Cut_Parameters cut_p);

    bool pass_analysis_cut(Cut_Parameters cut_p, string cut_name);

    int get_num_primary_razzled_muon_like_candidates(Cut_Parameters cut_p,  bool check_primary);


    bool pass_containment_cut(string cut_type, double fv_distance, Cut_Parameters cut_p);

    double get_max_razzled_score(int PDG, Cut_Parameters cut_p, bool check_if_primary);
    int get_num_razzled_particle(int PDG, Cut_Parameters cut_p, bool check_if_primary, bool check_if_track, bool check_quality, bool do_proton_rejection);
    TH1* build_dQdx_hist(std::vector<RecoVertexWire> reco_vertex, int target_plane_ID ,double histogram_low_edge, double histogram_high_edge, int id, int num_signal_points);

  dQdxSegmentInfo get_dQdx_segment_info(TH1* h_dQdx, double segment_min_distance, double segment_max_distance);

  TH1 *build_track_dQdx_hist(vector<RecoVertexWire> reco_vertex, int target_plane_ID, double histogram_low_edge,
                             double histogram_high_edge, int id, int num_signal_points, int pfp_ID);

  Reco_Particle get_pandora_p_by_ID(int id);
};

#endif //MUPIPROJECT_SLICE_H
