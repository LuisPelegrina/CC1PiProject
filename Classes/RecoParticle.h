//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_RECOPARTICLE_H
#define MUPIPROJECT_RECOPARTICLE_H

#include "TVector3.h"
#include "RecoVertexWire.h"
#include "Chi2.h"
#include "Razzled.h"
#include "Hit.h"
#include "Cut_Parameters.h"
#include "CaloPoint.h"

#include <iostream>
using namespace std;

class Reco_Particle {
public:

  Reco_Particle()
  :matched_pdg(-1)
  , purity(-1)
  , completeness(-1)
  ,is_pandora_primary(false)
  ,track_score(-1)
  ,track_lenght(-1)
  ,track_kinetic_energy(-1)
  ,track_visible_energy(-1)
  ,track_charge(-1)
  ,true_track_id(-1)
  ,true_energy(-1)
  ,true_P(TVector3(-1,-1,-1))
  ,track_start(TVector3(-1,-1,-1))
  ,start_direction(TVector3(-1,-1,-1))
  ,track_start_wires({RecoVertexWire(), RecoVertexWire(), RecoVertexWire()})
  ,track_end_wires({RecoVertexWire(), RecoVertexWire(), RecoVertexWire()})
  ,track_end(TVector3(-1,-1,-1))
  ,chi2_score(Chi2())
  ,razzled_score(Razzled())
  ,ID(-1)
  ,parent_ID(-1)
  ,has_track(false)
  ,pandora_PDG(-1)
  ,track_theta(-1)
  ,track_phi(-1)
  ,calo_point_vec({CaloPoint()})
  ,shower_direction(TVector3())
  ,shower_length(-1)
  ,shower_opening_angle(-1)
  ,shower_energy(-1)
  ,shower_dEdx(-1)
  {}

    int ID;
    int parent_ID;
    bool is_pandora_primary;
    double track_score;
    bool has_track;
    int pandora_PDG;

    double matched_pdg;
    double true_track_id;
    double true_energy;
    TVector3 true_P;
    double purity;
    double completeness;

    double track_lenght;
    double track_kinetic_energy;
    double track_visible_energy;
    double track_charge;

    TVector3 track_start;
    TVector3 track_end;
    double track_theta;
    double track_phi;
    std::vector<RecoVertexWire> track_start_wires;
     std::vector<RecoVertexWire> track_end_wires;
    TVector3 start_direction;

    TVector3 shower_direction;
    double shower_length;
    double shower_opening_angle;
    double shower_energy;
    double shower_dEdx;


    Chi2 chi2_score;
    Razzled razzled_score;
    std::vector<CaloPoint> calo_point_vec;

    void print(TVector3 vertex_cm);

    bool is_contained_in_fv(Cut_Parameters cut_p);

    bool is_track(Cut_Parameters cut_p);
    bool is_shower(Cut_Parameters cut_p, bool check_energy);

    bool pass_quality_cuts(Cut_Parameters cut_p, vector<Hit> slice_hits);
    bool is_primary(TVector3 vertex_cm, Cut_Parameters cut_p);
    bool is_razzled_muon_like_candidate(Cut_Parameters cut_p);
    bool is_razzled_muon_like(Cut_Parameters cut_p);
    bool is_not_razzled_proton(Cut_Parameters cut_p);
    bool is_muon_like_chi2_candidate(Cut_Parameters cut_p);
};

#endif //MUPIPROJECT_RECOPARTICLE_H
