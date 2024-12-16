//
// Created by Luis Pelegrina Gutiérrez on 9/10/24.
//

#include "ReducedVertexWire.h"
#include "TColor.h"
#include "TCanvas.h"
#include "../Hit.h"
#include "../Vertex.h"
#include "TGraph.h"
#include "Cluster.h"
#include "reco_config.h"


#ifndef CC1PIPROJECT_REFINEDSLICE_H
#define CC1PIPROJECT_REFINEDSLICE_H
class RegionOfInterest {
public:
  vector<int> associated_cluster_ID;
  int ID;

  RegionOfInterest()
  : associated_cluster_ID({-1})
  ,ID(-1)
  {}

  // You can also define a function to display the vector contents
  void print() const {
    std::cout << "Clusters inside the ROI: ";
    for (int value : associated_cluster_ID) {
      std::cout << value << " ";
    }
    std::cout << std::endl;
  }
};

class RefinedSlice {

public:
  vector<int> colors = {kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kMagenta+2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5};
  vector<int> colors_origins = {40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35};

  RefinedSlice()
    : cluster_vec({Cluster()})
    , primary_vertex_true(Vertex())
    , hit_vec({Hit()})
    , ROI_vec({RegionOfInterest()})
    , vertex_vec_C({ReducedVertexWire()})
  {}

  Vertex primary_vertex_true;
  std::vector<Hit> hit_vec;
  std::vector<Cluster> cluster_vec;
  std::vector<RegionOfInterest> ROI_vec;

  std::vector<ReducedVertexWire> vertex_vec_C;

  bool check_hit_conservation();
  bool check_ROI_hit_conservation();

  void remove_isolated_hits(double max_distance_btw_hits, int min_num_neibourghs);
  void add_missing_hits_to_nearest_cluster();
  void add_single_hits_to_closer_cluster();
  void remove_redundant_vertexes(double min_distance_b_ivtx);
  void set_initial_clusters();
  void separate_clusters_by_TPC();
  void add_vertexes_for_isolated_clusters(double min_distance_to_vtx);
  void graph(TCanvas *c, bool plot_associations, bool plot_ROI);
  void set_hits(vector<Hit> hit_vec);
  void set_vtx_vec(vector<Vertex>);
  void print_cluster_info();
  void print_ROI_info();
  void print_vtx_vec();
  void reindex_cluster_vector();
  Cluster merge_clusters(vector<int> cluster_IDs_to_merge, int final_cluster_ID);
  double get_distance_to_closer_cluster(Cluster cluster);

  void remove_isolated_hits_in_clusters(bool verbose);
  void remove_isolated_clusters_with_few_hits();

  void match_hit_vec();

  Hit get_hit_by_ID(int ID);
};


#endif //CC1PIPROJECT_REFINEDSLICE_H
