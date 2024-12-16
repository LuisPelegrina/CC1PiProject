//
// Created by Luis Pelegrina Gutiérrez on 9/10/24.
//
#include "../Slice.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "reco_config.h"


#ifndef CC1PIPROJECT_TPCPLOT_H
#define CC1PIPROJECT_TPCPLOT_H
class TPCPlot {
  public:
    vector<int> colors = {kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed, kGreen-5, kOrange, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5};
    vector<int> colors_origins = {40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35};

  void fill_gr_hits( std::vector<TGraph*> gr_vec,Slice* slice, int TPC, bool set_primary_hits);
  void fill_gr_hits_by_PDG( std::vector<TGraph*> gr_vec, Slice* slice, int TPC, int PDG, bool set_primary_hits);
  void fill_gr_hits_by_ID( std::vector<TGraph*> gr_vec, Slice* slice, int TPC, int ID, bool set_primary_hits);
  void plot_pandora_interaction(Slice* slice, std::vector<TCanvas*> c, bool plot_only_C);
  void fill_gr_vertex(std::vector<TGraph*> gr_vec, vector<RecoVertexWire> vtx);

};


#endif //CC1PIPROJECT_TPCPLOT_H
