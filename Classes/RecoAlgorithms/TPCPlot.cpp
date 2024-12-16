//
// Created by Luis Pelegrina Gutiérrez on 9/10/24.
//

#include "TPCPlot.h"

void TPCPlot::fill_gr_hits( std::vector<TGraph*> gr_vec,Slice* slice, int TPC, bool set_primary_hits) {
  //bucle para los hits
  for(Hit& hit: slice->hits) {
    if (hit.TPC_ID != TPC) continue;
    Reco_Particle associated_pandora_p = slice->get_pandora_p_by_ID(hit.associated_pfp_ID);
    if(set_primary_hits && !associated_pandora_p.is_pandora_primary) continue;

    gr_vec[hit.Plane_ID]->SetPoint(gr_vec[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
  }

  for(int i_plane = 0; i_plane < 3; i_plane++) {
    gr_vec[i_plane]->SetTitle(";Wire ID; Time [#mu s]");
    gr_vec[i_plane]->SetMarkerSize(0.5);
    gr_vec[i_plane]->SetMarkerColor(kGray+2);
    gr_vec[i_plane]->SetMarkerStyle(8);
  }
}

void TPCPlot::fill_gr_hits_by_PDG( std::vector<TGraph*> gr_vec, Slice* slice, int TPC, int PDG, bool set_primary_hits) {

  //bucle para los hits
  for(Hit& hit: slice->hits) {
    if (hit.TPC_ID != TPC) continue;
    Reco_Particle associated_pandora_p = slice->get_pandora_p_by_ID(hit.associated_pfp_ID);
    if(set_primary_hits && !associated_pandora_p.is_pandora_primary) continue;

    if(associated_pandora_p.matched_pdg != PDG) continue;

    gr_vec[hit.Plane_ID]->SetPoint(gr_vec[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
  }

  for(int i_plane = 0; i_plane < 3; i_plane++) {
    gr_vec[i_plane]->SetMarkerSize(0.5);
    gr_vec[i_plane]->SetMarkerColor(PDG);
    gr_vec[i_plane]->SetMarkerStyle(8);
  }
}

void TPCPlot::fill_gr_hits_by_ID( std::vector<TGraph*> gr_vec, Slice* slice, int TPC, int ID, bool set_primary_hits) {
  //bucle para los hits

  for(Hit& hit: slice->hits) {
    if(hit.TPC_ID != TPC) continue;
    Reco_Particle associated_pandora_p = slice->get_pandora_p_by_ID(hit.associated_pfp_ID);
    if(set_primary_hits && !associated_pandora_p.is_pandora_primary) continue;

    if(hit.associated_pfp_ID != ID) continue;

    gr_vec[hit.Plane_ID]->SetPoint(gr_vec[hit.Plane_ID]->GetN(), hit.Wire_ID, hit.drift_t);
  }

  for(int i_plane = 0; i_plane < 3; i_plane++) {
    gr_vec[i_plane]->SetMarkerSize(0.5);
    gr_vec[i_plane]->SetMarkerColor(colors.at(ID));
    gr_vec[i_plane]->SetMarkerStyle(8);
  }
}

void TPCPlot::fill_gr_vertex(std::vector<TGraph*> gr_vec, vector<RecoVertexWire> vtx) {
  for(int i_p = 0; i_p < 3; i_p++) {
    double vertex_drift_t = vtx.at(i_p).drift_t;
    int vertex_wire_ID = vtx.at(i_p).Wire_ID;
    gr_vec[i_p]->SetPoint(gr_vec[i_p]->GetN(), vertex_wire_ID, vertex_drift_t);
  }
}


void TPCPlot::plot_pandora_interaction(Slice* slice, std::vector<TCanvas*> c, bool plot_only_C) {
  std::vector<TGraph*> gr_general_hits = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<std::vector<TGraph*>> gr_pfp_hits;
  std::vector<std::vector<TGraph*>> gr_primary_pfp_hits;

  std::vector<TGraph*> gr_primary_vertex = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_primary_vertex_truth = {new TGraph(), new TGraph(), new TGraph()};
  std::vector<TGraph*> gr_secondary_vertex = {new TGraph(), new TGraph(), new TGraph()};

  int TPC = 0;
  if(slice->primary_vertex_true.vertex_cm.X() > 0) TPC = 1;
  cout << "TPC Chosen" << endl;

  fill_gr_hits(gr_general_hits, slice, TPC, false);


  for(int i_p = 0; i_p < slice->pandora_particle.size();i_p++) {
    gr_pfp_hits.push_back({new TGraph(), new TGraph(), new TGraph()});
    gr_primary_pfp_hits.push_back({new TGraph(), new TGraph(), new TGraph()});
    fill_gr_hits_by_ID(gr_pfp_hits.at(i_p), slice, TPC, slice->pandora_particle.at(i_p).ID, false);
    fill_gr_hits_by_ID(gr_primary_pfp_hits.at(i_p), slice, TPC, slice->pandora_particle.at(i_p).ID, true);
  }
  cout << "hits done" << endl;
  fill_gr_vertex(gr_primary_vertex, slice->primary_vertex_reco.vertex_wire);
  fill_gr_vertex(gr_primary_vertex_truth, slice->primary_vertex_true.vertex_wire);

  for(Vertex& vtx: slice->vertex_vec) {
    fill_gr_vertex(gr_primary_vertex, vtx.vertex_wire);
  }

  for(int i_p = 0; i_p < 3; i_p++) {
    gr_primary_vertex[i_p]->SetMarkerSize(1);
    gr_primary_vertex[i_p]->SetMarkerColor(kBlack);
    gr_primary_vertex[i_p]->SetMarkerStyle(8);

    gr_primary_vertex_truth[i_p]->SetMarkerSize(1);
    gr_primary_vertex_truth[i_p]->SetMarkerColor(kBlue + 2);
    gr_primary_vertex_truth[i_p]->SetMarkerStyle(22);

    gr_secondary_vertex[i_p]->SetMarkerSize(1);
    gr_secondary_vertex[i_p]->SetMarkerColor(kBlack);
    gr_secondary_vertex[i_p]->SetMarkerStyle(34);
  }

  //Dibuja las 3 vistas
  for(int i_p = 0; i_p < 3; i_p++) {
    if((plot_only_C) && (i_p != 2)) continue;
    cout << i_p << endl;
    if(plot_only_C) {
      c[0]->cd();
    } else {
      c[i_p]->cd();
    }
    cout << "Start plotting" << endl;

    gr_secondary_vertex[i_p]->SetMarkerStyle(34);
    TMultiGraph *mg = new TMultiGraph();
    if(gr_general_hits[i_p]->GetN() > 0) mg->Add(gr_general_hits[i_p],"AP");
    for(int i_part = 0; i_part < slice->pandora_particle.size();i_part++) if(gr_primary_pfp_hits[i_part][i_p]->GetN() > 0) mg->Add(gr_primary_pfp_hits[i_part][i_p],"AP");
    if(gr_primary_vertex_truth[i_p]->GetN() > 0)mg->Add(gr_primary_vertex_truth[i_p],"AP");
    if(gr_primary_vertex[i_p]->GetN() > 0)mg->Add(gr_primary_vertex[i_p],"AP");
    if(gr_secondary_vertex[i_p]->GetN() > 0)mg->Add(gr_secondary_vertex[i_p],"AP");

    mg->SetTitle(";Wire ID; Time [#mu s]");
    mg->Draw("A");
    if(plot_only_C) {
      c[0]->Update();
    } else {
      c[i_p]->Update();
    }

  }

}