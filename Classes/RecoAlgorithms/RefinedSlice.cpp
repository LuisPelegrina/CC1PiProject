//
// Created by Luis Pelegrina Gutiérrez on 9/10/24.
//

#include "RefinedSlice.h"

void RefinedSlice::reindex_cluster_vector() {
  vector<Cluster> dummy_cluster_vector;

  int id = 0;
  for(Cluster& cluster: this->cluster_vec) {
    cluster.ID = id;
    dummy_cluster_vector.push_back(cluster);
    id++;
  }
  this->cluster_vec = dummy_cluster_vector;

}

void RefinedSlice::separate_clusters_by_TPC() {
  //Create two cluster
  vector<Cluster> dummy_cluster_vector;
  for(int i_c = 0; i_c < this->cluster_vec.size();i_c++) {
    Cluster cluster_0;
    Cluster cluster_1;

    cluster_0.ID = i_c;
    cluster_0.TPC_ID = 0;
    cluster_1.ID = i_c + this->cluster_vec.size();
    cluster_1.TPC_ID = 1;

    if(this->cluster_vec.at(i_c).ID == -1) {
      cluster_0.ID = -1;
      cluster_1.ID = -2;
    }

    vector<Hit> dummy_hit_vec_0;
    vector<Hit> dummy_hit_vec_1;
    for(Hit &hit: this->cluster_vec.at(i_c).hit_vec) {
      if(hit.TPC_ID == 0) dummy_hit_vec_0.push_back(hit);
      if(hit.TPC_ID == 1) dummy_hit_vec_1.push_back(hit);
    }
    cluster_0.hit_vec = dummy_hit_vec_0;
    cluster_1.hit_vec = dummy_hit_vec_1;

    if(dummy_hit_vec_0.size() != 0) dummy_cluster_vector.push_back(cluster_0);
    if(dummy_hit_vec_1.size() != 0) dummy_cluster_vector.push_back(cluster_1);
  }

  this->cluster_vec = dummy_cluster_vector;

}

void RefinedSlice::add_missing_hits_to_nearest_cluster() {

  //Look to closer cluster and assing the hit to it
  map<int, vector<Hit>> missing_cluster_hit_vec_map;
  vector<int> processed_hits_ID;

  for(Cluster &unassigned_cluster: this->cluster_vec) {
    if (unassigned_cluster.ID >= 0) continue;

    for (Hit &unassigned_hit: unassigned_cluster.hit_vec) {
      double min_distance = 1000000;
      int closest_cluster_id = -9;

      for (Cluster &cluster: this->cluster_vec) {
        if (cluster.ID  < 0) continue;

        for (Hit &hit: cluster.hit_vec) {
          if(hit.TPC_ID != unassigned_cluster.TPC_ID) continue;
          double distance = hit.plane_point.get_distance_to_point_w_erry(unassigned_hit.plane_point);
          if (distance < min_distance) {
            min_distance = distance;
            closest_cluster_id = cluster.ID;
          }
        }
      }
      if(closest_cluster_id != -9) processed_hits_ID.push_back(unassigned_hit.ID);
      missing_cluster_hit_vec_map[closest_cluster_id].push_back(unassigned_hit);
    }
  }

  vector<Cluster> new_cluster_vector;
  for(Cluster cluster: this->cluster_vec) {
    Cluster new_cluster = cluster;
    if(cluster.ID < 0) {
       vector<Hit> new_hit_vec;
       for(Hit hit: cluster.hit_vec) {
         if (std::find(processed_hits_ID.begin(), processed_hits_ID.end(), hit.ID) == processed_hits_ID.end()) {
           new_hit_vec.push_back(hit);
         }
       }
       new_cluster.hit_vec = new_hit_vec;
    } else {
      new_cluster.hit_vec.insert(new_cluster.hit_vec.end(), missing_cluster_hit_vec_map[cluster.ID].begin(), missing_cluster_hit_vec_map[cluster.ID].end());
    }
    if(new_cluster.hit_vec.size() != 0) new_cluster_vector.push_back(new_cluster);
  }

  this->cluster_vec = new_cluster_vector;
}


void RefinedSlice::remove_redundant_vertexes(double min_distance_b_ivtx) {

  vector<ReducedVertexWire> dummy_vertex_vec;
  for (ReducedVertexWire& vtx: this->vertex_vec_C) {
    bool save = true;
    for (ReducedVertexWire& vtx_2: dummy_vertex_vec) {
      if(vtx.TPC_ID != vtx_2.TPC_ID) continue;
      if(vtx.get_distance_to(vtx_2) < min_distance_b_ivtx) save = false;
    }
    if(save) dummy_vertex_vec.push_back(vtx);
  }
  this->vertex_vec_C = dummy_vertex_vec;

}

void RefinedSlice::graph(TCanvas* c, bool plot_associations, bool plot_ROI) {
  TGraph* gr_general_hits = new TGraph();
  TGraph* gr_primary_vertex_truth = new TGraph();
  TGraph* gr_secondary_vertex = new TGraph();
  vector<TGraph*> gr_clusters_hits;
  int num_graphs = this->cluster_vec.size();
  if(plot_ROI) num_graphs = this->ROI_vec.size();
  for(int i_g = 0; i_g < num_graphs; i_g++) {
    gr_clusters_hits.push_back(new TGraph());
  }

  int TPC = 0;
  if(this->primary_vertex_true.vertex_cm.X() > 0) TPC = 1;

  //bucle para los hits
  for (Hit& hit : this->hit_vec) {
    if(hit.TPC_ID != TPC) continue;
    gr_general_hits->SetPoint(gr_general_hits->GetN(), hit.Wire_ID, hit.drift_t);
  }

  if(plot_ROI) {
    for(int i_g = 0; i_g < num_graphs; i_g++) {
      for(Cluster cluster: cluster_vec) {
        if (std::find(this->ROI_vec.at(i_g).associated_cluster_ID.begin(), this->ROI_vec.at(i_g).associated_cluster_ID.end(), cluster.ID) !=
          this->ROI_vec.at(i_g).associated_cluster_ID.end()) {
          for (Hit &hit: cluster.hit_vec) {
            if (hit.TPC_ID != TPC) continue;
            gr_clusters_hits.at(i_g)->SetPoint(gr_clusters_hits.at(i_g)->GetN(), hit.Wire_ID, hit.drift_t);
          }
          gr_clusters_hits.at(i_g)->SetMarkerSize(0.5);
          gr_clusters_hits.at(i_g)->SetMarkerColor(colors.at(i_g));
          gr_clusters_hits.at(i_g)->SetMarkerStyle(8);
        }
      }
    }

  } else {
    for(int i_g = 0; i_g < num_graphs; i_g++) {
      for (Hit& hit :this->cluster_vec.at(i_g).hit_vec) {
        if(hit.TPC_ID != TPC) continue;
        gr_clusters_hits.at(i_g)->SetPoint(gr_clusters_hits.at(i_g)->GetN(), hit.Wire_ID, hit.drift_t);
      }
      gr_clusters_hits.at(i_g)->SetMarkerSize(0.5);
      gr_clusters_hits.at(i_g)->SetMarkerColor(colors.at(i_g));
      gr_clusters_hits.at(i_g)->SetMarkerStyle(8);
    }
  }
  gr_primary_vertex_truth->SetPoint(gr_primary_vertex_truth->GetN(), this->primary_vertex_true.vertex_wire.at(2).Wire_ID, this->primary_vertex_true.vertex_wire.at(2).drift_t);
  for(int i_v = 0; i_v < this->vertex_vec_C.size(); i_v++) {
    if(this->vertex_vec_C.at(i_v).TPC_ID != TPC) continue;
    gr_secondary_vertex->SetPoint(gr_secondary_vertex->GetN(), this->vertex_vec_C.at(i_v).Wire_ID, this->vertex_vec_C.at(i_v).drift_t);
  }

  for(int i_p = 0; i_p < 3; i_p++) {
    gr_primary_vertex_truth->SetMarkerSize(1);
    gr_primary_vertex_truth->SetMarkerColor(kBlue + 2);
    gr_primary_vertex_truth->SetMarkerStyle(22);

    gr_secondary_vertex->SetMarkerSize(1);
    gr_secondary_vertex->SetMarkerColor(kBlack);
    gr_secondary_vertex->SetMarkerStyle(34);
  }

  gr_general_hits->SetMarkerSize(0.5);
  gr_general_hits->SetMarkerColor(kGray+1);
  gr_general_hits->SetMarkerStyle(8);

  gr_primary_vertex_truth->SetMarkerSize(1);
  gr_primary_vertex_truth->SetMarkerColor(kBlue + 2);
  gr_primary_vertex_truth->SetMarkerStyle(22);

  gr_secondary_vertex->SetMarkerSize(1);
  gr_secondary_vertex->SetMarkerColor(kBlack);
  gr_secondary_vertex->SetMarkerStyle(34);

  TMultiGraph *mg = new TMultiGraph();
  if(gr_general_hits->GetN() > 0) mg->Add(gr_general_hits,"AP");
  if(gr_primary_vertex_truth->GetN() > 0)mg->Add(gr_primary_vertex_truth, "AP");
  if(gr_secondary_vertex->GetN() > 0)mg->Add(gr_secondary_vertex,"AP");
  for(int i_g = 0; i_g < num_graphs; i_g++){
    mg->Add(gr_clusters_hits.at(i_g), "AP");
  }
  mg->SetTitle(";Wire ID; Time [#mu s]");
  mg->Draw("A");

  c->Update();
}


void RefinedSlice::remove_isolated_hits(double max_distance_btw_hits, int min_num_neibourghs) {
  std::vector<Hit> kept_hit_vec;

  for (Hit& hit : this->hit_vec) {
    int num_neighbour_hits = 0;
    for (Hit& hit2 : this->hit_vec) {
      if (hit.ID == hit2.ID) continue;
      if (hit.TPC_ID != hit2.TPC_ID) continue;
      if (hit.Plane_ID != hit2.Plane_ID) continue;

      double distance = hit.plane_point.get_distance_to_point_w_erry(hit2.plane_point);
      if (distance <= max_distance_btw_hits) {
        num_neighbour_hits++;
      }

    }

    if (num_neighbour_hits >= min_neighbours_hits) {
      kept_hit_vec.push_back(hit);
    }
  }

  this->hit_vec = kept_hit_vec;
}


void RefinedSlice::add_vertexes_for_isolated_clusters(double min_distance_fc_to_vtx) {
  for(Cluster& cluster: this->cluster_vec) {
    double min_distance_to_vtx = 1000000;

    for(ReducedVertexWire& vtx: this->vertex_vec_C) {
      if(vtx.TPC_ID != cluster.TPC_ID) continue;
      double distance = cluster.get_min_distance_to(vtx.Wire_ID, vtx.drift_t/4, vtx.TPC_ID);
      if(distance < min_distance_to_vtx) min_distance_to_vtx = distance;
    }

    if(min_distance_to_vtx > min_distance_fc_to_vtx) {
     ReducedVertexWire new_vtx;
     new_vtx.TPC_ID = cluster.TPC_ID;
     new_vtx.drift_t = cluster.hit_vec.at(0).drift_t;
     new_vtx.Wire_ID = cluster.hit_vec.at(0).Wire_ID;
     new_vtx.ID = this->vertex_vec_C.size();
     new_vtx.Plane_ID = 2;
     new_vtx.Channel_ID = cluster.hit_vec.at(0).channel_ID;

     this->vertex_vec_C.push_back(new_vtx);
    }
  }
}

void RefinedSlice::set_hits(vector<Hit> hit_vec) {
  vector<Hit> dummy_hit_vec;
  for (Hit& hit : hit_vec) {
    if(hit.Plane_ID == 2) {
      dummy_hit_vec.push_back(hit);
    }
  }
  this->hit_vec = dummy_hit_vec;

}

void RefinedSlice::set_vtx_vec(vector<Vertex> vertex_vec) {
  vector<ReducedVertexWire> dummy_vertex_vec;
  int i_v = 0;
  for (Vertex& vtx : vertex_vec) {
    ReducedVertexWire dummy_vtx;
    dummy_vtx.drift_t = vtx.vertex_wire.at(2).drift_t;
    dummy_vtx.Wire_ID =  vtx.vertex_wire.at(2).Wire_ID;
    dummy_vtx.Channel_ID =  vtx.vertex_wire.at(2).Channel_ID;
    dummy_vtx.Plane_ID = 2;

    dummy_vtx.TPC_ID = 0;
    if(vtx.vertex_cm.X() > 0) dummy_vtx.TPC_ID = 1;
    dummy_vtx.ID = i_v;
    i_v++;

    dummy_vertex_vec.push_back(dummy_vtx);
  }
  this->vertex_vec_C = dummy_vertex_vec;
}

void RefinedSlice::set_initial_clusters() {
  //Create a Cluster object for each group of pfps, and also 1 for the unassociated hits (-1)
  map<int, int> number_of_hits_in_cluster_map;

  for (Hit& hit : this->hit_vec) {
    number_of_hits_in_cluster_map[hit.associated_pfp_ID]++;
  }

  vector<Cluster> dummy_cluster_vec;
  for (std::map<int, int>::iterator it = number_of_hits_in_cluster_map.begin(); it != number_of_hits_in_cluster_map.end(); ++it) {
    Cluster cluster;
    cluster.ID = it->first;
    dummy_cluster_vec.push_back(cluster);
  }

  //Fill the Clusters
  for(int i_c = 0; i_c < dummy_cluster_vec.size(); i_c++){
    Cluster cluster = dummy_cluster_vec.at(i_c);
    vector<Hit> dummy_cluster_hit_vector;
    for (Hit& hit : this->hit_vec) {
      if(hit.associated_pfp_ID == cluster.ID) dummy_cluster_hit_vector.push_back(hit);
    }
    dummy_cluster_vec.at(i_c).hit_vec = dummy_cluster_hit_vector;
  }
  this->cluster_vec = dummy_cluster_vec;
}

void RefinedSlice::print_cluster_info() {
  int total_hits = 0;
  std::cout << "Number of hits in each cluster:" << std::endl;
  for(Cluster &cluster: this->cluster_vec) {
    cout << "Cluster ID: " << cluster.ID << ", TPC_ID: " << cluster.TPC_ID <<  ", Number of hits: " << cluster.hit_vec.size() << endl;
    total_hits += cluster.hit_vec.size();
  }
  cout << "Cluster total hits: " << total_hits << ", Total hits: " << this->hit_vec.size() << endl;
}
void RefinedSlice::print_vtx_vec() {
  cout << "Vertex of the interaction: " << endl;
  for (ReducedVertexWire& vtx : this->vertex_vec_C) {
    cout << "ID: "<< vtx.ID << " Wire ID: "  << vtx.Wire_ID << " Drift t: " << vtx.drift_t << " Channel ID: " << vtx.Channel_ID << " TPC ID: " << vtx.TPC_ID << endl;
  }
}

bool RefinedSlice::check_hit_conservation() {
  int total_hits = 0;
  bool conserves_hits = true;
  for(Cluster &cluster: this->cluster_vec) {
    total_hits += cluster.hit_vec.size();
  }
  if(total_hits != this->hit_vec.size()) conserves_hits = false;

  return conserves_hits;
}

Cluster RefinedSlice::merge_clusters(vector<int> cluster_IDs_to_merge, int final_cluster_ID) {
  Cluster final_cluster;
  int new_TPC_ID;
  vector<Hit> final_hit_vec;

  bool second = false;
  for(Cluster cluster: this->cluster_vec) {
    if(std::find(cluster_IDs_to_merge.begin(), cluster_IDs_to_merge.end(), cluster.ID) != cluster_IDs_to_merge.end()) {
      if(second && (new_TPC_ID !=  cluster.TPC_ID)) {
        throw std::runtime_error("Clusters to merge are in different TPCs!");
      } else {
        new_TPC_ID = cluster.TPC_ID;
        second = true;
      }
      for(Hit hit: cluster.hit_vec) {
        final_hit_vec.push_back(hit);
      }
    }
  }

  final_cluster.hit_vec = final_hit_vec;
  final_cluster.TPC_ID = new_TPC_ID;
  final_cluster.ID = final_cluster_ID;


  return final_cluster;
}

void RefinedSlice::print_ROI_info() {
  int total_hits = 0;
  for(RegionOfInterest ROI: this->ROI_vec) {
   cout << "ROI ID: " << ROI.ID << endl;
   for(Cluster cluster: this->cluster_vec) {
     if(std::find(ROI.associated_cluster_ID.begin(), ROI.associated_cluster_ID.end(), cluster.ID) != ROI.associated_cluster_ID.end()) {
       cout << "    Cluster ID: " << cluster.ID << ", TPC_ID: " << cluster.TPC_ID <<  ", Number of hits: " << cluster.hit_vec.size() << endl;
       total_hits += cluster.hit_vec.size();
     }
   }
  }
  cout << "ROI total hits: " << total_hits << ", Total hits: " << this->hit_vec.size() << endl;
}


bool RefinedSlice::check_ROI_hit_conservation() {
  int total_hits = 0;
  bool conserves_hits = true;
  for(RegionOfInterest ROI: this->ROI_vec) {
    for (Cluster cluster: this->cluster_vec) {
      if (std::find(ROI.associated_cluster_ID.begin(), ROI.associated_cluster_ID.end(), cluster.ID) !=
      ROI.associated_cluster_ID.end()) {
        total_hits += cluster.hit_vec.size();
      }
    }
  }
  if(total_hits != this->hit_vec.size()) conserves_hits = false;

  return conserves_hits;
}


Hit RefinedSlice::get_hit_by_ID(int ID) {
  for(Hit hit: this->hit_vec) {
    if(hit.ID == ID) return hit;
  }
  return Hit();
}



double RefinedSlice::get_distance_to_closer_cluster(Cluster main_cluster) {
  double min_distance = 100000;
  for(Cluster cluster: this->cluster_vec) {
    if(cluster.ID == main_cluster.ID) continue;
    if(cluster.TPC_ID != main_cluster.TPC_ID) continue;

    for(Hit hit: cluster.hit_vec) {
      for(Hit main_hit: main_cluster.hit_vec) {
        double distance = main_hit.plane_point.get_distance_to_point_w_erry(hit.plane_point);
        if(min_distance > distance) min_distance = distance;
      }
    }
  }
  return min_distance;
}


void RefinedSlice::remove_isolated_hits_in_clusters(bool verbose) {

  vector<Cluster> new_cluster_vec;
  vector<Hit> new_hit_vec;
  //Remove the hits 5 sigmas above average connectedness
  for(Cluster cluster: this->cluster_vec) {

    if(cluster.hit_vec.size() == 1) {
      new_cluster_vec.push_back(cluster);
      new_hit_vec.push_back(cluster.hit_vec.at(0));
      continue;
    }

    map<int, double> connectedness_map = cluster.calculate_connectedness_map();

    //Get mean connectedness
    double conn_mean = 0;
    for (std:: map<int, double>::iterator it = connectedness_map.begin(); it != connectedness_map.end(); ++it) {
      conn_mean += it->second;
    }
    conn_mean /= cluster.hit_vec.size();

    //Get connectedness standard deviation
    double std_dev = 0;
    for(Hit hit: cluster.hit_vec) {
      std_dev += (connectedness_map[hit.ID] - conn_mean)*(connectedness_map[hit.ID] - conn_mean);
    }
    std_dev = sqrt(std_dev/(cluster.hit_vec.size()-1));
    if(verbose) cout << "Cluster " << cluster.ID << ", Mean: " << conn_mean << " std_dev: " << std_dev << endl;

    Cluster new_cluster = cluster;
    vector<Hit> new_cluster_hit_vec;
    for(Hit hit: cluster.hit_vec) {
      if(connectedness_map[hit.ID] <= conn_mean + sigma_to_consider_in_cluster*std_dev) {
        new_cluster_hit_vec.push_back(hit);
        new_hit_vec.push_back(hit);
      }
    }
    new_cluster.hit_vec = new_cluster_hit_vec;

    if(new_cluster_hit_vec.size() != 0) new_cluster_vec.push_back(new_cluster);
  }
  this->cluster_vec = new_cluster_vec;
  this->hit_vec = new_hit_vec;

};

void RefinedSlice::remove_isolated_clusters_with_few_hits() {
  vector<Cluster> new_cluster_vec;
  vector<Hit> new_hit_vec;

  //Loop over all clusters
  for(Cluster cluster: this->cluster_vec) {
    //If it has enought hits save it
    if(cluster.hit_vec.size() > hough_cluster_min_hits) {
      new_cluster_vec.push_back(cluster);
      for(Hit hit: cluster.hit_vec) {
        new_hit_vec.push_back(hit);
      }
    } else if(cluster.hit_vec.size() == 1) {
      //If it has only 1 hit don't save it
      continue;
    } else {
      //If not only save it if it is close to another cluster
      bool  pass = false;
      if(this->get_distance_to_closer_cluster(cluster) <= min_distance_to_closer_hough_cluster) {
        pass = true;
      }
      if(pass) {
        new_cluster_vec.push_back(cluster);
        for(Hit hit: cluster.hit_vec) {
          new_hit_vec.push_back(hit);
        }
      }
    }
  }

  this->cluster_vec = new_cluster_vec;
  this->hit_vec = new_hit_vec;

}

void RefinedSlice::match_hit_vec() {
  vector<Hit> new_hit_vec;
  //Loop over all clusters
  for(Cluster cluster: this->cluster_vec) {
    for(Hit hit: cluster.hit_vec) {
      new_hit_vec.push_back(hit);
    }
  }
  this->hit_vec = new_hit_vec;

}