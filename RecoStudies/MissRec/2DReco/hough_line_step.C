#include "../../../Includes.h"
#include "../Graphs_utils.cpp"
#include "pre_processing.cpp"


#include <algorithm>  // For std::sort

vector<int> colors2 = {kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kRed+2, kGreen+2, kBlue +2 , kViolet+2, kOrange+2, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5};

bool compareByValue(const std::pair<int, int>& a, const std::pair<int, int>& b) {
  return a.second > b.second;  // Descending order
}

bool compareByWireID(const Hit& a, const Hit& b) {
  return a.Wire_ID < b.Wire_ID;
}

double get_hough_score(LineEquation line, std::vector<Hit> hit_list) {
  double weight = 0;
  for (auto hit : hit_list) {
    double distance_to_line = line.get_distance(hit.plane_point);
    weight += pow(distance_to_line,2);
  }

  if(weight>0) {
    return hit_list.size()*hit_list.size()/weight;
  } else{
    return 0;
  }
}

vector<Hit>  get_hits_from_line(std::vector<Hit> hit_list, LineEquation line, double max_distance) {
  std::vector<Hit> line_hit_list;
  for (auto hit : hit_list) {

    double d_to_hypo_line = line.get_distance(hit.plane_point);
    if (d_to_hypo_line < max_distance_to_hough_line) {
      line_hit_list.push_back(hit);
    }
  }
  return line_hit_list;
}

vector<Hit>  get_hough_hits_from_line(std::vector<Hit> hit_list, HoughLine hough_line) {
  return get_hits_from_line(hit_list, hough_line.get_line_equation(), max_distance_to_hough_line);
}


ReducedVertexWire get_vtx_by_ID(vector<ReducedVertexWire> vertex_vec, int ID) {
  ReducedVertexWire found_vtx;
  for(ReducedVertexWire vtx: vertex_vec) {
    if(vtx.ID == ID) found_vtx = vtx;
  }
  return found_vtx;
}


HoughLine get_best_hough_line(std::vector<Hit> hit_list)//, double thetaRes, double maxRadius, double maxDistance, bool displayDirectionHypo = false) {
{
  vector<Hit> pivot_hit_vec = hit_list;
  std::cout<<"  Pivot hits: "<<pivot_hit_vec.size()<<std::endl;

  HoughLine best_Hough_line;

  double theta_step = M_PI*hough_theta_resolution/180;

  for(Hit pivot_hit: pivot_hit_vec) {
    SPoint pivot_p = pivot_hit.plane_point;

    for (double theta = 0; theta < M_PI; theta +=theta_step) {
      LineEquation hypo_line = ComputeLineEquationFromHough(theta, pivot_p);

      std::vector<Hit> hit_Hough_list;
      int num_hits_in_tube = 0;
      for (auto hit : hit_list) {
        if (hit.ID == pivot_hit.ID) continue;

        double d_to_hypo_line = hypo_line.get_distance(hit.plane_point);

        if (d_to_hypo_line < max_distance_to_hough_line) {
          num_hits_in_tube++;
          hit_Hough_list.push_back(hit);
        }
      }
      // get a score to the hough line
      double hypo_line_w = get_hough_score(hypo_line, hit_Hough_list);
      if(num_hits_in_tube==0) hypo_line_w=0;
      if(num_hits_in_tube < 5) hypo_line_w=0;

      if (hypo_line_w > best_Hough_line.get_score()) {
        best_Hough_line.set_line_equation(hypo_line);
        best_Hough_line.set_score(hypo_line_w);
        best_Hough_line.set_nhits(hit_Hough_list.size() + 1); //+1 to include the pivot hit

        cout << best_Hough_line.get_score() << " " << best_Hough_line.get_nhits() <<  " " << best_Hough_line.get_line_equation().get_intercept() << "+" << best_Hough_line.get_line_equation().get_slope() <<  "#cdot x"<< endl;

      }
    }

  }

  /*
  cout << best_Hough_line.get_score() << " " << best_Hough_line.get_nhits() <<  " " << best_Hough_line.get_line_equation().get_intercept() << "+" << best_Hough_line.get_line_equation().get_slope() <<  "#cdot x"<< endl;
  if(c==nullptr)
    c = new TCanvas("c", "c",  0, 0, 1000, 800);
  c->cd();

  if(hit_list.size()>0) {
    TGraphErrors *g = new TGraphErrors();
    for(Hit hit: hit_list) {
      g->SetPoint(g->GetN(), hit.Wire_ID, hit.drift_t*drift_to_wire_conv);
      g->SetPointError(g->GetN(), 0, hit.drift_t_sigma);
    }
    g->SetMarkerColorAlpha(kGray + 2, 0.5);
    g->Draw("ap");
  }

  std::vector<Hit> hit_Hough_list;
  for (auto hit : hit_list) {
    double d_to_best_line= best_Hough_line.get_line_equation().get_distance(SPoint{hit.Wire_ID, hit.drift_t/4});

    if (d_to_best_line < max_distance_to_hough_line) {
      hit_Hough_list.push_back(hit);
    }
  }

  if(hit_Hough_list.size()>0) {
    TGraphErrors *g = new TGraphErrors();
    for(Hit hit: hit_Hough_list) {
      g->SetPoint(g->GetN(), hit.Wire_ID, hit.drift_t*drift_to_wire_conv);
      g->SetPointError(g->GetN(), 0, hit.drift_t_sigma);
    }
    g->SetMarkerColorAlpha(kRed+2, 0.5);
    g->Draw("p same");
  }

  std::vector<double> x;
  for (int i_h = 0; i_h < hit_Hough_list.size(); i_h++) {
    x.push_back(hit_Hough_list[i_h].Wire_ID);
  }
  double x_min = *std::min_element(x.begin(), x.end());
  double x_max = *std::max_element(x.begin(), x.end());

  // Draw a horizontal line from x = 1 to x = 4 at y = 2
  double y1 = best_Hough_line.get_line_equation().evaluate_x(x_min);
  double y2 = best_Hough_line.get_line_equation().evaluate_x(x_max);

  TLine *horizontal_line = new TLine(x_min, y1, x_max, y2);
  horizontal_line->SetLineColor(kRed+2); // Set line color
  horizontal_line->SetLineWidth(2);     // Set line width
  horizontal_line->Draw("same");

  c->Update();
  bool clicked = false;
  while(!clicked) {
    if(c->WaitPrimitive() == nullptr) {
      std::cout << "Canvas 1 clicked! Proceeding..." << std::endl;
      clicked = true;  // Exit the loop if the correct canvas is clicked
    }
  }
*/
  return best_Hough_line;
}

map<int,int> generate_ROI_cluster_association(vector<Cluster> original_cluster_vec, vector<ReducedVertexWire> original_vertex_vec, bool verbose) {

  //Get number of clusters next to each vertex, do it by checking if they are closer to the vertex than a fixed distance
  map<int, int> num_cluster_next_to_vtx;
  for(ReducedVertexWire vtx: original_vertex_vec) {
    for(Cluster cluster: original_cluster_vec) {
      if(cluster.TPC_ID != vtx.TPC_ID) continue;
      if(cluster.get_min_distance_to(vtx.Wire_ID, vtx.drift_t/4, vtx.TPC_ID) < min_distance_to_consider_cluster_isolated) num_cluster_next_to_vtx[vtx.ID]++;
    }
  }
  if(verbose) {
    std::cout << "Contents of the map using iterator:" << std::endl;
    for (std::map<int, int>::iterator it = num_cluster_next_to_vtx.begin(); it != num_cluster_next_to_vtx.end(); ++it) {
      std::cout << "Vtx: " << it->first << ", num cluster nearby: " << it->second << std::endl;
    }
  }

  //Sort the map to get the most populated vertex first
  std::vector<std::pair<int, int>> mapVector(num_cluster_next_to_vtx.begin(), num_cluster_next_to_vtx.end());
  std::sort(mapVector.begin(), mapVector.end(), compareByValue);
  if(verbose) {
    std::cout << "Vertex sorted by number of close clusters in descending order:" << std::endl;
    for (const auto& pair : mapVector) {
      std::cout << "Vtx: " << pair.first << ", num cluster nearby: " << pair.second << std::endl;
    }
  }

  //Create the cluster_ROI_association by looping through  all cluster and checking they are close to the corresponding ROI, if they are clasify them as used and create the association
  map<int, int> cluster_ROI_association;
  vector<int> used_clusters;
  int i_R = 0;

  for (const auto& pair : mapVector) {
    //Get the vertex associated to N number of clusters
    ReducedVertexWire vtx = get_vtx_by_ID(original_vertex_vec, pair.first);
    bool vtx_has_cluster = false;

    //Loop over the clusters
    for(Cluster cluster: original_cluster_vec) {
      if(cluster.TPC_ID != vtx.TPC_ID) continue;

      //Check if the cluster has already been used
      bool used = false;
      for(int i_used = 0; i_used < used_clusters.size() ;i_used++) {
        if(cluster.ID == used_clusters.at(i_used)) used = true;
      }
      if(used) continue;

      //IF the cluster is compatible to be in a ROI associated with the vertex push it into the used clusters and create the association
      if(cluster.get_min_distance_to(vtx.Wire_ID, vtx.drift_t/4, vtx.TPC_ID) < min_distance_to_consider_cluster_isolated) {
        used_clusters.push_back(cluster.ID);
        cluster_ROI_association[cluster.ID] = i_R;
        vtx_has_cluster = true;
      }
    }
    //IF the vtx has any cluster that were not in use increase the int that creates the ROI ID
    if(vtx_has_cluster) i_R++;
  }
  if(verbose) {
    for (std::map<int, int>::iterator it = cluster_ROI_association.begin(); it != cluster_ROI_association.end(); ++it) {
      std::cout << "Cluster: " << it->first << ", Associated ROI: " << it->second << std::endl;
    }
  }
  return cluster_ROI_association;
}


vector<RegionOfInterest> create_ROI_vec(map<int, int> cluster_ROI_association, vector<Cluster> cluster_vec) {

  vector<RegionOfInterest> new_region_of_interest_vec;

  int num_ROIS = 0;
  for (std::map<int, int>::iterator it = cluster_ROI_association.begin(); it != cluster_ROI_association.end(); ++it) {
    if (it->second > num_ROIS) {
      num_ROIS = it->second;
    }
  }
  num_ROIS++;

  for(int i_ROI = 0; i_ROI < num_ROIS; i_ROI++) {
    cout << i_ROI << endl;
    RegionOfInterest new_ROI;
    vector<int> cluster_in_ROI_vec;
    for(Cluster cluster: cluster_vec) {
      if(cluster_ROI_association[cluster.ID] == i_ROI) cluster_in_ROI_vec.push_back(cluster.ID);
    }
    new_ROI.associated_cluster_ID = cluster_in_ROI_vec;
    new_ROI.ID = i_ROI;
    if(cluster_in_ROI_vec.size() != 0)  new_region_of_interest_vec.push_back(new_ROI);
  }

  return new_region_of_interest_vec;
}


vector<Cluster> get_hough_cluster_vector(vector<RegionOfInterest> ROI_vec, vector<Cluster> cluster_vec, bool verbose) {
  vector<Cluster> new_cluster_vec;

  int n_clusters = 0;
  for(RegionOfInterest ROI: ROI_vec) {

    //PICK One cluster and start the hough
    for(Cluster cluster: cluster_vec) {
      if(cluster.get_hit_mean_ocuppation() >= max_occupation_for_Hough) {
        new_cluster_vec.push_back(cluster);
        continue;
      }

      //See if the cluster is inside the ROI
      if (std::find(ROI.associated_cluster_ID.begin(), ROI.associated_cluster_ID.end(), cluster.ID) != ROI.associated_cluster_ID.end()) {
        Cluster cluster_to_hough = cluster;

        //If the occupation per wire is too high skip the hough step
        double mean_occupation = cluster_to_hough.get_hit_mean_ocuppation();
        if(verbose) std::cout<<"Analyzing cluster "<<cluster.ID<<" with "<<cluster.hit_vec.size()<<" hits, MeanOccupation: "<< mean_occupation<<std::endl;
        if(mean_occupation > max_occupation_for_Hough) {
          new_cluster_vec.push_back(cluster);
          continue;
        }

        //Start the hough algorithm
        vector<Hit> hit_vec_to_hough = cluster_to_hough.hit_vec;
        vector<int> used_hits_ID;
        int hough_iter = 0;

        while(abs((int)used_hits_ID.size() - (int)cluster.hit_vec.size()) != 0) {
          if(verbose) std::cout<<" **** Hough iteration: "<< hough_iter<< " # remaining hough hits:"<< hit_vec_to_hough.size()<< std::endl;

          //Get the best Hough line with the available hits and create a Cluster with that
          HoughLine hough_line = get_best_hough_line(hit_vec_to_hough);
          vector<Hit> hough_used_hits = get_hough_hits_from_line(hit_vec_to_hough, hough_line);
          Cluster new_cluster(n_clusters, cluster.TPC_ID, hough_used_hits);
          new_cluster_vec.push_back(new_cluster);

          //Get the ID for the used hits and remove them for the nex iteration hit_vector
          for(Hit hit: hough_used_hits) {
            used_hits_ID.push_back(hit.ID);
          }
          vector<Hit> next_iter_hit_vector;
          for(Hit hit: hit_vec_to_hough) {
            if (std::find(used_hits_ID.begin(), used_hits_ID.end(), hit.ID) == used_hits_ID.end()) {
              next_iter_hit_vector.push_back(hit);
            }
          }
          hit_vec_to_hough = next_iter_hit_vector;

          //Increase the number of clusters and the hough iteration
          hough_iter++;
          n_clusters++;
        }

      } //End of hough algorithm

    } //End of cluster loop

  }//END ROI LOOP

  return new_cluster_vec;
}

vector<Cluster> perform_width_tolerance_declustering(vector<Cluster> cluster_vec, bool verbose){

  vector<Cluster> new_cluster_vec;
  map<int,vector<Hit>>  cluster_hit_vec_associations;
  int current_cluster_id = -1;

  for(Cluster cluster: cluster_vec) {
    vector<Hit> cluster_hits = cluster.hit_vec;

    //get Mean width of hits in cluster
    double mean_width = cluster.get_mean_hit_width();

    // Sort the hits by the X attribute in ascending order
    std::sort(cluster_hits.begin(), cluster_hits.end(), compareByWireID);

    //Prepare the cluster associations for the sorting
    map<int,int>  hits_cluster_associations;
    for(Hit hit: cluster.hit_vec) hits_cluster_associations[hit.ID] = -1;

    //Start the sorting (basic sorting, pick the leftmost hit and change hit id to current_cluster_id, then loop over all the hits and change their ID to the hit_ID if they are connected to it and don't have an ID
    for(Hit hit: cluster_hits) {
      if(hits_cluster_associations[hit.ID] == -1){
        current_cluster_id++;
        hits_cluster_associations[hit.ID] = current_cluster_id;
      }

      for(Hit hit2: cluster_hits) {
        if(hit.ID == hit2.ID) continue;
        if (hit.plane_point.get_distance_to_point_center(hit2.plane_point) < mean_width * connectedness_width_tolerance) {
          if(hit2.ID != -1) hits_cluster_associations[hit2.ID] = hits_cluster_associations[hit.ID];
        }
      }
    }

   //push_back each hit to a map containing in the first the cluster ID and the second the hit verctor associated with that ID
    for(Hit hit: cluster_hits) {
      cluster_hit_vec_associations[hits_cluster_associations[hit.ID]].push_back(hit);
      if(verbose) cout << "ID: " <<   hit.ID << " ,Wire_ID: "<< hit.Wire_ID << " ,drift_t: "<< hit.drift_t << " ,drift_t_sigma: "<< hit.drift_t_sigma<< " ,Point: " << hit.plane_point.get_x() << " " << hit.plane_point.get_y() << " " << hit.plane_point.get_err_y()  << " ,Cluster: " << hits_cluster_associations[hit.ID]  << " " << endl;
    }
  }

  //Fill the new_cluster_vector using the associations
  for (std::map<int, vector<Hit>>::iterator it = cluster_hit_vec_associations.begin(); it != cluster_hit_vec_associations.end(); ++it) {
    Cluster new_cluster(it->first, it->second.at(0).TPC_ID, it->second);
    new_cluster_vec.push_back(new_cluster);
  }
  return new_cluster_vec;
};


vector<Cluster> perform_DBSCAN_clustering(vector<Cluster> cluster_vec, bool verbose) {
  vector<Cluster> new_cluster_vec;
  map<int,vector<Hit>>  cluster_hit_vec_associations;

  int last_cluster_id = -1;
  for(Cluster cluster: cluster_vec) {


    //Get epsilon parameter for dbScan
    double mean_hit_width = cluster.get_mean_hit_width();

    map<int, double> connectedness_map = cluster.calculate_connectedness_map();
    //Get mean connectedness, everything with a connectedness higher than connectedness widht tolerance* widht is skipped
    double conn_mean = 0;
    for (std:: map<int, double>::iterator it = connectedness_map.begin(); it != connectedness_map.end(); ++it) {
      if(it->second > connectedness_width_tolerance * mean_hit_width) continue;
      conn_mean += it->second;
    }
    conn_mean /= cluster.hit_vec.size();
    //Get connectedness standard deviation
    double conn_std_dev = 0;
    for(Hit hit: cluster.hit_vec) {
      if(connectedness_map[hit.ID] > connectedness_width_tolerance * mean_hit_width) continue;
      conn_std_dev += (connectedness_map[hit.ID] - conn_mean)*(connectedness_map[hit.ID] - conn_mean);
    }
    conn_std_dev = sqrt(conn_std_dev/(cluster.hit_vec.size()-1));
    if(verbose) cout << "Cluster " << cluster.ID << ", Mean: " << conn_mean << " std_dev: " << conn_std_dev << endl;

    double final_epsilon = conn_mean + sigma_to_consider_in_cluster *conn_std_dev;

    std::vector<Hit> data = cluster.hit_vec;
    std::sort( data.begin(), data.end(), [&](Hit& h1, Hit& h2) {return h1.plane_point.get_x() < h2.plane_point.get_x();} );

    std::vector<SPoint> data_point_vec;
    for(Hit hit: data) {
      data_point_vec.push_back(hit.plane_point);
    }

    DBSCAN dbscan(final_epsilon, min_dbscan_cluster_hits);
    dbscan.setDistanceFunction(DBSCANHitWidthDistance);
    dbscan.fit(data_point_vec);

    std::vector<int> clusterAssignment = dbscan.getClusterAssignment();
    int max_cluster_id = 0;
    for (size_t i = 0; i < data.size(); ++i) {
      int cluster = clusterAssignment[i] + last_cluster_id;
      if(clusterAssignment[i] != -1) {
        cluster_hit_vec_associations[cluster].push_back(data.at(i));
      }
      if(clusterAssignment[i] > max_cluster_id) max_cluster_id = clusterAssignment[i];
      if (clusterAssignment[i] == -1) {
        std::cout << "SPoint " << data[i].ID << " with position: (" << data[i].plane_point.get_x() << ", " << data[i].plane_point.get_y() << ") is noise" << std::endl;
      } else {
        std::cout << "SPoint " << data[i].ID << " with position: (" << data[i].plane_point.get_x() << ", " << data[i].plane_point.get_y() << ") belongs to cluster " << last_cluster_id + clusterAssignment[i] << std::endl;
      }

    }
    last_cluster_id += max_cluster_id;

  }


  for (std::map<int, vector<Hit>>::iterator it = cluster_hit_vec_associations.begin(); it != cluster_hit_vec_associations.end(); ++it) {
    Cluster new_cluster(it->first, it->second.at(0).TPC_ID, it->second);
    new_cluster_vec.push_back(new_cluster);
  }
  return new_cluster_vec;
}



vector<Cluster> recover_unassigned_hits_using_PCA_line(vector<Cluster> cluster_vec, vector<Hit> unassigned_hit_vec) {
  vector<Cluster> new_cluster_vec;
  vector<int> recovered_hit_vec;

  // Sort the vector of vectors by the size of each inner vector
  /*
  std::sort(cluster_vec.begin(), cluster_vec.end(), [](const Cluster& a, const Cluster& b) {
    return a.hit_vec.size() > b.hit_vec.size();
  });
  */

  for(Cluster cluster: cluster_vec) {
   new_cluster_vec.push_back(cluster);
  }

  for(Hit missed_hit: unassigned_hit_vec) {

    /*
    Cluster best_PCA_cluster;

    double closer_PCA_distance = 1000;
    for(Cluster cluster: cluster_vec) {
      if(cluster.get_min_distance_to(missed_hit.plane_point.get_x(), missed_hit.plane_point.get_y()) > distance_to_recover_hit) continue;
      double distance_to_PCA = cluster.PCA_line_equation.get_distance(missed_hit.plane_point);
      if(distance_to_PCA < closer_PCA_distance)  {
        closer_PCA_distance = distance_to_PCA;
        best_PCA_cluster = cluster;
      }
    }
     */
    vector<Cluster> compatible_cluster_vec;
    for(Cluster cluster: cluster_vec) {
      if (cluster.get_min_distance_to(missed_hit.plane_point.get_x(), missed_hit.plane_point.get_y(), missed_hit.TPC_ID) > distance_to_recover_hit) continue;

      double distance_to_PCA = cluster.PCA_line_equation.get_distance(missed_hit.plane_point);
      if(distance_to_PCA < connectedness_width_tolerance*missed_hit.plane_point.get_err_y()) compatible_cluster_vec.push_back(cluster);
    }

    Cluster best_PCA_cluster;
    double max_size = 0;
   for(Cluster compatible_cluster: compatible_cluster_vec) {
     if(compatible_cluster.hit_vec.size() > max_size) {
       max_size = compatible_cluster.hit_vec.size();
       best_PCA_cluster = compatible_cluster;
     }
   }

    for(int i_c= 0; i_c < new_cluster_vec.size();i_c++) {
      if(new_cluster_vec.at(i_c).ID == best_PCA_cluster.ID) new_cluster_vec.at(i_c).hit_vec.push_back(missed_hit);
    }

    /*
    if(closer_PCA_distance < connectedness_width_tolerance*missed_hit.plane_point.get_err_y()) {
      for(int i_c= 0; i_c < new_cluster_vec.size();i_c++) {
        if(new_cluster_vec.at(i_c).ID == best_PCA_cluster.ID) new_cluster_vec.at(i_c).hit_vec.push_back(missed_hit);
      }
    }
     */

  }

  for(int i_c= 0; i_c < new_cluster_vec.size();i_c++) {
    new_cluster_vec.at(i_c).PCA_line_equation =  new_cluster_vec.at(i_c).calculate_PCA_line_equation();
  }


  /*
  for(Cluster cluster: cluster_vec) {
    cout << "Cluster: " << cluster.ID << endl;
    Cluster new_cluster = cluster;
    vector<Hit> new_cluster_hit_vec = cluster.hit_vec;

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

    for(Hit missed_hit: unassigned_hit_vec) {

      //If the hit has already been assigned continue
      if(std::find(recovered_hit_vec.begin(), recovered_hit_vec.end(), missed_hit.ID) != recovered_hit_vec.end()) continue;

      //Calculate how close the hit is to the PCA, if it is to far go to the next hit
      double distance_to_PCA = cluster.PCA_line_equation.get_distance(missed_hit.plane_point);
      if(distance_to_PCA > recover_distance_to_PCA) continue;

      //Calculate the connectedness off the hit to the cluster we are trying to complete
      double hit_conn = 1000;
      for(Hit hit: cluster.hit_vec) {
        double distance_to_hit = hit.plane_point.get_distance_to_point_w_erry(missed_hit.plane_point);
        if( distance_to_hit < hit_conn ) hit_conn = distance_to_hit;
      }

      //If the hit connectedness is compatible with the cluster and the widht of the hit is compatible with the
      //hit distance to the line recover it

      TVector2 dir_cluster(1,cluster.PCA_line_equation.get_slope());
      dir_cluster = dir_cluster.Unit();

      bool is_close_to_PCA = cluster.PCA_line_equation.get_distance(missed_hit.plane_point) < connectedness_width_tolerance*missed_hit.plane_point.get_err_y();
      if((hit_conn <= hit_recover_tolerance*cluster.get_mean_hit_width()) && is_close_to_PCA) {
        recovered_hit_vec.push_back(missed_hit.ID);
        new_cluster_hit_vec.push_back(missed_hit);
        cout << hit_conn << " " << distance_to_PCA << "  "<< missed_hit.ID << endl;
      }
    } //End of loop looking in missing hits

    new_cluster.hit_vec = new_cluster_hit_vec;

    //Refresh the PCA line equation
    new_cluster.PCA_line_equation =  new_cluster.calculate_PCA_line_equation();
    new_cluster_vec.push_back(new_cluster);
  }
   */

  return new_cluster_vec;
}



RefinedSlice hough_line_step(Slice* slice, TCanvas *c, bool plot)
{
  RefinedSlice refined_slice = pre_processing(slice);
  if(verbose_pre_processing) cout << "-+-+-+-+-+-+- Pre processing finished -+-+-+-+-+-+-" << endl;
  refined_slice.print_vtx_vec();
  refined_slice.print_cluster_info();

  //Create Region of interest for each vertex (to start the algorithm per ROI)
  map<int, int> cluster_ROI_association = generate_ROI_cluster_association(refined_slice.cluster_vec, refined_slice.vertex_vec_C, true);
  vector<RegionOfInterest> ROI_vec = create_ROI_vec(cluster_ROI_association, refined_slice.cluster_vec);
  refined_slice.ROI_vec = ROI_vec;
  if(verbose_hough) cout << endl << "------- ROI Created -------" << endl;
  if(verbose_hough) refined_slice.print_ROI_info();
  if(!refined_slice.check_ROI_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false
  if(plot) {
    refined_slice.graph(c, true, true);
    c->Update();
    bool clicked = false;
    while(!clicked) {
      if(c->WaitPrimitive() == nullptr) {
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }
  }

  //Start the hough line algorithm (Create hough lines for each ROI)
  refined_slice.cluster_vec = get_hough_cluster_vector(refined_slice.ROI_vec, refined_slice.cluster_vec, verbose_hough);
  if(verbose_hough) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false
  if(plot) {
    refined_slice.graph(c, true, false);
    c->Update();
    bool clicked = false;
    while(!clicked) {
      if(c->WaitPrimitive() == nullptr) {
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }
  }

  //Do some declustering using DBSCAN, separate each hough cluster by grouping all hit closer to each other
  refined_slice.cluster_vec = perform_DBSCAN_clustering(refined_slice.cluster_vec, verbose_hough);
  refined_slice.match_hit_vec();
  refined_slice.reindex_cluster_vector();
  if(verbose_hough) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false
  if(plot) {
    refined_slice.graph(c, true, false);
    c->Update();
    bool clicked = false;
    while(!clicked) {
      if(c->WaitPrimitive() == nullptr) {
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }
  }

  //Remove all hits 5 sigmas above average connectedness from the clusters
  refined_slice.remove_isolated_hits_in_clusters(true);
  refined_slice.reindex_cluster_vector();
  if(verbose_hough) cout << endl << "------- Cluster hits more than avrg+5*std away from other hits removed -------" << endl;
  if(verbose_hough) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false
  if(plot) {
    refined_slice.graph(c, true, false);
    c->Update();
    bool clicked = false;
    while(!clicked) {
      if(c->WaitPrimitive() == nullptr) {
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }
  }


  //Remove clusters away from others and with few hits and clusters with size 1
  refined_slice.remove_isolated_clusters_with_few_hits();
  refined_slice.reindex_cluster_vector();
  if(verbose_hough) cout << endl << "------- Cluster reindexed- ------" << endl;
  if(verbose_hough) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false
  if(plot) {
    refined_slice.graph(c, true, false);
    c->Update();
    bool clicked = false;
    while(!clicked) {
      if(c->WaitPrimitive() == nullptr) {
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }
  }



  //Calculate PCA Line equation for all clusters
  for(int i_c = 0; i_c < refined_slice.cluster_vec.size();i_c++) {
    refined_slice.cluster_vec.at(i_c).PCA_line_equation =  refined_slice.cluster_vec.at(i_c).calculate_PCA_line_equation();
  }

  //Try to recover hits that are unasigned
  //Current_hits_IDS
  vector<int> used_hit_ID;
  for(Hit hit: refined_slice.hit_vec) used_hit_ID.push_back(hit.ID);
  //get Unassigned hits
  vector<Hit> unassigned_hit_vec;
  for(Hit hit: slice->hits) {
    if(hit.Plane_ID != 2) continue;
    if (std::find(used_hit_ID.begin(), used_hit_ID.end(), hit.ID) == used_hit_ID.end()) {
      unassigned_hit_vec.push_back(hit);
    }
  }
  refined_slice.cluster_vec = recover_unassigned_hits_using_PCA_line(refined_slice.cluster_vec, unassigned_hit_vec);
  refined_slice.match_hit_vec();
  refined_slice.reindex_cluster_vector();
    if(verbose_hough) cout << endl << "------- Missed hits recovered -------" << endl;
  if(verbose_hough) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false
  if(plot) {
    refined_slice.graph(c, true, false);
    c->Update();
    bool clicked = false;
    while(!clicked) {
      if(c->WaitPrimitive() == nullptr) {
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }
  }

  return refined_slice;
}


