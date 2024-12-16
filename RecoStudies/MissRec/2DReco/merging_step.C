#include "../../../Includes.h"
#include "hough_line_step.C"


#include <algorithm>  // For std::sort

TCanvas *c;
bool plot = false;

std::pair<double,double> get_conn_mean_std_dev(Cluster cluster) {
  std::pair<double,double> mean_std_dev_pair;
  map<int, double> connectedness_map_1 = cluster.calculate_connectedness_map();

  double conn_mean = 0;
  for (std:: map<int, double>::iterator it = connectedness_map_1.begin(); it != connectedness_map_1.end(); ++it) {
    conn_mean += it->second;
  }
  conn_mean /= cluster.hit_vec.size();

  //Get connectedness standard deviation
  double conn_std_dev = 0;
  for(Hit hit: cluster.hit_vec) {
    conn_std_dev += (connectedness_map_1[hit.ID] - conn_mean)*(connectedness_map_1[hit.ID] - conn_mean);
  }
  conn_std_dev = sqrt(conn_std_dev/(cluster.hit_vec.size()-1));

  mean_std_dev_pair = {conn_mean, conn_std_dev};

  return mean_std_dev_pair;
}



bool is_neighbour(Cluster main_cluster, Cluster testing_cluster) {
  bool is_neighbour = false;
  /*

  cout << main_cluster.ID << " " << testing_cluster.ID << " "<< main_cluster.get_min_distance_to_cluster_w_err(testing_cluster) << "  " << cos_c1c2 << " " << cos(ang_to_merge_cluster*TMath::Pi()/180)<< endl;
  */

  //Order by PCA axis
  double d1 = main_cluster.pca_start_point.get_distance_to_point_center(testing_cluster.pca_start_point);
  double d2 = main_cluster.pca_start_point.get_distance_to_point_center(testing_cluster.pca_end_point);
  double d3 = main_cluster.pca_end_point.get_distance_to_point_center(testing_cluster.pca_start_point);
  double d4 = main_cluster.pca_end_point.get_distance_to_point_center(testing_cluster.pca_end_point);
  double min_distance = std::min({d1, d2, d3, d4});

  bool is_close = min_distance < distance_to_merge_cluster;

  TVector2 dir_main_cluster_general(1,main_cluster.PCA_line_equation.get_slope());
  dir_main_cluster_general = dir_main_cluster_general.Unit();
  TVector2 dir_testing_cluster_general(1, testing_cluster.PCA_line_equation.get_slope());
  dir_testing_cluster_general = dir_testing_cluster_general.Unit();

  cout << "---" <<  dir_main_cluster_general.X() << " " << dir_main_cluster_general.Y() << " " << dir_testing_cluster_general.X() << " " << dir_testing_cluster_general.Y() << endl;
  double general_cos_c1c2 = dir_main_cluster_general*dir_testing_cluster_general;
  bool is_generaly_connected = false;
  if (general_cos_c1c2 > abs(cos(ang_to_merge_cluster*TMath::Pi()/180))){
    is_generaly_connected = true;
  }


  if(is_close) {
    if(is_generaly_connected) {
      is_neighbour = true;
    } else {
      //few hits connection
      SPoint main_cluster_connection_point(0,0,0,0);
      SPoint testing_cluster_connection_point(0,0,0,0);
      if(min_distance == d1) {
        main_cluster_connection_point = main_cluster.pca_start_point;
        testing_cluster_connection_point = testing_cluster.pca_start_point;
      } else  if(min_distance == d2) {
        main_cluster_connection_point = main_cluster.pca_start_point;
        testing_cluster_connection_point = testing_cluster.pca_end_point;
      } else  if(min_distance == d3) {
        main_cluster_connection_point = main_cluster.pca_end_point;
        testing_cluster_connection_point = testing_cluster.pca_start_point;
      } else  if(min_distance == d4) {
        main_cluster_connection_point = main_cluster.pca_end_point;
        testing_cluster_connection_point = testing_cluster.pca_end_point;
      }

      //get the x points closer to each end and do a PCA
      vector<Hit> sorted_main_cluster_vec = main_cluster.hit_vec;
      vector<Hit> sorted_testing_cluster_cluster_vec = testing_cluster.hit_vec;

      std::sort(sorted_main_cluster_vec.begin(), sorted_main_cluster_vec.end(),
                [&main_cluster_connection_point](const Hit& hit1, const Hit& hit2) {
                  return main_cluster_connection_point.get_distance_to_point_center(hit1.plane_point) < main_cluster_connection_point.get_distance_to_point_center(hit2.plane_point);
                });
      std::sort(sorted_testing_cluster_cluster_vec.begin(), sorted_testing_cluster_cluster_vec.end(),
                [&testing_cluster_connection_point](const Hit& hit1, const Hit& hit2) {
                  return testing_cluster_connection_point.get_distance_to_point_center(hit1.plane_point) < testing_cluster_connection_point.get_distance_to_point_center(hit2.plane_point);
                });


      //num_points_for_PCA_merging
      vector<Hit> closer_main_x_hits;
      if(sorted_main_cluster_vec.size() >= short_num_points_for_PCA_merging) {
        for(int i_h = 0; i_h < short_num_points_for_PCA_merging; i_h++) {
          closer_main_x_hits.push_back(sorted_main_cluster_vec.at(i_h));
        }
      } else {
        for(int i_h = 0; i_h <sorted_main_cluster_vec.size(); i_h++) {
          closer_main_x_hits.push_back(sorted_main_cluster_vec.at(i_h));
        }
      }

      vector<Hit> closer_testing_x_hits;
      if(sorted_testing_cluster_cluster_vec.size() >= short_num_points_for_PCA_merging) {
        for(int i_h = 0; i_h < short_num_points_for_PCA_merging; i_h++) {
          closer_testing_x_hits.push_back(sorted_testing_cluster_cluster_vec.at(i_h));
        }
      } else {
        for(int i_h = 0; i_h <sorted_testing_cluster_cluster_vec.size(); i_h++) {
          closer_testing_x_hits.push_back(sorted_testing_cluster_cluster_vec.at(i_h));
        }
      }

      cout << "main_connection_point" << main_cluster_connection_point.get_x() << " " << main_cluster_connection_point.get_y() << endl;
      for(Hit hit: sorted_main_cluster_vec) {
        cout << hit.plane_point.get_x() << " "<< hit.plane_point.get_y() << endl;
      }
      for(Hit p: closer_main_x_hits) {
        cout << "----"<< p.plane_point.get_x()<< " " << p.plane_point.get_y() << endl;
      }

      cout << "testing_connection_point" << testing_cluster_connection_point.get_x() << " " << testing_cluster_connection_point.get_y() << endl;
      for(Hit hit: sorted_testing_cluster_cluster_vec) {
        cout << hit.plane_point.get_x() << " "<< hit.plane_point.get_y() << endl;
      }
      for(Hit p: closer_testing_x_hits) {
        cout << "----"<< p.plane_point.get_x()<< " " << p.plane_point.get_y() << endl;
      }

      //Perform both PCAs
      Cluster reduced_main_cluster = main_cluster;
      Cluster reduced_testing_cluster = testing_cluster;

      reduced_main_cluster.hit_vec = closer_main_x_hits;
      reduced_testing_cluster.hit_vec = closer_testing_x_hits;

      LineEquation main_line_eq =  reduced_main_cluster.calculate_PCA_line_equation();
      LineEquation testing_line_eq =  reduced_testing_cluster.calculate_PCA_line_equation();


      TVector2 dir_main_cluster(1,main_line_eq.get_slope());
      dir_main_cluster = dir_main_cluster.Unit();
      TVector2 dir_testing_cluster(1, testing_line_eq.get_slope());
      dir_testing_cluster = dir_testing_cluster.Unit();


      bool is_alligned = false;
      double cos_c1c2 = dir_main_cluster*dir_testing_cluster;
      if(cos_c1c2 > abs(cos(ang_to_merge_cluster*TMath::Pi()/180))) is_alligned = true;

      cout << "Cos: " << cos_c1c2 << " " << abs(cos(ang_to_merge_cluster*TMath::Pi()/180)) << endl;
      if(is_alligned) is_neighbour = true;

      /*
      if((main_cluster.pca_start_point.get_y() < 1650./4) &&(main_cluster.pca_start_point.get_x() > 600)) {
        //if(false) {
        if(c==nullptr)
          c = new TCanvas("c", "c",  0, 0, 1000, 800);
        c->cd();

        TGraphErrors *g_hits = new TGraphErrors();
        if(main_cluster.hit_vec.size()>0) {
          for(Hit hit: main_cluster.hit_vec) {
            int n = g_hits->GetN();
            g_hits->SetPoint(n, hit.plane_point.get_x(), hit.plane_point.get_y());
            g_hits->SetPointError(n, 0, hit.plane_point.get_err_y());
          }
        }

        if(testing_cluster.hit_vec.size()>0) {
          for(Hit hit: testing_cluster.hit_vec) {
            int n = g_hits->GetN();
            g_hits->SetPoint(n, hit.plane_point.get_x(), hit.plane_point.get_y());
            g_hits->SetPointError(n, 0, hit.plane_point.get_err_y());
          }
        }
        g_hits->SetLineColorAlpha(kGray + 2, 0.5);
        g_hits->SetMarkerColorAlpha(kGray + 2, 0.5);
        g_hits->Draw("ap");

        if(reduced_main_cluster.hit_vec.size()>0) {
          TGraphErrors *g = new TGraphErrors();
          for(Hit hit: reduced_main_cluster.hit_vec) {
            int n = g->GetN();
            g->SetPoint(n, hit.plane_point.get_x(), hit.plane_point.get_y());
            g->SetPointError(n, 0, hit.plane_point.get_err_y());
          }
          g->SetMarkerStyle(kFullCircle);
          g->SetMarkerColorAlpha(kBlue + 2, 1);
          g->SetLineColorAlpha(kBlue + 2, 1);
          g->Draw("p same");
        }

        if(reduced_testing_cluster.hit_vec.size()>0) {
          TGraphErrors *g = new TGraphErrors();
          for(Hit hit: reduced_testing_cluster.hit_vec) {
            int n = g->GetN();
            g->SetPoint(n, hit.plane_point.get_x(), hit.plane_point.get_y());
            g->SetPointError(n, 0, hit.plane_point.get_err_y());
          }
          g->SetMarkerStyle(kFullCircle);
          g->SetMarkerColorAlpha(kRed + 2, 1);
          g->SetLineColorAlpha(kRed + 2, 1);
          g->Draw("p same");
        }

        std::vector<double> x;
        for (int i_h = 0; i_h < reduced_main_cluster.hit_vec.size(); i_h++) {
          x.push_back(reduced_main_cluster.hit_vec[i_h].Wire_ID);
        }
        double x_min = *std::min_element(x.begin(), x.end());
        double x_max = *std::max_element(x.begin(), x.end());

        // Draw a horizontal line from x = 1 to x = 4 at y = 2
        double y1 = main_line_eq.evaluate_x(x_min);
        double y2 = main_line_eq.evaluate_x(x_max);

        TLine *horizontal_line = new TLine(x_min, y1, x_max, y2);
        horizontal_line->SetLineColor(kBlue); // Set line color
        horizontal_line->SetLineWidth(2);     // Set line width
        horizontal_line->Draw("same");



        std::vector<double> x2;
        for (int i_h = 0; i_h < reduced_testing_cluster.hit_vec.size(); i_h++) {
          x2.push_back(reduced_testing_cluster.hit_vec[i_h].Wire_ID);
        }
        double x2_min = *std::min_element(x2.begin(), x2.end());
        double x2_max = *std::max_element(x2.begin(), x2.end());

        // Draw a horizontal line from x = 1 to x = 4 at y = 2
        y1 = testing_line_eq.evaluate_x(x2_min);
        y2 = testing_line_eq.evaluate_x(x2_max);

        TLine *horizontal_line_2 = new TLine(x2_min, y1, x2_max, y2);
        horizontal_line_2->SetLineColor(kRed); // Set line color
        horizontal_line_2->SetLineWidth(2);     // Set line width
        horizontal_line_2->Draw("same");


        c->Update();
        bool clicked = false;
        while(!clicked) {
          if(c->WaitPrimitive() == nullptr) {
            std::cout << "Canvas 1 clicked! Proceeding..." << std::endl;
            clicked = true;  // Exit the loop if the correct canvas is clicked
          }
        }
      }
       */
    }

  }

  return is_neighbour;
}



std::vector<int> regionQuery( std::vector<Cluster>& cluster_vector,  Cluster& main_cluster) {
  std::vector<int> neighbors;
  for (size_t i = 0; i < cluster_vector.size(); ++i) {
    if(cluster_vector.at(i).ID == main_cluster.ID) continue;
    if(cluster_vector.at(i).TPC_ID != main_cluster.TPC_ID) continue;

    if(is_neighbour(main_cluster, cluster_vector.at(i))) {
      neighbors.push_back(i);
    }

  }
  return neighbors;
}


std::vector<int> expandCluster( std::vector<Cluster>& cluster_vec, int clusterIdx, int merged_clusterIdx, std::vector<int>  cluster_assignment) {

  Cluster cluster = cluster_vec.at(clusterIdx);
  //Gets the neighbours
  std::vector<int> seeds = regionQuery(cluster_vec, cluster);

  //Assings the point that started the cluster the corresponding ID
  cluster_assignment[clusterIdx] = merged_clusterIdx;

  //Do a loop through all the neighbours
  for (size_t i = 0; i < seeds.size(); ++i) {
    //Get the point ID of the seed
    int currentIdx = seeds[i];

    //If seed is not asigned...
    if (cluster_assignment[currentIdx] == 0) {
      //Assing it to the cluster we are expanding
      cluster_assignment[currentIdx] = merged_clusterIdx;

      //If the neighbours are enought expand the current cluster
      std::vector<int> currentSHitNeighbors = regionQuery(cluster_vec, cluster);
      if ((int)currentSHitNeighbors.size() >= 1) {
        cluster_assignment = expandCluster(cluster_vec, currentIdx, merged_clusterIdx, cluster_assignment);
      }
    }
  }
  return cluster_assignment;
}


std::vector<int> fit_cluster_vector( std::vector<Cluster>& cluster_vec) {
  std::vector<int> clusterAssignment;

  clusterAssignment.assign(cluster_vec.size(), 0);  // 0 indicates unassigned

  int merged_cluster_ID = 0;
  for (size_t i = 0; i < cluster_vec.size(); ++i) {

    Cluster cluster = cluster_vec.at(i);

    if (clusterAssignment[i] == 0) {  // Not assigned to any cluster
      std::vector<int> neighbors = regionQuery(cluster_vec, cluster);

      if ((int)neighbors.size() > 0) {
        ++merged_cluster_ID;
        clusterAssignment = expandCluster(cluster_vec, i, merged_cluster_ID, clusterAssignment);
      } else {
        ++merged_cluster_ID;
        clusterAssignment[i] = merged_cluster_ID;
      }
    }
  }

  return clusterAssignment;
}




std::vector<Cluster> perform_high_density_merging(vector<Cluster> cluster_vec) {
  vector<Cluster> new_cluster_vec;
  map<int, int> map_merging_cluster_ID;

  int merging_ID = -1;
  for(Cluster main_cluster: cluster_vec) {
    map_merging_cluster_ID[main_cluster.ID] = merging_ID;
  }

  for(Cluster main_cluster: cluster_vec) {
    if(main_cluster.hit_vec.size() > max_hits_for_high_density_merging) continue;
    double occupation = main_cluster.get_hit_mean_ocuppation();
    if(occupation < min_occupation_for_high_density_merging) continue;

    if( map_merging_cluster_ID[main_cluster.ID] == -1)  {
      map_merging_cluster_ID[main_cluster.ID] = merging_ID;
      merging_ID++;
    }

//Look for candidates for merging nearby
    for(Cluster candidate_cluster: cluster_vec) {
      if(candidate_cluster.hit_vec.size() > max_hits_for_high_density_merging) continue;
      double distance = candidate_cluster.get_min_distance_to_cluster_w_err(main_cluster);
      if(distance < max_distance_for_high_density_merging) map_merging_cluster_ID[candidate_cluster.ID] = merging_ID;
    }
  }

  for(Cluster main_cluster: cluster_vec) {
    int minX = 1e6;
    int maxX = -1e6;
    for (Hit &h: main_cluster.hit_vec) {
      if (h.plane_point.get_x() < minX) minX = h.plane_point.get_x();
      if (h.plane_point.get_x() > maxX) maxX = h.plane_point.get_x();
    }

    cout << "Cluster " << main_cluster.ID << " Merging ID: " << map_merging_cluster_ID[main_cluster.ID] << ", hit density: " << main_cluster.hit_vec.size() * 1.0 / (maxX - minX + 1)
         << ", mean occupation: " << main_cluster.get_hit_mean_ocuppation() << ", num hits: "
         << main_cluster.hit_vec.size() << " Cluster Start :" << main_cluster.pca_start_point.get_x() << " "
         << main_cluster.pca_start_point.get_y() << endl;
  }

  for(Cluster cluster: cluster_vec) {
    if(map_merging_cluster_ID[cluster.ID] == -1) new_cluster_vec.push_back(cluster);
  }

  for(int i_m = 0; i_m<=merging_ID ;i_m++) {
    Cluster new_cluster;
    for(Cluster cluster: cluster_vec) {
      if(map_merging_cluster_ID[cluster.ID] == i_m) {
        new_cluster.hit_vec.insert( new_cluster.hit_vec.end(),  cluster.hit_vec.begin(), cluster.hit_vec.end());
        new_cluster.TPC_ID = cluster.TPC_ID;
        new_cluster.ID = cluster.ID;
      }
    }
    new_cluster.PCA_line_equation = new_cluster.calculate_PCA_line_equation();

    new_cluster_vec.push_back(new_cluster);
  }


  cout << "NEW clusters charac: " << endl;
  for(Cluster cluster: new_cluster_vec) {
    cout << "Cluster " << cluster.ID << " Merging ID: " << map_merging_cluster_ID[cluster.ID]
         << ", mean occupation: " << cluster.get_hit_mean_ocuppation() << ", num hits: "
         << cluster.hit_vec.size() << " Cluster Start :" << cluster.pca_start_point.get_x() << " "
         << cluster.pca_start_point.get_y() << endl;
  }

  return new_cluster_vec;

}




vector<Cluster> perform_linear_cluster_refinement(vector<Cluster> cluster_vec) {
  vector<Cluster> new_cluster_vec;
  for(Cluster cluster: cluster_vec) {
    if((cluster.get_hit_mean_ocuppation() > max_occupation_for_linear_refinement) || (cluster.get_hit_mean_ocuppation() == 1)) {
      new_cluster_vec.push_back(cluster);
      continue;
    }
    //Get the hits with occupatiom > 1
    map<int,int> wire_ocupation;
    for(Hit hit: cluster.hit_vec) {
      wire_ocupation[hit.Wire_ID]++;
    }

    vector<int> degenerated_wires_ID;
    for (std::map<int,int>::const_iterator it = wire_ocupation.begin(); it != wire_ocupation.end(); ++it) {
      std::cout << "Wire: " << it->first << ", Occupatiom: " << it->second << std::endl;
      if(it->second > 1) degenerated_wires_ID.push_back(it->first);
    }

    cout << "Wires with more than 2 occupancy: "<< endl;
    for(int i_hd = 0; i_hd < degenerated_wires_ID.size(); i_hd++) {
      cout << degenerated_wires_ID.at(i_hd) << endl;
    }

    //Create the ROI for the degenerency fixing
    map<int,int> region_hit_ID_map;
    for(Hit hit: cluster.hit_vec) {
      region_hit_ID_map[hit.ID] = -1;
    }

    //Create ROIS of +- 5 wires sourounding the hit, if two ROIS are touhching expand the ROIS
    int current_ROI_ID = -1;
    for(int i_hd = 0; i_hd < degenerated_wires_ID.size(); i_hd++) {
      int current_wire_ID = degenerated_wires_ID.at(i_hd);
      bool first = true;
      int seed_wire_ID = 0;

      for(Hit hit: cluster.hit_vec) {
        if(hit.Wire_ID != current_wire_ID) continue;

        if(region_hit_ID_map[hit.ID] == -1) {
          //Look for the region of the hits in the ends of the ROI we are going to generate
          int origin_wire_ID = -1;
          for(Hit hit2: cluster.hit_vec) {
            if ((hit2.Wire_ID == current_wire_ID - linear_refinement_ROI_parameter) ||
                (hit2.Wire_ID == current_wire_ID + linear_refinement_ROI_parameter)) {
              if (region_hit_ID_map[hit2.ID] != -1) origin_wire_ID = region_hit_ID_map[hit2.ID];
            }
          }

          if(origin_wire_ID != -1) {
            region_hit_ID_map[hit.ID] = origin_wire_ID;
          } else{
            if(first) {
              current_ROI_ID++;
              first = false;
            }
            region_hit_ID_map[hit.ID] = current_ROI_ID;
          }

        }
        seed_wire_ID = region_hit_ID_map[hit.ID];
      }

      for(int i_wid = current_wire_ID - linear_refinement_ROI_parameter; i_wid <= current_wire_ID + linear_refinement_ROI_parameter;i_wid++) {
        for(Hit hit3: cluster.hit_vec) {
          if((hit3.Wire_ID == i_wid) && (region_hit_ID_map[hit3.ID] == -1)) region_hit_ID_map[hit3.ID] = seed_wire_ID;
        }
      }
    }

    cout << "Region hit ID map: "<< endl;
    for(Hit hit: cluster.hit_vec) {
      cout << hit.Wire_ID << " " <<  region_hit_ID_map[hit.ID] << endl;
    }

    cout << "ROI stuff: " << endl;


    map<int,vector<Hit>> ROI_hit_map;
    //Get the wire_IDs with degenerancys in each ROI
    for(int i_ROI = 0; i_ROI <= current_ROI_ID;i_ROI++) {
      cout << "ROI "<< i_ROI << endl;
      map<int,int> ROI_wire_ocupation;
      for (Hit hit: cluster.hit_vec) {
        if (region_hit_ID_map[hit.ID] != i_ROI) continue;
        ROI_wire_ocupation[hit.Wire_ID]++;
      }

      vector<int> ROI_degenerated_wires_ID;
      vector<int> ROI_single_wires_ID;
      for (std::map<int,int>::const_iterator it = ROI_wire_ocupation.begin(); it != ROI_wire_ocupation.end(); ++it) {
        std::cout << "Wire: " << it->first << ", Occupatiom: " << it->second << std::endl;
        if(it->second > 1) ROI_degenerated_wires_ID.push_back(it->first);
        if(it->second == 1) ROI_single_wires_ID.push_back(it->first);
      }

      vector<Hit> base_hit_vec;
      for (Hit hit: cluster.hit_vec) {
        if (std::find(ROI_single_wires_ID.begin(), ROI_single_wires_ID.end(), hit.Wire_ID) != ROI_single_wires_ID.end()) {
          base_hit_vec.push_back(hit);
        }
      }

      cout << endl << "base hit vector: " << endl;
      for(Hit hit: base_hit_vec) {
        cout << "Hit ID: " << hit.ID << " Wire ID: " << hit.Wire_ID << endl;
      }


      vector<vector<Hit>> degenerated_hit_vec;
      for(int i_d = 0;i_d < ROI_degenerated_wires_ID.size();i_d++) {
        int current_wire_ID = ROI_degenerated_wires_ID.at(i_d);
        vector<Hit> current_wire_ID_hit_ID_vec;
        for (Hit hit: cluster.hit_vec) {
          if(hit.Wire_ID == current_wire_ID) current_wire_ID_hit_ID_vec.push_back(hit);
        }
        degenerated_hit_vec.push_back(current_wire_ID_hit_ID_vec);
      }

      cout << endl << "degenerate hit vector: " << endl;
      for(int i = 0; i < degenerated_hit_vec.size(); i++) {
        cout << "Wire: " << ROI_degenerated_wires_ID.at(i) << ", hit IDs:  ";
        for(int j = 0; j < degenerated_hit_vec.at(i).size(); j++) {
            cout << degenerated_hit_vec.at(i).at(j).ID << " ";
        }
        cout << endl;
      }



      vector<string> code = {""};
      vector<string> position_code = {""};
      for(int i = 0; i < degenerated_hit_vec.size(); i++) {
        cout << "ITERATION " << i << endl;
        vector<string> new_code;
        vector<string> new_position_code;
        //vector<Hit> new_hit_vec;
        for(int j = 0; j < degenerated_hit_vec.at(i).size(); j++) {

          for(int k = 0; k < code.size(); k++) {
            new_code.push_back(code.at(k) + to_string(degenerated_hit_vec.at(i).at(j).ID) + "_");
          }

          for(int k = 0; k < position_code.size(); k++) {
            new_position_code.push_back(position_code.at(k) + to_string(j) + "_");
          }
        }

        //hit_vector.push_back(degenerated_hit_vec.at(i).at(j));
        for(int k = 0; k < new_code.size(); k++) {
          cout << new_code.at(k) << endl;
        }
        code = new_code;
        position_code = new_position_code;
      }



      cout << "Final_code: " << endl;
      map<string, vector<Hit>> degenerated_hit_vec_combination;
      for(int k = 0; k < code.size(); k++) {
        cout << code.at(k) << " " <<  position_code.at(k) << endl;
        std::vector<int> positions;
        std::stringstream ss(position_code.at(k));
        std::string segment;
        vector<Hit> hit_vec = base_hit_vec;
        // Read each part of the string separated by '_'
        while (std::getline(ss, segment, '_')) {
          if (!segment.empty()) { // Ignore empty segments (in case the string ends with '_')
            positions.push_back(std::stoi(segment));
          }
        }

        cout << "Positions: ";
        for(int i_p = 0; i_p < positions.size(); i_p++) {
          cout << positions.at(i_p) << " ";
          hit_vec.push_back(degenerated_hit_vec.at(i_p).at(positions.at(i_p)));
        }
        cout << endl;

        degenerated_hit_vec_combination[code.at(k)] = hit_vec;
      }



      for (std::map<string, vector<Hit>>::const_iterator it = degenerated_hit_vec_combination.begin(); it != degenerated_hit_vec_combination.end(); ++it) {
        std::cout << "Combination: " << it->first << endl;
        for(Hit hit: it->second) {
          cout << "Hit ID: " << hit.ID << " Wire ID: " << hit.Wire_ID << " drift_t " << hit.drift_t << endl;
        }
        cout << endl;
      }

      //Get the PCA autovalue for each one
      map<string, double> pca_score_map;
      for (std::map<string, vector<Hit>>::const_iterator it = degenerated_hit_vec_combination.begin(); it != degenerated_hit_vec_combination.end(); ++it) {
        vector<Hit> hit_vec = it->second;

        // first get the data points
        std::vector<SPoint> data;
        for (Hit& hit : hit_vec) {
          data.push_back(hit.plane_point);
        }

        /*
        // Initialize TPrincipal for 2 dimensions
        TPrincipal principal(2, "D");

        // Add each point to TPrincipal
        for (const auto& point : data) {
          Double_t row[2] = {point.get_x(), point.get_y()};
          principal.AddRow(row);
        }

        // Perform PCA analysis
        principal.MakePrincipals();

        // Retrieve eigenvalues and eigenvectors
        const TMatrixD* eigenVectors = principal.GetEigenVectors();
        const TVectorD* eigenValues = principal.GetEigenValues();
        eigenValues->Print();
        eigenVectors->Print();

        cout << "eigenvalue: " << (*eigenValues)[0] << endl;
        pca_score_map[it->first] = (*eigenValues)[0];
        */

        //Get the mean
        double mean_x = 0.0, mean_y = 0.0;
        for(SPoint& p : data) {
          mean_x += p.get_x();
          mean_y += p.get_y();
        }
        mean_x /= data.size();
        mean_y /= data.size();

        //Calculate the covariance matrix
        double cov_xx = 0.0, cov_yy = 0.0, cov_xy = 0.0;
        for(SPoint& p : data) {
          double delta_x = p.get_x() - mean_x;
          double delta_y = p.get_y() - mean_y;
          cov_xx += delta_x * delta_x;
          cov_yy += delta_y * delta_y;
          cov_xy += delta_x * delta_y;
        }
        cov_xx /= (data.size() - 1);
        cov_yy /= (data.size() - 1);
        cov_xy /= (data.size() - 1);

        //calculate the 1st eigenvalue
        double eigenvalue_1 = (cov_xx + cov_yy + std::sqrt((cov_xx + cov_yy) * (cov_xx + cov_yy) - 4. * (cov_xx * cov_yy - cov_xy*cov_xy))) / 2.0;
        //calculate the 2nd eigenvalue
        double eigenvalue_2= (cov_xx + cov_yy - std::sqrt((cov_xx + cov_yy) * (cov_xx + cov_yy) - 4 * (cov_xx * cov_yy - cov_xy*cov_xy))) / 2.0;
       //Normalize both to see which one contributes more
        cout << "normalized eigenvalue: " << endl;
        eigenvalue_1 = eigenvalue_1/(eigenvalue_1+eigenvalue_2);
        eigenvalue_2 = eigenvalue_2/(eigenvalue_1+eigenvalue_2);

        pca_score_map[it->first] = eigenvalue_1;
      }

      cout << "keys with eigenvalue score: " << endl;
      string best_key = "";
      double best_score = 0;
      for (map<string, double>::const_iterator it = pca_score_map.begin(); it != pca_score_map.end(); ++it) {
        cout << "key: " << it->first << " score: " << it->second << endl;
        if(it->second > best_score) {
          best_score = it->second;
          best_key = it->first;
        }
      }
      ROI_hit_map[i_ROI] = degenerated_hit_vec_combination[best_key];


      cout << "ROI " << i_ROI << " end" << endl;
    }


    //Get the non roi hits
    cout << "Region hit ID map: "<< endl;
    vector<Hit> new_hit_vec;
    for(Hit hit: cluster.hit_vec) {
      if( region_hit_ID_map[hit.ID] == -1) new_hit_vec.push_back(hit);
    }
    for(int i_ROI = 0; i_ROI <= current_ROI_ID;i_ROI++) {
      new_hit_vec.insert(new_hit_vec.end(), ROI_hit_map[i_ROI].begin(),ROI_hit_map[i_ROI].end());
    }
    Cluster new_cluster = cluster;
    new_cluster.hit_vec = new_hit_vec;
    new_cluster.calculate_PCA_line_equation();
    new_cluster_vec.push_back(new_cluster);




      if(c==nullptr)
      c = new TCanvas("c", "c",  0, 0, 1000, 800);
    c->cd();

      TGraphErrors *g_hits = new TGraphErrors();
      if(cluster.hit_vec.size()>0) {
        for (Hit hit: cluster.hit_vec) {
          int n = g_hits->GetN();
          g_hits->SetPoint(n, hit.plane_point.get_x(), hit.plane_point.get_y());
          g_hits->SetPointError(n, 0, hit.plane_point.get_err_y());
        }
      }
    g_hits->SetLineColorAlpha(kGray + 2, 0.5);
    g_hits->SetMarkerColorAlpha(kGray + 2, 0.5);
    g_hits->Draw("ap");

    vector<TGraphErrors*> g_ROI_vec;
    for(int i_ROI = 0; i_ROI <= current_ROI_ID;i_ROI++) {
      g_ROI_vec.push_back(new TGraphErrors);
    }

    for(int i_ROI = 0; i_ROI <= current_ROI_ID;i_ROI++) {
      for (Hit hit: cluster.hit_vec) {
        if(region_hit_ID_map[hit.ID] != i_ROI) continue;
        g_ROI_vec.at(i_ROI)->SetPoint(g_ROI_vec.at(i_ROI)->GetN(), hit.plane_point.get_x(), hit.plane_point.get_y());
      }
      g_ROI_vec.at(i_ROI)->SetMarkerColor(colors2.at(i_ROI+1));
    }
    cout << "ROI size: " << g_ROI_vec.size() << endl;

    for(int i_ROI = 0; i_ROI <= current_ROI_ID;i_ROI++) {
      g_ROI_vec.at(i_ROI)->Draw("p same");
    }


    vector<TGraphErrors*> g_ROI_selected_vec;
    for(int i_ROI = 0; i_ROI <= current_ROI_ID;i_ROI++) {
      g_ROI_selected_vec.push_back(new TGraphErrors);
    }

    for(int i_ROI = 0; i_ROI <= current_ROI_ID;i_ROI++) {
      vector<Hit> asd_hit_vec = ROI_hit_map[i_ROI];

      for (Hit hit: asd_hit_vec) {
        g_ROI_selected_vec.at(i_ROI)->SetPoint(g_ROI_selected_vec.at(i_ROI)->GetN(), hit.plane_point.get_x(), hit.plane_point.get_y());
      }
      g_ROI_selected_vec.at(i_ROI)->SetMarkerColor(colors2.at(i_ROI+1));
      g_ROI_selected_vec.at(i_ROI)->SetMarkerStyle(kFullCircle);
    }
    cout << "ROI size: " << g_ROI_selected_vec.size() << endl;

    for(int i_ROI = 0; i_ROI <= current_ROI_ID;i_ROI++) {
      g_ROI_selected_vec.at(i_ROI)->Draw("p same");
    }


    c->Update();
    bool clicked = false;
    while(!clicked) {
      if(c->WaitPrimitive() == nullptr) {
        std::cout << "Canvas 1 clicked! Proceeding..." << std::endl;
        clicked = true;  // Exit the loop if the correct canvas is clicked
      }
    }


  }

  return new_cluster_vec;
}


void merging_step()
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
  GenerateDictionaries();
  TTree *tree;
  TFile *input;

  if(c==nullptr)
    c = new TCanvas("c", "c",  0, 0, 1000, 800);
  c->cd();

  //Declare the variables
  Slice *slice = 0;

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_83k.root";
  input = new TFile(strRuta.c_str());
  tree =(TTree*)input->Get("tree");
  tree->SetBranchAddress("slice", &slice);
  int nEntries = tree->GetEntries();

  int event_to_plot_index = 3449;
  //int event_to_plot_index = 7999;
  //int event_to_plot_index = 123;
  //int event_to_plot_index = 2;
  //int event_to_plot_index = 52;
  //int event_to_plot_index = 4242;

  tree->GetEntry(event_to_plot_index);
  RefinedSlice refined_slice = hough_line_step(slice, c, plot);
  if(verbose_pre_processing) cout << "-+-+-+-+-+-+- Hough Step finished -+-+-+-+-+-+-" << endl;
  refined_slice.print_cluster_info();


  //Check i
  std::vector<int> cluster_assigments = fit_cluster_vector(refined_slice.cluster_vec);
  int max_cluster_ID = 0;
  for(int i = 0; i < cluster_assigments.size();i++) {
    if(max_cluster_ID < cluster_assigments.at(i)) max_cluster_ID = cluster_assigments.at(i);
    cout << "Cluster " << refined_slice.cluster_vec.at(i).ID << " is assigned to: "  << cluster_assigments[i] << endl;
  }
  cout << "MAX CLUSTER ID IS: "<< max_cluster_ID << endl;
  vector<Cluster> new_cluster_vec;
  for(int i_c = 0; i_c < max_cluster_ID; i_c++) {
    Cluster new_cluster;
    new_cluster.ID = i_c;

    for(int i = 0; i < cluster_assigments.size();i++) {
      if(cluster_assigments[i] == i_c+1) {
        Cluster current_cluster = refined_slice.cluster_vec.at(i);
        new_cluster.hit_vec.insert(new_cluster.hit_vec.end(), current_cluster.hit_vec.begin(), current_cluster.hit_vec.end());
        new_cluster.TPC_ID = current_cluster.TPC_ID;
      }
    }
    new_cluster_vec.push_back(new_cluster);
  }

  refined_slice.cluster_vec = new_cluster_vec;
  for(int i_c = 0; i_c < refined_slice.cluster_vec.size();i_c++) {
    refined_slice.cluster_vec.at(i_c).PCA_line_equation =  refined_slice.cluster_vec.at(i_c).calculate_PCA_line_equation();
  }
  if(verbose_merging) refined_slice.print_cluster_info();
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

  //high density, low hits  cluster merging
  refined_slice.cluster_vec = perform_high_density_merging(refined_slice.cluster_vec);
  refined_slice.reindex_cluster_vector();
  if(verbose_merging) refined_slice.print_cluster_info();
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
  cout << "------ HIGH DENSITY LOW HITS MERGING DONE ------- " << endl;



  //Try to make linear clusters better
  //Identify the problematic clusters
  refined_slice.cluster_vec = perform_linear_cluster_refinement( refined_slice.cluster_vec);
  refined_slice.match_hit_vec();
  refined_slice.reindex_cluster_vector();
  if(verbose_merging) refined_slice.print_cluster_info();
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
  cout << "------ HIGH DENSITY LOW HITS MERGING DONE ------- " << endl;







  //Remove clusters away from others and with few hits and clusters with size 1
  /*
  refined_slice.remove_isolated_clusters_with_few_hits();
  refined_slice.reindex_cluster_vector();
  if(verbose_merging) cout << endl << "------- LOW HITS ISOLATED CLUSTER REMOVED ------" << endl;
  if(verbose_merging) refined_slice.print_cluster_info();
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
  */


  //REPEAT THE REMERGING, EXAMPLE CASE 123
  //Recover hits with 1 wire limitation

  //Final graph code
  c->cd();
  refined_slice.graph(c, true, false);

  TCanvas *c2 = new TCanvas("c2");
  c2->cd();
  TGraph *gr_general_hits = new TGraph();

  int TPC = 0;
  if(refined_slice.primary_vertex_true.vertex_cm.X() > 0) TPC = 1;
  for (Hit& hit : slice->hits) {
    if(hit.TPC_ID != TPC) continue;
    if(hit.Plane_ID != 2) continue;
    gr_general_hits->SetPoint(gr_general_hits->GetN(), hit.Wire_ID, hit.drift_t/4);
  }
  gr_general_hits->Draw("AP");
  gr_general_hits->SetMarkerSize(0.5);
  gr_general_hits->SetMarkerColor(kGray+1);
  gr_general_hits->SetMarkerStyle(8);
  for(Cluster cluster: refined_slice.cluster_vec) {
    std::vector<double> x;
    for (Hit& h : cluster.hit_vec) {
      x.push_back(h.plane_point.get_x());
    }
    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());

    // Draw a horizontal line from x = 1 to x = 4 at y = 2
    double y1 =  cluster.PCA_line_equation.evaluate_x(x_min);
    double y2 =  cluster.PCA_line_equation.evaluate_x(x_max);


    TLine *horizontal_line = new TLine(x_min, y1, x_max, y2);
    horizontal_line->SetLineColor(kRed); // Set line color
    horizontal_line->SetLineWidth(2);     // Set line width
    horizontal_line->Draw("same");
  }


  c2->Update();
}


