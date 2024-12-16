//
// Created by Luis Pelegrina Gutiérrez on 9/10/24.
//

#include "Cluster.h"

double Cluster::get_min_distance_to(int x, double y, int TPC_ID){
  if(this->TPC_ID !=TPC_ID) return 100000000;

  SPoint p(x, y,0,0);
  double min_distance = 10000;
  for(Hit& hit: this->hit_vec) {
    double distance = p.get_distance_to_point_center(hit.plane_point);
    if(distance < min_distance) min_distance = distance;
  }
  return min_distance;
}

double Cluster::get_min_distance_to_cluster(Cluster cluster){
  if(this->TPC_ID != cluster.TPC_ID) return 100000000;

  double min_distance = 10000;
  for(Hit& hit: this->hit_vec) {
    for(Hit hit2: cluster.hit_vec) {
      double distance = hit.plane_point.get_distance_to_point_center(hit2.plane_point);
      if(distance < min_distance) min_distance = distance;
    }
  }
  return min_distance;
}


double Cluster::get_min_distance_to_cluster_w_err(Cluster cluster){
  if(this->TPC_ID != cluster.TPC_ID) return 100000000;

  double min_distance = 10000;
  for(Hit& hit: this->hit_vec) {
    for(Hit hit2: cluster.hit_vec) {
      double distance = hit.plane_point.get_distance_to_point_w_erry(hit2.plane_point);
      if(distance < min_distance) min_distance = distance;
    }
  }
  return min_distance;
}

double Cluster::get_hit_mean_ocuppation() {

  std::map<int, int> occupied_IDs;
  for(Hit &h:this->hit_vec){

    if(occupied_IDs.find(h.Wire_ID)==occupied_IDs.end()){
      occupied_IDs[h.Wire_ID]=1;
    }
    else{
      occupied_IDs[h.Wire_ID]++;
    }
  }

  double mean_occupation = 0;
  for(auto &occ:occupied_IDs){
    mean_occupation+=occ.second;
  }
  mean_occupation /= occupied_IDs.size();

  return mean_occupation;
}

double Cluster::get_mean_hit_width() {
  double mean_width = 0;
  for(Hit hit: this->hit_vec) {
    mean_width += hit.plane_point.get_err_y();
  }
  mean_width /= this->hit_vec.size();

  return mean_width;
}


std::map<int, double> Cluster::calculate_connectedness_map() {
  std::map<int, double> connectedness_map;

  if(this->hit_vec.size() == 1) {
    std::map<int,double> indef_map;
    return indef_map;
  }

  for(Hit hit: this->hit_vec) {
    double min_distance = 1000;
    for(Hit hit2: this->hit_vec) {
      if(hit.ID == hit2.ID) continue;
      double d = hit2.plane_point.get_distance_to_point_w_erry(hit.plane_point);
      if(d < min_distance) min_distance = d;
    }
    connectedness_map[hit.ID] = min_distance;
  }
  return connectedness_map;
}


double projectOntoAxis(const SPoint& point, const TVector2& axis) {
  double dotProduct = point.get_x() * axis.X() + point.get_y()*axis.Y();
  return dotProduct; // Scalar projection
}

LineEquation Cluster::calculate_PCA_line_equation() {
  //cout << endl << endl << endl << "Cluster: " << this->ID << endl;
  // first get the data points
  std::vector<SPoint> data;
  for (Hit& hit : this->hit_vec) {
    data.push_back(hit.plane_point);
  }

  //Get the mean
  double mean_x = 0.0, mean_y = 0.0;
  for(SPoint& p : data) {
    mean_x += p.get_x();
    mean_y += p.get_y();
  }
  mean_x /= data.size();
  mean_y /= data.size();
  //std::cout << "Own method: " << endl;
  //std::cout << "data size: " << data.size() << std::endl;
  //std::cout << "data mean: " << mean_x << " " << mean_y << std::endl;

  //calculate the covariances
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
  //cout << "cov xx: " << cov_xx  << ", cov xy: " << cov_xy  << ", cov xx: " << cov_yy << endl;

  //calculate the 1st eigenvalue
  double eigenvalue_1 = (cov_xx + cov_yy + std::sqrt((cov_xx + cov_yy) * (cov_xx + cov_yy) - 4 * (cov_xx * cov_yy - cov_xy*cov_xy))) / 2.0;

  double eigenvector_x1 = cov_xy/(eigenvalue_1 - cov_xx);
  double eigenvector_y1 = 1;
  double slope_1 = eigenvector_y1/eigenvector_x1;
  double intercept_1 = mean_y - slope_1*mean_x;
  //cout << "eigenvalue 1: " << eigenvalue_1 << " slope 1: " << slope_1 << " intercept 1: " << intercept_1 << endl;

  //calculate the 2nd eigenvalue
  double eigenvalue_2= (cov_xx + cov_yy - std::sqrt((cov_xx + cov_yy) * (cov_xx + cov_yy) - 4 * (cov_xx * cov_yy - cov_xy*cov_xy))) / 2.0;
  double eigenvector_x2 = cov_xy/(eigenvalue_2 - cov_xx);
  double eigenvector_y2 = 1;
  double slope_2 = eigenvector_y2/eigenvector_x2;
  double intercept_2 = mean_y - slope_2*mean_x;
  //cout << "eigenvalue 2: " << eigenvalue_2 << " slope 2: " << slope_2 << " intercept 2: " << intercept_2 << endl;

  //Get the first eigenvalue line equation
  double mae = 0.0;
  for (SPoint& p : data) {
    double predictedY = slope_1 * p.get_x() + intercept_1;
    mae += std::abs(p.get_y() - predictedY);
  }
  mae = mae / data.size();
  double avg_width = 0.0;
  for (SPoint& p : data) {
    avg_width += p.get_err_y();
  }
  avg_width = avg_width / data.size();
  double correlation = mae/avg_width;
  LineEquation first_eigen_value_line_eq(slope_1, intercept_1, correlation);

  //Get the second eigenvalue line equation
  mae = 0.0;
  for (SPoint& p : data) {
    double predictedY = slope_2 * p.get_x() + intercept_2;
    mae += std::abs(p.get_y() - predictedY);
  }
  mae = mae / data.size();
  correlation = mae/avg_width;
  LineEquation second_eigen_value_line_eq(slope_2, intercept_2, correlation);


  std::vector<SPoint> sorted_data = data;
  TVector2 PCA_axis(eigenvector_x1, eigenvector_y1);
  PCA_axis = PCA_axis.Unit();
  cout << PCA_axis.X() << " " << PCA_axis.Y() << endl;

  std::sort(sorted_data.begin(), sorted_data.end(),
            [&PCA_axis](const SPoint& p1, const SPoint& p2) {
              return projectOntoAxis(p1, PCA_axis) < projectOntoAxis(p2, PCA_axis);
            });

  SPoint new_pca_start_point = sorted_data[0];
  SPoint new_pca_end_point = sorted_data[sorted_data.size() - 1];
  this->pca_end_point = new_pca_end_point;
  this->pca_start_point = new_pca_start_point;


  //for(SPoint p: sorted_data) {
  //  cout << p.get_x() << "  " << p.get_y() << " " <<  projectOntoAxis(p, PCA_axis) << endl;
  //}
  //cout << "Start: " << new_pca_start_point.get_x() << " " << new_pca_start_point.get_y() << endl;
  //cout << "End: " << new_pca_end_point.get_x() << " " << new_pca_end_point.get_y() << endl;

/*
  //Using TPrincipal Method
  TPrincipal* pca = new TPrincipal(2, ""); // 3D y la opción "" es importante (ver manual)
  for(SPoint& p : data) {
    Double_t row[2] = {p.get_x(), p.get_y()};
    pca->AddRow(row);
  }
  //Evaluate the PCA
  cout << "TPrincipal method: " << endl;
  pca->MakePrincipals();

  //Get the Eigenvectors and eigenvalues.
  const TMatrixD* Eigenvectors = pca->GetEigenVectors();
  const TVectorD* Eigenvalues  = pca->GetEigenValues();

  double eigenvector_x = (*Eigenvectors)[0][0];
  double eigenvector_y = (*Eigenvectors)[1][0];
  float Eigenvalue = (*Eigenvalues)[0];

  Eigenvectors->Print();
  Eigenvalues->Print();
  double PCA_slope =eigenvector_y/eigenvector_x  ;
  double PCA_intercept = mean_y - eigenvector_y/eigenvector_x*mean_x;
  cout << "PCA slope 1: " << eigenvector_y/eigenvector_x << ", PCA intercept 1: " <<  mean_y - eigenvector_y/eigenvector_x*mean_x << endl;
  eigenvector_x = (*Eigenvectors)[0][1];
  eigenvector_y = (*Eigenvectors)[1][1];
  Eigenvalue = (*Eigenvalues)[1];
  cout << "PCA slope 2: " << eigenvector_y/eigenvector_x << ", PCA intercept 2: " <<  mean_y - eigenvector_y/eigenvector_x*mean_x << endl;

  delete pca;


  cout << first_eigen_value_line_eq.get_slope()<< " " << first_eigen_value_line_eq.get_intercept() << endl;
*/
  //Plotting

  /*
  TCanvas *c1 = new TCanvas();
  TGraphErrors *g = new TGraphErrors();
  for (SPoint& p : data) {
    int n = g->GetN();
    g->SetPoint(n, p.get_x(), p.get_y());
    g->SetPointError(n, p.get_err_x(), p.get_err_y());
  }
  g->SetMarkerColorAlpha(kGray + 2, 0.5);
  g->SetLineColorAlpha(kGray + 2, 0.5);
  g->Draw("ap");


  TGraph* graph = new TGraph();
  // Add first and last points to the graph
  graph->SetPoint(0, new_pca_start_point.get_x(), new_pca_start_point.get_y()); // First point
  graph->SetPoint(1, new_pca_end_point.get_x(), new_pca_end_point.get_y()); // Last point

  // Draw the points
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kRed);

  graph->Draw("p same");


  std::vector<double> x;
  for (SPoint& p : data) {
    x.push_back(p.get_x());
  }
  double x_min = *std::min_element(x.begin(), x.end());
  double x_max = *std::max_element(x.begin(), x.end());

  // Draw a horizontal line from x = 1 to x = 4 at y = 2
  double y1 = first_eigen_value_line_eq.evaluate_x(x_min);
  double y2 = first_eigen_value_line_eq.evaluate_x(x_max);


  TLine *horizontal_line = new TLine(x_min, y1, x_max, y2);
  horizontal_line->SetLineColor(kRed+2); // Set line color
  horizontal_line->SetLineWidth(2);     // Set line width
  horizontal_line->Draw("same");


  // Draw a horizontal line from x = 1 to x = 4 at y = 2
  y1 = second_eigen_value_line_eq.evaluate_x(mean_x -1);
  y2 = second_eigen_value_line_eq.evaluate_x(mean_x + 1);
  cout << x_min << " " << y1 << "  " << x_max << "  " << y2 << endl;

  TLine *horizontal_line_2 = new TLine(x_min, y1, x_max, y2);
  horizontal_line_2->SetLineColor(kBlue+2); // Set line color
  horizontal_line_2->SetLineWidth(2);     // Set line width
  horizontal_line_2->Draw("same");



  c1->Update();
  bool clicked = false;
  while(!clicked) {
    if(c1->WaitPrimitive() == nullptr) {
      std::cout << "Canvas 1 clicked! Proceeding..." << std::endl;
      clicked = true;  // Exit the loop if the correct canvas is clicked
    }
  }


  */

  return first_eigen_value_line_eq;
}