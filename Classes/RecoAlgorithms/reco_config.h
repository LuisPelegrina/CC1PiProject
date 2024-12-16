//
// Created by Luis Pelegrina Gutiérrez on 9/10/24.
//

#ifndef CC1PIPROJECT_CONFIGRECO_H
#define CC1PIPROJECT_CONFIGRECO_H
//Pre proccessing parameters
bool verbose_pre_processing = true;
bool verbose_interaction = true;
bool verbose_hough = true;
bool verbose_merging = true;

//Vertex association logic
double min_distance_to_consider_cluster_isolated = 40;
double min_distance_between_isolated_vtx = 20;

//Isolated hits removal parameters
double min_neighbours_hits = 2;
double max_distance_to_neighbours = 10;


//Hough Parameters
double max_occupation_for_Hough = 2;
double hough_theta_resolution = 20;
//double max_distance_to_hough_line = 3;
double max_distance_to_hough_line = 5;

double hough_cluster_min_hits = 5;
double min_distance_to_closer_hough_cluster = 4;

double connectedness_width_tolerance = 2;
double sigma_to_consider_in_cluster = 7;
double hit_recover_tolerance = 2;

double distance_to_recover_hit = 5;
double recover_distance_to_PCA = 5;

//DBScan parameters
double min_dbscan_cluster_hits = 3;

//Mergin cluster parameter
double distance_to_merge_cluster = 3;
double ang_to_merge_cluster = 165;
int short_num_points_for_PCA_merging = 4;
int max_hits_for_high_density_merging = 10;
double min_occupation_for_high_density_merging = 1.5;
double max_distance_for_high_density_merging = 3;


double max_occupation_for_linear_refinement = 1.5;
int linear_refinement_ROI_parameter = 3;
#endif