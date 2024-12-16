#include "../../../Includes.h"

RefinedSlice pre_processing(Slice* slice)
{
  Cut_Parameters cut_p;
/*
  std::vector<TCanvas*> c;
  bool plot_only_C = true;
  if(!plot_only_C) {
    c = {new TCanvas("c1", "My Canvas", 500, 400), new TCanvas("c2", "My Canvas", 500, 400), new TCanvas("c3", "My Canvas", 500, 400)};
    c[0]->SetWindowPosition(10, 50);
    c[1]->SetWindowPosition(700, 50);
    c[2]->SetWindowPosition(10, 500);
  } else {
    c = {new TCanvas("c1", "My Canvas", 800, 600)};
    c[0]->SetWindowPosition(10, 50);
  }


  if(verbose_interaction) slice->print(cut_p, true, true, true, true,false, true);

  TPCPlot tpc_plot;
  tpc_plot.plot_pandora_interaction(slice, c, plot_only_C);
  */

  RefinedSlice refined_slice;
  //Set primary vertex true
  refined_slice.primary_vertex_true = slice->primary_vertex_true;

  //Set Hits
  refined_slice.set_hits(slice->hits);
  if(verbose_pre_processing) cout << endl << "------- Hits created ------" << endl;
  if(verbose_pre_processing) cout << "Hit vector size: " << refined_slice.hit_vec.size() << endl;

  //Remove isolated hits
  refined_slice.remove_isolated_hits(max_distance_to_neighbours, min_neighbours_hits);
  if(verbose_pre_processing) cout << endl << "------- Isolated Hits removed ------" << endl;
  cout << "Hit vector size after removal: " << refined_slice.hit_vec.size() << endl;

  //Create vertex objects
  refined_slice.set_vtx_vec(slice->vertex_vec);
  if(verbose_pre_processing)  {
    cout << endl << "------- Vertexes created ------" << endl;
    refined_slice.print_vtx_vec();
  }
  //Remove extra vertexes
  refined_slice.remove_redundant_vertexes(min_distance_between_isolated_vtx);
  if(verbose_pre_processing)  {
    cout << endl << "------- Redundant vertexes removed  ------" << endl;
    refined_slice.print_vtx_vec();
  }

  //Fill The clusters
  refined_slice.set_initial_clusters();
  if(verbose_pre_processing) cout << endl << "------- Cluster created ------" << endl;
  if(verbose_pre_processing) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false

  //Separate them by TPCs
  refined_slice.separate_clusters_by_TPC();
  if(verbose_pre_processing) cout<< endl  << "------- Cluster Separated by TPC ------" << endl;
  if(verbose_pre_processing) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false

  //Show the output
  //TCanvas* c_reduced = new TCanvas("c_reduced", "My Canvas", 800, 600);
  //c_reduced->cd();
  //c_reduced->SetWindowPosition(900, 50);
  //refined_slice.graph(c_reduced, true, false);

  refined_slice.add_missing_hits_to_nearest_cluster();
  if(verbose_pre_processing) cout << endl << "------- Missing hits added to closest slice- ------" << endl;
  if(verbose_pre_processing) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false

  refined_slice.reindex_cluster_vector();
  if(verbose_pre_processing) cout << endl << "------- Cluster reindexed- ------" << endl;
  if(verbose_pre_processing) refined_slice.print_cluster_info();
  if(!refined_slice.check_hit_conservation()) throw std::runtime_error("Hits in cluster are not equal to total hits!");  // Throw an exception if condition is false

  refined_slice.add_vertexes_for_isolated_clusters(min_distance_to_consider_cluster_isolated);
  if(verbose_pre_processing)  {
    cout << endl << "------- Vertexes added in isolated clusters  ------" << endl;
    refined_slice.print_vtx_vec();
  }

  //TCanvas* c_reduced2 = new TCanvas("c_reduced2", "My Canvas", 800, 600);
  //refined_slice.graph(c_reduced2, true, false);

  return refined_slice;
}