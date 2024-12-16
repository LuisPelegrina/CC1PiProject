#include "../../../Includes.h"
#include "../Graphs_utils.cpp"

void SetPoints2D(TGraph2D* gr2D, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z) {

  gr2D->SetPoint(gr2D->GetN(),max_x+20,max_y+20,max_z+20);
  gr2D->SetPoint(gr2D->GetN(),min_x-20,min_y-20,min_z-20);
}

void PlotEvent()
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
  GenerateDictionaries();
  TTree *tree;
  TFile *input;
    
  Cut_Parameters cut_p;
  cut_p.min_track_lenght = 0;
  cut_p.min_track_score = 0;


  //Declare the variables
  Slice *slice = 0;

  string strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_83k.root";
  input = new TFile(strRuta.c_str());
  tree =(TTree*)input->Get("tree");
  tree->SetBranchAddress("slice", &slice);
  int nEntries = tree->GetEntries();


  int event_to_plot_index = 3449;
  tree->GetEntry(event_to_plot_index);

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

  slice->print(cut_p, true, true, true, true,false, true);

  TPCPlot tpc_plot;
  tpc_plot.plot_pandora_interaction(slice, c, plot_only_C);

}