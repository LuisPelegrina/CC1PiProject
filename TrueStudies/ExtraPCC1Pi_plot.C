
#include "../Includes.h"
bool save_graphs = true;
bool Normalize = false;


bool is_inside_AV(double v_x, double v_y, double v_z){
    bool is_inside = true;

    if (v_x > 200) is_inside = false;
    if (v_x < -200) is_inside = false;

    if (v_y > 200) is_inside = false;
    if (v_y < -200) is_inside = false;

    if (v_z > 500) is_inside = false;
    if (v_z < 0) is_inside = false;

    return is_inside;
}

void ExtraPCC1Pi_plot()
{
  GenerateDictionaries();

    Cut_Parameters cut_p;

    //Declare the variables
    string strRuta;

    const int num_trees = 1;

    TTree *tree[num_trees];
    TFile *input[num_trees];

    Slice* slice = nullptr;

    for(int i_t = 0; i_t < num_trees; i_t++) {
        if (i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
        input[i_t] = new TFile(strRuta.c_str());
        tree[i_t] = (TTree *) input[i_t]->Get("tree");
        tree[i_t]->SetBranchAddress("slice", &slice);
    }



    map<string, int> interaction_name_map;
    for(int i_t = 0; i_t < num_trees; i_t++) {
        int num_entries = tree[i_t]->GetEntries();

        for (int i_e = 0; i_e < num_entries; ++i_e) {
          tree[i_t]->GetEntry(i_e);
          if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;
          if(slice->slice_ID != 0) continue;
          TrueInteraction true_interaction = slice->true_interaction;
          //if((tm_purity < 0.5)||(tm_completeness < 0.5)) continue;

          double w = 1;
          if (Normalize) w = slice->weight;

          int num_mu = true_interaction.get_generator_num_primary_p(13);
          int num_pi = true_interaction.get_generator_num_primary_p(211);
          if((num_mu != 1) ||( num_pi != 1)) continue;

          interaction_name_map["CC1Pi_old"]+= w;

          int num_primaries = true_interaction.get_generator_num_primaries();
          int num_p = true_interaction.get_generator_num_primary_p(2212);
          int num_n = true_interaction.get_generator_num_primary_p(2112);
          int num_pi0 = true_interaction.get_generator_num_primary_p(111);
          int num_O = num_primaries - num_mu - num_pi;
          if(((num_O - num_n - num_p) == 0)) {
            interaction_name_map["CC1Pi_new"]+= w;
            continue;
          } else {
            interaction_name_map["CC1Pi_extra"]+= w;
            if(num_pi0 > 0) {
              interaction_name_map["CC1Pi_pi0"]+= w;
              continue;
            }
          }


          string name = "-";
          for(GeneratorParticle g_p: true_interaction.gen_primary_particles) {
            if((abs(g_p.PDG) == 211) || (abs(g_p.PDG) == 13) || (abs(g_p.PDG) == 2212) || (abs(g_p.PDG) == 2112)) continue;
            name += "_" + to_string(g_p.PDG);
          }
          interaction_name_map["CC1Pi_extra" + name]+= w;
        }
    }

  for (std::map<string, int>::iterator it = interaction_name_map.begin(); it != interaction_name_map.end(); ++it) {
    std::cout << "Interaction: " << it->first << ", number: " << it->second << std::endl;
  }

    //for(int i_t = 0; i_t < num_trees; i_t++) input[i_t]->Close();
}