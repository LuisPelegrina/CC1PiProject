#include "tree_utils.cpp"
#include "../Includes.h"

void read_tree()
{
    //Definition of the file from which the tree is readed and the folder inside the ".root" file where the tree is stored
    TFile *input_file;
    TDirectoryFile *tree_dir;

    //Definition of the tree that stores information about each LArSoft event
    TTree *event_tree;

    //Definition of the tree that stores information about each LArSoft subrun (number of events generated in the subrun, POT...)
    TTree *subrun_tree;

    //Definition of the path where the file containing the truth tree is stored, along with the code that opens the file and both the subrun and the event by event tree
    string file_name = "83k";
    string path_to_tree = "../Data/analysis_output_" + file_name + ".root";
    input_file = new TFile(path_to_tree.c_str());
    tree_dir = (TDirectoryFile*)input_file->Get("ana");
    event_tree =(TTree*)tree_dir->Get("tree");
    subrun_tree =(TTree*)tree_dir->Get("subrun_tree");

    //Code that sets the Branches fot the subrun and event trees in order to use the information they store, please visit "tree_utils.cpp" to see how the functions are defined
    set_branch(event_tree);
    set_branch_subtree(subrun_tree);

    //Example of how the subrun_tree is used. In this case the number of entries in the tree is readed. After that a loop calculates the total POT simulated to create the root tree and prints it on the screen.
    //Declare a counter for the total POT
    double total_pot = 0;
    //Calculate the number of entries of the subrun_tree
    int num_subrun_entries = subrun_tree->GetEntries();
    //Print the number of subruns needed to simulate the truth_tree
    cout << "Total number of subruns: " << num_subrun_entries << endl;
    //Start the loop on each subrun (each subrun_tree entry)
    for(int i_e = 0; i_e < num_subrun_entries; i_e++) {
        //Start reading the entry with index i_e
        subrun_tree->GetEntry(i_e);
        //Add the POT of that entry to the total POT
        total_pot += POT;
        //Print the subrun index and the POT simulated in that Subrun
        cout << "The subrun " << i_e << " simulated a total of " << POT  << " POT"<< endl;
    }
    //Print the total POT simulated in the whole tree
    cout << "Total POT: " << total_pot << endl;

    //Example of how the event_tree is used. In this case the number of entries in the tree is readed. After that a loop prints on screen all the information in an human-readable way. Please refer to "tree_utils.cpp" for a definition of each variable printed.
    //Calculate the number of entries of the event_tree
    int num_entries = event_tree->GetEntries();
    //Print the number of events on the truth_tree
    cout << "Total number of events: " << num_entries << endl;
    //Start the loop on each event (each event_tree entry)
    for(int i_e = 0; i_e < num_entries; i_e++) {
        //Start reading the entry with index i_e
        event_tree->GetEntry(i_e);

        //Print relevant information about the tree
        cout << "---- run " << run_ID << " ----" << endl;
        cout << "---- run " << subrun_ID << " ----" << endl;
        cout << "---- event " << event_ID << " ----" << endl;

        //Print relevant informationa about the neutrino that generated the interaction
        cout << "Neutrino properties: " << endl;
        cout << "PDG: " << nu_PDG << " Energy: " << nu_E0 << " weight: " << nu_weight << endl;
        cout << "Interaction mode: " << nu_interaction_mode << " Interaction type: " <<  nu_interaction_type << " CCorNC: " << nu_CC_NC << endl;
        cout << "Target: " << nu_target << " HitNuc: " << nu_HitNuc << " HitQuark: "<< nu_HitQuark << endl;
        cout << "W: " <<  nu_W << " X: " << nu_X << " Y: "<< nu_Y << " Qsqr: "<< nu_QSqr << endl;
        cout << endl;

        //Print relevant information about the generator particles in the interaction
        cout << "Generator Particles: " << endl;
        //The number of particles is given by the size of any gen_part_* vector
        int num_gen_particles = gen_part_E0->size();
        for(int i_p = 0; i_p < num_gen_particles; i_p++) {
            cout << "ID: " << gen_part_trackID->at(i_p) << " Status code: "<< gen_part_statusCode->at(i_p) <<  " Mother: " << gen_part_mother->at(i_p) << endl;
            cout << "   PDGCode: " << gen_part_PDGcode->at(i_p) << " Mass: " << gen_part_mass->at(i_p) << " Energy: " << gen_part_E0->at(i_p)<< endl;
            cout << "   Start Pos (x,y,z) : (" << gen_part_start_pos_X->at(i_p) << ", " << gen_part_start_pos_Y->at(i_p) << ", " << gen_part_start_pos_Z->at(i_p)  << ")" << endl;
            cout << "   Start Momentum (x,y,z) : (" << gen_part_P0_X->at(i_p) << ", " << gen_part_P0_Y->at(i_p) << ", " << gen_part_P0_Z->at(i_p)  << ")" << endl;
        }
        cout << endl;

        //Print relevant information about the particles in the interaction once they have been propagated by Geant4
        cout << "G4 Particles: " << endl;
        //The number of particles is given by the size of any g4_part_* vector
        int num_g4_particles = g4_part_E0->size();
        for(int i_p = 0; i_p < num_g4_particles; i_p++) {
            cout << "ID: " << g4_part_trackID->at(i_p) << " Mother: " << g4_part_mother->at(i_p) << endl;
            cout << "   PDGCode: " << g4_part_PDGcode->at(i_p) << " Mass: " << g4_part_mass->at(i_p) << " TL: " << g4_part_TL->at(i_p)  << " Initial Energy: " << g4_part_E0->at(i_p) << " Final Energy: " << g4_part_Ef->at(i_p) << endl;
            cout << "   Start Pos (x,y,z,t) : (" << g4_part_start_pos_X->at(i_p) << ", " << g4_part_start_pos_Y->at(i_p) << ", " << g4_part_start_pos_Z->at(i_p) << ", " << g4_part_start_T->at(i_p) <<  ")" << endl;
            cout << "   End Pos (x,y,z,t) : (" << g4_part_end_pos_X->at(i_p) << ", " << g4_part_end_pos_Y->at(i_p) << ", " << g4_part_end_pos_Z->at(i_p)  << ", " << g4_part_end_T->at(i_p) << ")" << endl;
            cout << "   Start Momentum (x,y,z) : (" << g4_part_P0_X->at(i_p) << ", " << g4_part_P0_Y->at(i_p) << ", " << g4_part_P0_Z->at(i_p)  << ")" << endl;
            cout << "   End Momentum (x,y,z) : (" << g4_part_Pf_X->at(i_p) << ", " << g4_part_Pf_Y->at(i_p) << ", " << g4_part_Pf_Z->at(i_p)  << ")" << endl;
            cout << "   Process: " << g4_part_process->at(i_p) << endl;
            cout << "   End_process: " << g4_part_end_process->at(i_p) << endl;
        }
        cout << endl;
        cout << endl << endl;
    }
}