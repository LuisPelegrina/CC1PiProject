
#include "../Includes.h"
#include "Resonance_utils.cpp"
#include "MultiTH1_utils.cpp"
bool save_graphs = true;
bool Normalize = true;


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

void Double_diff_xsec()
{
    GenerateDictionaries();

    Cut_Parameters cut_p;

    gSystem->Exec("rm -fr /Users/luispelegrinagutierrez//Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MuPiTrue/MultiTH1/*");

    //Declare the variables
    string strRuta;

    const int num_trees = 1;

    TTree *tree[num_trees];
    TFile *input[num_trees];


    Slice *slice = nullptr;

    for(int i_t = 0; i_t < num_trees; i_t++) {
        if (i_t == 0) strRuta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Data/processed_data/processed_data_81k.root";
        input[i_t] = new TFile(strRuta.c_str());
        tree[i_t] = (TTree *) input[i_t]->Get("tree");
        tree[i_t]->SetBranchAddress("slice", &slice);
    }

    //STRUCTS CALCULATIONS
    const int kNumMultiTH1 = 2;
    int kNTypes = 1;

    struct DrawSettings DsP;
    struct DrawSettings DsTheta;
    DsTheta.N_divisionsX = 507;
    DsTheta.YScaling = 1.5;
    DsTheta.LegendXPos = 0.6;

    struct MultiBinInformation mBin_Thetamu = {"ThetaMu", "#theta_{#mu}","#circ", 0, 180};
    struct MultiBinInformation mBin_Pmu = {"PMu", "P_{#mu}", "GeV" ,0, 5};

    struct LocalBinInformation lBin_Thetamu = {"ThetaMu", "#theta_{#mu} [#circ]", 90, 0, 180};
    struct LocalBinInformation lBin_Pmu = {"PMu", "P_{#mu} [GeV]", 40, -0.20, 5};
    struct LocalBinInformation lBin_Pnu = {"PNu", "P_{#nu} [GeV]", 55, -0.25, 7};

    struct MultiTH1 mTH1[kNumMultiTH1] = {
            {3, 3, true, false, "CC1Pi", mBin_Thetamu, lBin_Pmu, DsP, true},
            {3, 3, true, false, "CC1Pi", mBin_Thetamu, lBin_Pnu, DsP, true}
    };

    //Map for data filling
    map<string, double> Data;


    double CurrentFillData;
    double CurrentBinData;

    //Inizialize Histograms
    //Inizialize Histograms
    std::vector<TH1*> h[kNumMultiTH1][kNTypes];
    double BinDividerPas[kNumMultiTH1];
    int nMultiBins[kNumMultiTH1];
    for (int i = 0; i < kNTypes; i++) {
        for (int iMulti = 0; iMulti<kNumMultiTH1;iMulti++) {
            nMultiBins[iMulti] = mTH1[iMulti].nMultiBinx * mTH1[iMulti].nMultiBiny;
            BinDividerPas[iMulti] = (mTH1[iMulti].MultiBinInfo.UpBin - mTH1[iMulti].MultiBinInfo.LowBin) / nMultiBins[iMulti];
            InitializeTH1Vec(h[iMulti][i], iMulti+100*i, nMultiBins[iMulti], mTH1[iMulti]);
        }
    }

    int ii = 0;
    for(int i_t = 0; i_t < num_trees; i_t++) {
        int num_entries = tree[i_t]->GetEntries();

        for (int i_e = 0; i_e < num_entries; ++i_e) {
            tree[i_t]->GetEntry(i_e);
            if (i_e % 100 == 0) cout << "Entry:" << i_e << endl;
            double w = 1;
            TrueInteraction true_interaction = slice->true_interaction;

            if(slice->slice_ID != 0) continue;
            if(!true_interaction.is_selected_final_state("CC1Pi")) continue;

            std::vector<GeneratorParticle> muon_vec = true_interaction.get_generator_primary_p(13);
            if(muon_vec.size() != 1) cout << "ERROR" << endl;
            GeneratorParticle muon = muon_vec.at(0);
            Data["PMu"] = muon.P0.Mag();
            Data["PNu"] = true_interaction.E0_nu;
            Data["ThetaMu"] = muon.P0.Angle(TVector3(0,0,01))*360/(2*TMath::Pi());

            for (int iMulti = 0; iMulti<kNumMultiTH1;iMulti++) {

                if(mTH1[iMulti].Normalize) {
                    w = slice->weight;
                }

                CurrentFillData = Data[mTH1[iMulti].LocalBinInfo.FillDataType];
                CurrentBinData = Data[mTH1[iMulti].MultiBinInfo.BinDataType];

                //Calculate the index
                int IndexBin = int(CurrentBinData/BinDividerPas[iMulti]);
                if(IndexBin >= nMultiBins[iMulti]) IndexBin = nMultiBins[iMulti]-1;

                h[iMulti][ii].at(IndexBin)->Fill(CurrentFillData, w);
            }

        }
    }

    for (int iMulti = 0; iMulti<kNumMultiTH1;iMulti++) {
        PlotMultiTH1(h[iMulti],  mTH1[iMulti], 1, true);

        cout << "Graph: " << iMulti <<" Done "<< endl;
    }
}




