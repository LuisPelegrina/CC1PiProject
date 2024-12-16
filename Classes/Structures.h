//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//



#ifndef MUPIPROJECT_STRUCTURES_H
#define MUPIPROJECT_STRUCTURES_H

#include <iostream>
#include "TColor.h"
using namespace std;


struct Cut {
    string containment_cut = "not_def";
    string final_state_cut = "not_def";
};


struct DrawSettings{
    double lMargin = 0.12;
    double lMargin_NotShared = 0.06;
    double lMargin_NotShared_Norm = 0.06;
    double rMargin = 0.05;
    double bMargin = 0.12;
    double tMargin = 0.05;
    double vSpacing = 0.01;
    double hSpacing = 0.01;

    double PadRightMargin = 0.02;
    double PadLeftMargin = 0.02;
    double PadLeftMargin_NotShared = 0.16;
    double PadLeftMargin_NotShared_Norm = 0.23;
    double PadLeftMarginSharedCorrection = 0.04;
    double PadToptMargin = 0.04;
    double PadBottomMargin = 0.04;
    double PadBottomMargin_NotShared = 0.12;
    double PadBottomMarginSharedCorrection = 0.06;

    int N_divisionsX = 510;
    int N_divisionsY = 507;

    double YScaling = 1.2;
    double LegendXPos = 0.65;
};

struct MultiBinInformation {
    string BinDataType = "ThetaMu";
    string Title = "#theta_{#mu}";
    string Unit = "[#circ]";
    double LowBin = 0;
    double UpBin = 180;
};

struct LocalBinInformation {
    string FillDataType = "ThetaMu";
    string Title = "#theta_{#mu}";
    double nBins = 90;
    double LowBin = 0;
    double UpBin = 180;
};

struct MultiTH1 {
    int nMultiBinx = 3;
    int nMultiBiny = 3;
    bool XAxisShared = false;
    bool YAxisShared = false;
    string FinalStateCut = "CC1Pi";
    struct MultiBinInformation MultiBinInfo;
    struct LocalBinInformation LocalBinInfo;
    struct DrawSettings DrawSet;
    bool Normalize = true;
};



struct GeneratorInformation{
    string GeneratorType = "Genie";
    string FileName= "analysisOutput.root";
    EColor HistColor = kBlue;
};


#endif //MUPIPROJECT_STRUCTURES_H
//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//
