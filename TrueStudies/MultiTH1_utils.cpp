//
// Created by Luis Pelegrina Gutiérrez on 30/6/24.
//
#include "../Includes.h"
double XtoPad(double x);
double YtoPad(double x);

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin,
                     Float_t vSpacing, Float_t hSpacing,
                     Float_t PR, Float_t PL,
                     Float_t PT, Float_t PB,
                     Float_t PL_SCorr,
                     Float_t PB_SCorr,
                     bool XAxisShared, bool YAxisShared)
{
    if (!C) return;

    // Setup Pad layout:
    Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

    Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

    Float_t vposd,vposu,vmard,vmaru,vfactor;
    Float_t hposl,hposr,hmarl,hmarr,hfactor;

    double LAdjust = 1 - hStep/(hStep + PL_SCorr)*(1 - PL);
    double BAdjust = 1 - vStep/(vStep + PB_SCorr)*(1 - PB);

    for (Int_t i=0;i<Nx;i++) {

        if (i==0) {
            hposl = lMargin;
            if (YAxisShared) hposl = lMargin - PL_SCorr;
            hposr = lMargin + hStep;
        } else {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
        }

        for (Int_t j=0;j<Ny;j++) {

            if (j==0) {
                vposd = bMargin;
                if (XAxisShared) vposd = bMargin - PB_SCorr;
                vposu = bMargin + vStep;
            } else {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep;
            }

            C->cd(0);

            auto name = TString::Format("pad_%d_%d",i,j);
            auto pad = (TPad*) C->FindObject(name.Data());
            if (pad) delete pad;
            pad = new TPad(name.Data(),"",hposl,vposd,hposr,vposu);

            pad->SetLeftMargin(PL);
            if ((YAxisShared)&& (i == 0)) pad->SetLeftMargin(LAdjust);
            pad->SetRightMargin(PR);
            pad->SetBottomMargin(PB);
            if ((XAxisShared)&& (j == 0))  pad->SetBottomMargin(BAdjust);
            pad->SetTopMargin(PT);

            pad->SetFrameBorderMode(0);
            pad->SetBorderMode(0);
            pad->SetBorderSize(0);

            pad->Draw();
        }
    }
}


double XtoPad(double x)
{
    double xl,yl,xu,yu;
    gPad->GetPadPar(xl,yl,xu,yu);
    double pw = xu-xl;
    double lm = gPad->GetLeftMargin();
    double rm = gPad->GetRightMargin();
    double fw = pw-pw*lm-pw*rm;
    return (x*fw+pw*lm)/pw;
}

double YtoPad(double y)
{
    double xl,yl,xu,yu;
    gPad->GetPadPar(xl,yl,xu,yu);
    double ph = yu-yl;
    double tm = gPad->GetTopMargin();
    double bm = gPad->GetBottomMargin();
    double fh = ph-ph*bm-ph*tm;
    return (y*fh+bm*ph)/ph;
}

void PlotMultiTH1(std::vector<TH1*> TH1Vec[], struct MultiTH1 mTH1,  int NPlots, bool Show) {

    int Nx = mTH1.nMultiBinx;
    int Ny = mTH1.nMultiBiny;
    bool XAxisShared = mTH1.XAxisShared;
    bool YAxisShared = mTH1.YAxisShared;
    struct DrawSettings Ds = mTH1.DrawSet;

    TCanvas *cj = new TCanvas();

    double lMargin = Ds.lMargin;
    double rMargin = Ds.rMargin;
    double bMargin = Ds.bMargin;
    double tMargin = Ds.tMargin;
    double vSpacing = Ds.vSpacing;
    double hSpacing = Ds.hSpacing;
    double PR = Ds.PadRightMargin;
    double PL = Ds.PadLeftMargin;
    double PT = Ds.PadToptMargin;
    double PB = Ds.PadBottomMargin;

    double PL_NS = Ds.PadLeftMargin_NotShared ;
    double lMargin_NS = Ds.lMargin_NotShared;

    double PL_SCorr =  Ds.PadLeftMarginSharedCorrection;
    double PB_SCorr = Ds.PadBottomMarginSharedCorrection;

    if(mTH1.Normalize) {
        PL_NS = Ds.PadLeftMargin_NotShared_Norm;
        lMargin_NS = Ds.lMargin_NotShared_Norm;
    }
    if(!YAxisShared) PL = PL_NS;
    if(!XAxisShared) PB = Ds.PadBottomMargin_NotShared;
    if(!YAxisShared) lMargin = lMargin_NS;


    Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
    Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

    // Canvas setup
    CanvasPartition(cj,Nx,Ny,lMargin,rMargin,bMargin,tMargin, vSpacing, hSpacing, PR, PL, PT, PB, PL_SCorr, PB_SCorr, XAxisShared, YAxisShared);
    TPad *pad[Nx][Ny];



    pad[0][0] = (TPad*) cj->FindObject(TString::Format("pad_%d_%d",0,0).Data());
    double normalwidth = cj->GetWw()*pad[0][0]->GetWNDC();
    double normalheigth = cj->GetWh()*pad[0][0]->GetHNDC();
    if(Nx > 1) {
        pad[1][0] = (TPad*) cj->FindObject(TString::Format("pad_%d_%d",1,0).Data());
        normalwidth = cj->GetWw()*pad[1][0]->GetWNDC();
        normalheigth = cj->GetWh()*pad[1][0]->GetHNDC();
    }
    if(Ny > 1) {
        pad[0][1] = (TPad*) cj->FindObject(TString::Format("pad_%d_%d",0,1).Data());
        normalwidth = cj->GetWw()*pad[0][1]->GetWNDC();
        normalheigth = cj->GetWh()*pad[0][1]->GetHNDC();
    }
    if((Nx > 1) &&(Ny > 1)) {
        pad[1][1] = (TPad*) cj->FindObject(TString::Format("pad_%d_%d",1,1).Data());
        normalwidth = cj->GetWw()*pad[1][1]->GetWNDC();
        normalheigth = cj->GetWh()*pad[1][1]->GetHNDC();
    }

    for (int i = 0; i<Nx; i++) {
        for (int j = Ny -1; j >= 0; j--) {
            int Index = i+j*Nx;
            cj->cd(0);

            pad[i][j] = (TPad*) cj->FindObject(TString::Format("pad_%d_%d",i,j).Data());
            pad[i][j]->Draw();
            pad[i][j]->SetFillStyle(4000);
            pad[i][j]->SetFrameFillStyle(4000);
            pad[i][j]->cd();

            TH1F *hFrame[NPlots];
            for(int ik = 0; ik <NPlots; ik++) {
                hFrame[ik] = (TH1F*)  TH1Vec[ik].at(Index)->Clone(TString::Format("h_%d_%d",i,j).Data());
            }

            // y axis range
            double YMax = hFrame[0]->GetMaximum();
            if(YAxisShared) {
                for(int ik = 0; ik <NPlots; ik++) {
                    for (int ii = 0; ii<Nx; ii++) {
                        int Index2 = ii+j*Nx;
                        if(TH1Vec[ik].at(Index2)->GetMaximum() > YMax) YMax = TH1Vec[ik].at(Index2)->GetMaximum();
                    }
                }
            }



            int nMultiBins = mTH1.nMultiBinx * mTH1.nMultiBiny;
            double BinDividerPas = (mTH1.MultiBinInfo.UpBin - mTH1.MultiBinInfo.LowBin) / nMultiBins;

            for(int ik = 0; ik < NPlots; ik++) {
                pad[i][j]->cd();
                hFrame[ik]->SetMaximum(Ds.YScaling * YMax);
                hFrame[ik]->SetMinimum(0);

                double pij_height = cj->GetWh()*pad[i][j]->GetHNDC();
                double pij_width = cj->GetWw()*pad[i][j]->GetWNDC();

                // Format for y axis
                hFrame[ik]->GetYaxis()->SetLabelFont(43);
                hFrame[ik]->GetYaxis()->SetLabelOffset(0.02);
                hFrame[ik]->GetYaxis()->SetLabelSize(16*normalwidth/pij_width);
                hFrame[ik]->GetYaxis()->SetTitleSize(16*normalwidth/pij_width);
                hFrame[ik]->GetYaxis()->SetNdivisions(Ds.N_divisionsY);

                // Format for x axis

                hFrame[ik]->GetXaxis()->SetLabelFont(43);
                hFrame[ik]->GetXaxis()->SetLabelOffset(0.02);
                hFrame[ik]->GetXaxis()->SetLabelSize(16*normalheigth/pij_height);
                hFrame[ik]->GetXaxis()->SetTitleSize(16*normalheigth/pij_height);
                hFrame[ik]->GetXaxis()->SetNdivisions(Ds.N_divisionsX);
                // Draw cloned histogram with individual setting

                hFrame[ik]->SetStats(0);
                hFrame[ik]->SetTitle(";;");
                hFrame[ik]->SetFillColorAlpha(kBlue + 2, 0.1);
                hFrame[ik]->SetLineColor(kBlue +2);
                if(ik == 0) {
                    hFrame[ik]->Draw("hist");
                } else {
                    hFrame[ik]->Draw("hist same");
                }


                cj->cd(0);



                TBox* box = new TBox(i*(hSpacing + hStep) + lMargin + hStep*Ds.LegendXPos+0.06,  j*(vSpacing + vStep) + bMargin + vStep*(0.73-ik*0.07)  , i*(hSpacing + hStep) + lMargin + hStep*Ds.LegendXPos+0.08, j*(vSpacing + vStep) + bMargin + vStep*(0.73-ik*0.07) );
                box->SetLineColor(kBlue +2);
                box->SetFillStyle(0);
                box->Draw();

                TString sNE=  TString("Events: ");
                double Integral = hFrame[ik]->GetSumOfWeights();
                if(!mTH1.Normalize) sNE += TString::Format("%4.0f",Integral);
                if(mTH1.Normalize) sNE += TString::Format("%6.0f",Integral);
                TLatex *tNE = new TLatex(i*(hSpacing + hStep) + lMargin + hStep*Ds.LegendXPos,j*(vSpacing + vStep) + bMargin + vStep*(0.73-ik*0.1), sNE);
                tNE->SetTextSize(0.025);
                tNE->Draw();
                cj->Update();

            }

            TString sNom = " [" + TString::Format("%2.1f",BinDividerPas*Index) + ", " + TString::Format("%2.1f",BinDividerPas*(Index+1)) + "]";
            TLatex *tNom = new TLatex(i*(hSpacing + hStep) + lMargin + hStep*Ds.LegendXPos,j*(vSpacing + vStep) + bMargin + vStep*0.87, TString(mTH1.MultiBinInfo.Title)  + " #in " + sNom + TString(mTH1.MultiBinInfo.Unit));
            tNom->SetTextSize(0.03);
            tNom->Draw();
            cj->Update();

        }
    }
    cj->cd(0);

    TLatex *tPmu = new TLatex(lMargin + (hStep + hSpacing)*Nx/2 ,0.035, mTH1.LocalBinInfo.Title.c_str());
    tPmu->SetTextSize(0.06);
    tPmu->Draw();
    cj->Update();

    TLatex *tEvents = new TLatex(0.035,0.5+(bMargin-tMargin)/2,"Events");
    if(mTH1.Normalize) tEvents->SetX(0.035);
    tEvents->SetTextSize(0.06);
    tEvents->SetTextAngle(90);
    tEvents->Draw();
    cj->Update();


    //Save the file
    string FolderName = "MultiBin_" + mTH1.MultiBinInfo.BinDataType;
    string MkdirRuta = "mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MuPiTrue/MultiTH1/" + FolderName;
    string MkdirRutaFinal = "mkdir /Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MuPiTrue/MultiTH1/" + FolderName + "/" + mTH1.FinalStateCut;
    gSystem->Exec(MkdirRuta.c_str());
    gSystem->Exec(MkdirRutaFinal.c_str());



    string xShared = "xS";
    if(!mTH1.XAxisShared) xShared = "xNotS";

    string yShared = "yS";
    if(!mTH1.YAxisShared) yShared = "yNotS";

    string sNorm = "N";
    if(!mTH1.Normalize) sNorm = "NotN";

    string ruta = "/Users/luispelegrinagutierrez/Desktop/Doctorado/Neutrino_Physics_Analisis/CC1PiProject/Graphs/MuPiTrue/MultiTH1/" + FolderName + "/" + mTH1.FinalStateCut
                  + "/"+ mTH1.LocalBinInfo.FillDataType + "_"
                  + to_string(mTH1.nMultiBinx) + "x" + to_string(mTH1.nMultiBiny) + "_"
                  + sNorm + "_"
                  + xShared + "_" + yShared;

    string rutaPdf = ruta +".pdf";
    string rutaRoot = ruta +".root";
    cj->SaveAs(rutaPdf.c_str());
    cj->SaveAs(rutaRoot.c_str());

    if (!Show) cj->Close();
}

void InitializeTH1Vec(std::vector<TH1*>& h, int Identifier, int Size, struct MultiTH1 mTH){
    for (int j = 0; j<Size; j++) {
        string title = "h" + to_string(Identifier) + to_string(j);
        h.push_back(new TH1D(title.c_str()," ", mTH.LocalBinInfo.nBins, mTH.LocalBinInfo.LowBin, mTH.LocalBinInfo.UpBin));
    }
    return;
}
