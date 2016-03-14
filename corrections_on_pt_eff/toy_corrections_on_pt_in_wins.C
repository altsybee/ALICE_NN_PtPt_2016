#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>


#include "TObject.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TString.h"

#include "../SupplementaryClasses.cxx"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;

#include "../utils.C"


void toy_corrections_on_pt_in_wins() // TString inputFileName = "MergedOutput.root")
{
    // ########## Select and initiate windows:
    //    const int nEtaWins = 3;//2;

    TFile *fEff_HIJING = new TFile( "histEffFromHIJING_pt_02_20.root" );
    TH1D *histEff = (TH1D *) fEff_HIJING->Get( "fHist2D_ptDistrInEtaPhiVsCentr_winF_reco_eta_7_phi0_pt_0_proj0_eff_clone" );

    histEff->Draw();

//    int binEff = histEff->GetXaxis()->FindBin(0.99);
//    cout << binEff << endl;
//    cout << histEff->GetBinContent(binEff) << endl;
//    return;

    const int nRuns = 22;//18;
    TH1D *fHistPt_rec[nRuns];


    TGraphErrors *gr_bCorr_truth_point = new TGraphErrors;
    TGraphErrors *gr_bCorr = new TGraphErrors;
    TGraphErrors *gr_bCorr_rec = new TGraphErrors;
    TGraphErrors *gr_bCorr_recOnlyF = new TGraphErrors;
    TGraphErrors *gr_bCorr_recOnlyB = new TGraphErrors;

    for ( int run = 0; run < nRuns; run++ )
    {
        fHistPt_rec[run] = new TH1D("fHistPt", "p_{T} distribution", 100, 0.0, 2.0);

        WinPair winFB;
        WinPair winFBrec;
        WinPair winFBrecOnlyF;
        WinPair winFBrecOnlyB;

        winFB.init(0, 100, 0, 0);
        winFBrec.init(0, 100, 0, 0);
        winFBrecOnlyF.init(0, 100, 0, 0);
        winFBrecOnlyB.init(0, 100, 0, 0);

        double ineff = 0.25*run; //run*0.05;
        // { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };//
//        double ineffLowPt = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };//run*0.1;

        gRandom->SetSeed(0);
        for ( int ev = 0; ev < 150000; ev++ )
        {
            // 27.02.16: fill wins "by hand" using info from branches:
            UShort_t nF = 0;
            UShort_t nB = 0;
            Float_t PtF = 0;
            Float_t PtB = 0;

            UShort_t nFrec = 0;
            UShort_t nBrec = 0;
            Float_t PtFrec = 0;
            Float_t PtBrec = 0;

            int nPart = gRandom->Poisson(1000);
            double ptMean = gRandom->Gaus(0.45,0.03);
//            double ptMean = gRandom->Uniform( 0.42, 0.48 );

            int pCounter = 0;
            for ( int i = 0; i < nPart; i++ )
            {

                double eta = gRandom->Uniform( -1, 1 );
                double pt = gRandom->Exp( ptMean );

                if ( pt < 0.2 || pt > 2.0 )
                    continue;

                // #### truth level
                if ( eta > 0 )
                {
                    nF++;
                    PtF += pt;
                }
                else if ( eta < 0 )
                {
                    nB++;
                    PtB += pt;
                }

                // ### apply artificial ineff:
//                float additionalIneff = (pt < 0.4 ? (1-ineff)*0.2 : 0);
//                if ( ineff > 0 )
//                    if ( gRandom->Uniform() < ineff+additionalIneff )
//                        continue;

                int binEff = histEff->GetXaxis()->FindBin( pt );
                //cout << binEff << endl;
                //cout << histEff->GetBinContent(binEff) << endl;
                float effFromHist = histEff->GetBinContent(binEff);

//                cout << "ineff*(1-effFromHist)=" << ineff*(1-effFromHist) << endl;

                if ( gRandom->Uniform() < ineff*(1-effFromHist) )
                    continue;

//                if ( gRandom->Uniform() < 1-effFromHist )
//                    continue;
//                if ( gRandom->Uniform() < ineff*0.1 )
//                    continue;


                fHistPt_rec[run]->Fill(pt);

                pCounter++;

                // #### reco level
                if ( eta > 0 )
                {
                    nFrec++;
                    PtFrec += pt;
                }
                else if ( eta < 0 )
                {
                    nBrec++;
                    PtBrec += pt;
                }



            }
            //final for truth event
            if ( nF == 0 )
                PtF = -1;
            else
                PtF /= nF;
            if ( nB == 0 )
                PtB = -1;
            else
                PtB /= nB;

            winFB.fill( 50, nF, nB, PtF, PtB );


            //final for rec event
            if ( nFrec == 0 )
                PtFrec = -1;
            else
                PtFrec /= nFrec;
            if ( nBrec == 0 )
                PtBrec = -1;
            else
                PtBrec /= nBrec;


            winFBrec.fill( 50, nFrec, nBrec, PtFrec, PtBrec );

            winFBrecOnlyF.fill( 50, nFrec, nB, PtFrec, PtB );
            winFBrecOnlyB.fill( 50, nF, nBrec, PtF, PtBrec );

//            cout << pCounter << " ";
        }
//        cout << endl;

        winFB.calcCorrCoeffs();
        winFBrec.calcCorrCoeffs();
        winFBrecOnlyF.calcCorrCoeffs();
        winFBrecOnlyB.calcCorrCoeffs();

        double bCorr = winFB.PtPt_bCorr;
        double bCorr_rec = winFBrec.PtPt_bCorr;
        double bCorr_recOnlyF = winFBrecOnlyF.PtPt_bCorr;
        double bCorr_recOnlyB = winFBrecOnlyB.PtPt_bCorr;

        cout << "bCorr = " << bCorr << ", bCorr_rec = " << bCorr_rec << endl;

        int nMean_truth = winFB.NN_nF/winFBrec.NN_Nevents + winFB.NN_nB/winFBrec.NN_Nevents;
        int nMean_rec = winFBrec.NN_nF/winFBrec.NN_Nevents + winFBrec.NN_nB/winFBrec.NN_Nevents;
        int nMean_recOnlyF = winFBrec.NN_nF/winFBrec.NN_Nevents + winFB.NN_nB/winFB.NN_Nevents;
        int nMean_recOnlyB = winFB.NN_nF/winFB.NN_Nevents + winFBrec.NN_nB/winFBrec.NN_Nevents;

        gr_bCorr->SetPoint( run, nMean_rec/*nMean_truth*/, bCorr );
        gr_bCorr->SetPointError( run, 0, 0.01 );

        gr_bCorr_rec->SetPoint( run, nMean_rec, bCorr_rec );
        gr_bCorr_rec->SetPointError( run, 0, 0.01 );

        gr_bCorr_recOnlyF->SetPoint( run, nMean_recOnlyF, bCorr_recOnlyF );
        gr_bCorr_recOnlyF->SetPointError( run, 0, 0.01 );

        gr_bCorr_recOnlyB->SetPoint( run, nMean_recOnlyB, bCorr_recOnlyB );
        gr_bCorr_recOnlyB->SetPointError( run, 0, 0.01 );

        gr_bCorr_truth_point->SetPoint( run, nMean_truth, bCorr);


        double cFact_F = (winFB.PtPt_PtF/winFB.PtPt_Nevents) / (winFBrec.PtPt_PtF/winFBrec.PtPt_Nevents);
        double cFact_B = (winFB.PtPt_PtB/winFB.PtPt_Nevents) / (winFBrec.PtPt_PtB/winFBrec.PtPt_Nevents);
        double cFact_FB = (winFB.PtPt_PtF_PtB/winFB.PtPt_Nevents) / (winFBrec.PtPt_PtF_PtB/winFBrec.PtPt_Nevents);
        double cFact_F2 = (winFB.PtPt_PtF2/winFB.PtPt_Nevents) / (winFBrec.PtPt_PtF2/winFBrec.PtPt_Nevents);

    }


    TCanvas *canv_dep_on_ineff = new TCanvas("canv_dep_on_ineff","canv_dep_on_ineff",250,150,700,600 );
    tuneCanvas(canv_dep_on_ineff);

    TGraphErrors *grFrame = new TGraphErrors;
    grFrame->SetPoint(0,0,0);
    grFrame->SetPoint(1,0,1);
    grFrame->SetPoint(2,1000,1);
    grFrame->SetPoint(3,1000,0);
    grFrame->SetMarkerColor(kWhite);

    tuneGraphAxisLabels( grFrame );
    grFrame->Draw("AP");

    //gr truth
    gr_bCorr->SetMarkerStyle(24);
    gr_bCorr->SetMarkerColor(kBlue);
    gr_bCorr->SetLineColor(kBlue);
    gr_bCorr->DrawClone("PL");

    //gr rec
    gr_bCorr_rec->SetMarkerStyle(24);
    gr_bCorr_rec->SetMarkerColor(kRed);
    gr_bCorr_rec->SetLineColor(kRed);
    gr_bCorr_rec->DrawClone("PL");


    //gr recOnlyF
    drawGraph( gr_bCorr_recOnlyF, 24, kGreen+1, "PL" );
    //gr recOnlyB
    drawGraph( gr_bCorr_recOnlyB, 24, kPink+1, "PL" );


    gr_bCorr_truth_point->SetMarkerStyle(30);
    gr_bCorr_truth_point->SetMarkerColor(kMagenta);
    gr_bCorr_truth_point->DrawClone("P");



    TF1 *fFitFunc = new TF1("fitFunc","[0]*x/(1+[0]*x)",0,1000);
    fFitFunc->SetParameter(0,0.004);
    gr_bCorr_rec->Fit( fFitFunc );//,"Q");

//    fFitFunc->SetLineColor(kOrange-5+5);
    fFitFunc->DrawCopy("same");

    TF1 *fFitFuncMod = new TF1("fFitFuncMod","[0]*(x+[1])/(1+[0]*(x+[1]))",0,1000);
    fFitFuncMod->SetParameter(0,0.004);
    fFitFuncMod->SetParameter(1,0.01);


    //pt reco
    TCanvas *canv_pt_reco = new TCanvas("canv_pt_reco","canv_pt_reco",250,150,700,600 );
    tuneCanvas(canv_pt_reco);

    for ( int run = 0; run < nRuns; run++ )
    {
        fHistPt_rec[run]->SetLineColor( kOrange-5+run );
        fHistPt_rec[run]->DrawCopy( run==0 ? "" : "same" );
    }

    return;
}
