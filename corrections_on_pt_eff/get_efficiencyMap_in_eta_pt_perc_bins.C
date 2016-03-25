#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>


#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TRandom3.h"

#include "../SupplementaryClasses.h"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;

#include "../utils.C"


double* prepareArrayToBeUsedInRebinning( int nArr, double *arr, double firstValueInHist, double lastValueInHist )
{
    double *arrNew = new double[nArr+1];
    arrNew[0] = 0;
    for ( int i = 0; i < nArr; i++ )
        arrNew[i+1] = arr[i];
    //put last value by hand:
    arrNew[0] = firstValueInHist;
    arrNew[nArr] = lastValueInHist;
    return arrNew;
}

void stripHist1D( TH1D &hist )
{
//    int lowBin = -1;
//    int upBin = -1;
//    //find first filled bin and first empty bin after the filled bins:
//    for ( int bin = 0; bin < hist.GetNbinsX(); bin++ )
//    {
//        if ( hist.GetBinContent(bin+1) > 0 && lowBin == -1 )
//            lowBin = bin+1;
//        if ( hist.GetBinContent(bin+1) == 0 && lowBin > 0 && upBin == -1 )
//        {
//            upBin = bin+1;
//            break;
//        }
//    }

    // Right way to do it!:
    int lowBin = hist.FindFirstBinAbove(0);
    int upBin  =  hist.FindLastBinAbove(0)+1;

    double lowEdge  = hist.GetBinLowEdge(lowBin);
    double highEdge = hist.GetBinLowEdge(upBin);

    int nNewBins = upBin - lowBin;

//    cout << "lowEdge=" << lowEdge << ", highEdge=" << highEdge << ", nNewBins=" << nNewBins << endl;

    // hist without "left" and "right" empty bins
    TH1D *hist_stripped = new TH1D( Form("%s_stripped", hist.GetName()), Form("%s_stripped", hist.GetName()), nNewBins, lowEdge, highEdge );
    for ( int bin = 0; bin < hist.GetNbinsX(); bin++ )
    {
        if ( hist.GetBinContent(bin+1) == 0 )
            continue;

        //get bin content and error
        double value = hist.GetBinContent(bin+1);
        Double_t error = hist.GetBinError(bin+1);
//        cout << "value = " << value << ", error = " << error << endl;

        //set bin content and error
        hist_stripped->SetBinContent( bin+1 - lowBin + 1, value );
        hist_stripped->SetBinError( bin+1 - lowBin + 1, error );
    }

    hist = *hist_stripped;
    delete hist_stripped;

    //    return hist_stripped;

    //    double *arrNew = new double[nArr+1];
    //    arrNew[0] = 0;
    //    for ( int i = 0; i < nArr; i++ )
    //        arrNew[i+1] = arr[i];
    //    //put last value by hand:
    //    arrNew[nArr] = lastValueInHist;
    //    return arrNew;
}

void get_efficiencyMap_in_eta_pt_perc_bins() //int fileId)
{
    gStyle->SetOptStat(0);

    TFile *myFile;

    //HIJING LHC11a10a_bis 900k events AOD162
    //    myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_28_PbPb_MCAOD_LHC11a10a_bis_AOD162_TEST_Efficiency_kine_vs_reco_try5/MergedOutput.root" );
    //    myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_03_06_PbPb_MCAOD_LHC11a10a_bis_AOD162_TEST_Efficiency_kine_vs_reco/MergedOutput.root" );
    myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_03_09_PbPb_MCAOD_LHC11a10a_bis_AOD162_Efficiency_kine_vs_reco_fixedCentrBins_try2/MergedOutput.root" );



    if (!myFile)
    {
        cout << "No input file!" << endl;
        return;
    }
    myFile->ls();

    int listId = 0;//1;



    TList *listKeys = myFile->GetListOfKeys();
    cout << "going into list: " << listKeys->At(listId)->GetName() << endl;

    //return;
    myFile->cd( listKeys->At(listId)->GetName() );
    TList *listKeys2 = gDirectory->GetListOfKeys();
    TList *myTask = (TList*) gDirectory->Get( listKeys2->At(0)->GetName() );


    TCanvas *canv_QA_eff = new TCanvas("canv_QA_eff","canv_QA_eff",20,50,700,600 );
    tuneCanvas(canv_QA_eff);
    canv_QA_eff->SetGridx();

    TCanvas *canv_QA_diff = new TCanvas("canv_QA_diff","canv_QA_diff",350,50,700,600 );
    tuneCanvas(canv_QA_diff);

    const int nEtaWins = 8;//1;
    const int nPhiWins = 1;
    const int nPtBins = 1;


    TH2D *hist_kine[nPtBins][2*nEtaWins][nPhiWins];
    TH2D *hist_reco[nPtBins][2*nEtaWins][nPhiWins];

    const int nPtBinsNew = 20;
    //    double ptNewBounds[nPtBinsNew]; // array to contain the quantiles
    double ptNewBounds[nPtBinsNew] = // array to contain the quantiles BY HAND: bin accuracy up to 0.2 GeV/c !!! (=binning of pt distr!)
    { 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.38, 0.42, 0.46,
            0.50, 0.54, 0.58, 0.64, 0.70, 0.78, 0.88, 1.04, 1.30, 2.00 };
    double *ptNewBoundsMod = prepareArrayToBeUsedInRebinning( nPtBinsNew, ptNewBounds, 0.2, 2.0 );

    TGraphErrors *grCentrDep[2*nEtaWins][nPtBinsNew];

    const int nCentrBins = 20;
    const float cWidth = 5;

    const float etaMin = -0.8;
    const float etaStep = 0.1;

    // EFFICIENCY MAP
    double etaBinsFor3D[2*nEtaWins+1];
    for ( int etaW = 0; etaW < 2*nEtaWins+1; etaW++ )
        etaBinsFor3D[etaW] = etaMin + etaW*etaStep;

    double cBinsFor3D[nCentrBins+1];
    for ( int cBin = 0; cBin < nCentrBins+1; cBin++ )
        cBinsFor3D[cBin] = cBin*cWidth;

//    TH3D *hist3D_effMap = new TH3D( "hist3D_effMap", "hist3D_effMap", 2*nEtaWins, -0.8, 0.8, nPtBinsNew, 0.2, 2.0, nCentrBins, 0, 100 );
    TH3D *hist3D_effMap = new TH3D( "hist3D_effMap", "hist3D_effMap", 2*nEtaWins, etaBinsFor3D, nPtBinsNew, ptNewBoundsMod, nCentrBins, cBinsFor3D );

    for ( int etaW = 0; etaW < 2*nEtaWins; etaW++ )
    {
//        if ( etaW > 0 && etaW < 7 )
//            continue;
        for ( int phiW = 0; phiW < nPhiWins; phiW++ )
            for ( int ptW = 0; ptW < nPtBins; ptW++ )
            {
                int testWinF = -1;
                if ( etaW >= nEtaWins ) // F wins
                    testWinF = 2*nEtaWins-1-etaW;

                cout << " >>> eta win: " << (etaW < nEtaWins ? etaW : testWinF) << endl;
                TString namePostfix = Form("eta_%d_phi%d_pt_%d"
                                           , (etaW < nEtaWins ? etaW : testWinF), phiW, ptW );

                if ( etaW < nEtaWins ) // B wins
                {
                    // ### for kine:
                    hist_kine[ptW][etaW][phiW] = (TH2D*)  myTask->FindObject( Form("fHist2D_ptDistrInEtaPhiVsCentr_winB_kine_%s", namePostfix.Data() ) );
                    // ### for reco:
                    hist_reco[ptW][etaW][phiW] = (TH2D*)  myTask->FindObject( Form("fHist2D_ptDistrInEtaPhiVsCentr_winB_reco_%s", namePostfix.Data() ) );
                }
                else // F wins
                {
                    // ### for kine:
                    hist_kine[ptW][etaW][phiW] = (TH2D*)  myTask->FindObject( Form("fHist2D_ptDistrInEtaPhiVsCentr_winF_kine_%s", namePostfix.Data() ) );
                    // ### for reco:
                    hist_reco[ptW][etaW][phiW] = (TH2D*)  myTask->FindObject( Form("fHist2D_ptDistrInEtaPhiVsCentr_winF_reco_%s", namePostfix.Data() ) );
                }

                //graphs for each "rebinned" pt bin as func from centrality
                for ( int rebinnedPtBin = 0; rebinnedPtBin < nPtBinsNew; rebinnedPtBin++ )
                    grCentrDep[etaW][rebinnedPtBin] = new TGraphErrors;

                // ### prepare slices:
                for ( int slice = 0; slice < /*14*/ hist_kine[ptW][etaW][phiW]->GetNbinsX()-2; slice++ )
                {
                    //                                        if (slice>0 && slice<13)
                    //                                            continue;
                    TH1D *proj_kine = hist_kine[ptW][etaW][phiW]->ProjectionY( Form("%s_proj%d", hist_kine[ptW][etaW][phiW]->GetName(), slice ),
                                                                                         slice+1, slice+1 );
                    TH1D *proj_reco = hist_reco[ptW][etaW][phiW]->ProjectionY( Form("%s_proj%d", hist_reco[ptW][etaW][phiW]->GetName(), slice ),
                                                                                         slice+1, slice+1 );

                    //                proj_kine->DrawCopy();
                    stripHist1D( *proj_kine );
                    stripHist1D( *proj_reco );

                    if ( etaW < nEtaWins ) // B wins
                    {
//                        proj_reco->SetLineColor(kOrange-9+etaW);
                        //proj_reco->SetMarkerColor(kOrange-9+etaW);
                        proj_reco->SetLineColorAlpha( kOrange, 0.2 );
                        proj_reco->SetMarkerColorAlpha( kOrange, 0.2 );
                    }
                    else // F wins
                    {
//                        proj_reco->SetLineColor(kViolet-9+testWinF);
//                        proj_reco->SetMarkerColor(kViolet-9+testWinF);
                        proj_reco->SetLineColorAlpha( kViolet, 0.2 );
                        proj_reco->SetMarkerColorAlpha( kViolet, 0.2 );
                    }


                    if (0)
                    {
                        proj_reco->SetLineWidth(2);
                        proj_reco->DrawCopy( (slice==0 && etaW==0 && phiW==0 && ptW==0) ? "" : "same" );

                        proj_reco->SetLineColor(kRed+slice);
                        proj_reco->DrawCopy("same");

                        continue;
                    }

//                    if (etaW == 0)
//                    {
//                        proj_reco->SetLineColor(kRed);
//                        proj_reco->SetMarkerColor(kRed);
//                    }
//                    if (etaW == 7)
//                    {
//                        proj_reco->SetLineColor(kBlue);
//                        proj_reco->SetMarkerColor(kBlue);
//                    }


                    //                    proj_reco->SetLineColor(kOrange-5+slice);



                    // GET QUANTILES of pT distr:
                    if (0)
                    {
                        if(1)
                            getQuantiles( proj_kine, nPtBinsNew, ptNewBounds );

                        double *multBinCenters = new double[nPtBinsNew]; // array to contain means of 1D-histograms in bins
                        drawCanvasWithClasses(proj_kine, Form("proj_kine_%d_bins", nPtBinsNew)
                                              , nPtBinsNew, ptNewBounds, multBinCenters );
                        // !!!
                        ptNewBoundsMod = prepareArrayToBeUsedInRebinning( nPtBinsNew, ptNewBounds, 0.2, 2.0 );
                    }


                    proj_kine = (TH1D*)proj_kine->Rebin( nPtBinsNew, "proj_kine_rebinned", ptNewBoundsMod );
                    proj_reco = (TH1D*)proj_reco->Rebin( nPtBinsNew, "proj_reco_rebinned", ptNewBoundsMod );


                    //pt_bins_new1->Draw();
                    //return;

                    // ### eff:
                    canv_QA_eff->cd();
                    TH1D *histEff = (TH1D *)proj_reco->Clone( Form( "%s_eff_clone", proj_reco->GetName() ) );
                    histEff->Divide( proj_kine );

                    //                    histEff->SetBinError(slice);
                    histEff->SetMarkerStyle(20);
                    histEff->SetTitle(";p_{T}, GeV/c; efficiency");
                    tuneHist1D(histEff);
                    histEff->DrawCopy( (slice==0 && etaW==0 && phiW==0 && ptW==0) ? "" : "same" );


                    //graphs for each "rebinned" pt bin as func from centrality
                    for ( int rebinnedPtBin = 0; rebinnedPtBin < nPtBinsNew; rebinnedPtBin++ )
                    {
                        double value = histEff->GetBinContent(rebinnedPtBin+1);
                        double err = histEff->GetBinError(rebinnedPtBin+1);
                        grCentrDep[etaW][rebinnedPtBin]->SetPoint( slice, slice*cWidth + cWidth/2, value );
                        grCentrDep[etaW][rebinnedPtBin]->SetPointError( slice, 0, err );

                        // FILL EFF MAP:
                        hist3D_effMap->SetBinContent( etaW+1, rebinnedPtBin+1, slice+1, value );
                        hist3D_effMap->SetBinError( etaW+1, rebinnedPtBin+1, slice+1, err );
                    }




                    // ### diff
                    //                    canv_QA_diff->cd();
                    //                    TH1D *histDiff_reco = (TH1D *)proj_reco->Clone( Form( "%s_diff_clone", proj_reco->GetName() ) );
                    //                    TH1D *histDiff_kine = (TH1D *)proj_kine->Clone( Form( "%s_diff_clone", proj_kine->GetName() ) );
                    //                    histDiff_kine->Add( histDiff_reco, -1 );
                    //                    //                histDiff_kine->Divide( proj_kine );//, -1 );
                    //                    histDiff_kine->DrawCopy( (etaW==0 && phiW==0 && ptW==0) ? "" : "same" );
                } // end of slices
            } // end of ptW
    } // end of etaW



    //draw eff as func from centr for pt bins in several eta-wins:
    for ( int etaW = 0; etaW < 2*nEtaWins; etaW++ )
    {
        TString strCanvName = Form( "canv_eff_as_func_from_centr_etaW%d", etaW );
        TCanvas *canv_eff_as_func_from_centr = new TCanvas( strCanvName, strCanvName, 200+20*etaW,150,700,600 );
        tuneCanvas(canv_eff_as_func_from_centr);

        tuneGraphAxisLabels( grCentrDep[etaW][0] );
        grCentrDep[etaW][0]->SetTitle( ";centrality percentile;efficiency" );
        grCentrDep[etaW][0]->GetYaxis()->SetRangeUser(0,1);

        for ( int bin = 0; bin < nPtBinsNew; bin++ )
            drawGraph( grCentrDep[etaW][bin], 20, kOrange-9+bin, bin == 0 ? "APL" : "PL");

    }



    //    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistVz" );
    //    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistPt" );
    //    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistEta" );
    //TH1D *histQA = (TH1D*)  myTask->FindObject( "hist_kine_" );
    //histQA->SetName( Form( "%s_%d", histQA->GetName(), fileId ) );


    //histQA->DrawCopy();

    if (0)
    {
        TFile *fileEffMap = new TFile( "fileEffMap.root", "RECREATE" );
        hist3D_effMap->Write();
        fileEffMap->Close();
    }


    //    return;


}
