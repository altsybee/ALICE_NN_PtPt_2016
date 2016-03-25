#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>


#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1I.h"
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

//#include "SupplementaryClasses.h"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;

#include "../utils.C"

int getRunNumberFromBinLabel( TH1D *hist, int binId )
{
    TString str = Form( "%s", hist->GetXaxis()->GetBinLabel( binId ) );
    return str.Atoi();
}

TString stripRunNumber( TString s )
{
    TString runName = s;
    //    runName.Remove(0,16);
    //    runName.Remove( runName.Length()-5, 5);
    while ( !runName.IsDigit() )
        runName.Remove(0,1);
    return runName;
    //    return runName.Atoi();
}

void get_QA_plots_run_by_run_FROM_FB_ANALYSIS_OUTPUT() // TString inputFileName = "MergedOutput.root")
{
    //    gStyle->SetOptStat(0);

    const int nFiles = 3;
    TFile *files[nFiles];
    files[0] = new TFile( "output_FB_ANALYSIS_hist_mult_in_wins_run_by_run/histos_QA_mult_in_eta_wins_run-by-run_HIJING_reco_vertZ_8.root" );
    files[1] = new TFile( "output_FB_ANALYSIS_hist_mult_in_wins_run_by_run/histos_QA_mult_in_eta_wins_run-by-run_MFminus_vertZ_8.root" );
    files[2] = new TFile( "output_FB_ANALYSIS_hist_mult_in_wins_run_by_run/histos_QA_mult_in_eta_wins_run-by-run_MFplus_vertZ_8.root" );
    //    TFile *file = new TFile( "output_FB_ANALYSIS_hist_mult_in_wins_run_by_run/histos_QA_mult_in_eta_wins_run-by-run_vertZ_5_8.root" );

    int nRunsInFiles[nFiles];
    TH1D *histMultFileRun[nFiles][1000];
    TH1D *histVertZFileRun[nFiles][1000];
    TH1D *histNeventsInRuns[nFiles];

    TString strFilePrefix[nFiles] = { "MC", "Data1", "Data2" };

    for ( Int_t fId = 0; fId < nFiles; fId++ )
    {
        TFile *file = files[fId];
        if ( !file )
        {
            cout << "NO INPUT FILE!" << endl;
            return;
        }

        TList *listKeys = file->GetListOfKeys();
        listKeys->ls();

        const int nTrees = listKeys->GetSize()/2; //45;//30;
        cout << ">>> nTrees = " << nTrees << endl;


        //prepare canvas for QA
        TString str_canv_QA_hist = Form( "canv_QA_hist_file%d", fId );
        TCanvas *canv_QA_hist = new TCanvas( str_canv_QA_hist, str_canv_QA_hist, 20+100*fId,50+100*fId,700,600 );
        tuneCanvas(canv_QA_hist);
        gPad->SetGridx();


        //prepare hist with n events in runs
        TString str_name_histNumberOfEvents = Form( "fHistNumberOfEvents_file%d", fId );
        TH1D *fHistNumberOfEvents = new TH1D( str_name_histNumberOfEvents,";;N_{events}", nTrees, -0.5, nTrees-0.5 );
        histNeventsInRuns[fId] = fHistNumberOfEvents;

        int runCounter = 0;
        bool firstHist = true;
        TString name;
        TH1D *hist;
        TIter next(listKeys);
        while ( hist = (TH1D*)next() )
        {
            // mult hist
            TString strNameMult = hist->GetName();
            cout << strNameMult << endl;

            TString strRunNum = stripRunNumber(strNameMult);

            TH1D *hist1D_multInWin = (TH1D*)file->Get(strNameMult);
            hist1D_multInWin->SetName( Form( "histMult_file%s_run_%s", strFilePrefix[fId].Data(), strRunNum.Data() ) );
            histMultFileRun[fId][runCounter] = hist1D_multInWin;

            //vertex hist
            hist = (TH1D*)next();
            TString strNameVertex = hist->GetName();
            cout << strNameVertex << endl;
            TH1D *hist1D_vertexZ = (TH1D*)file->Get(strNameVertex);
            hist1D_vertexZ->SetName( Form( "hist1D_vertexZ_file%s_run_%s", strFilePrefix[fId].Data(), strRunNum.Data() ) );
            histVertZFileRun[fId][runCounter] = hist1D_vertexZ;

            //        continue;

            int nEventsInRun = hist1D_vertexZ->GetEntries();
            if(1)if ( nEventsInRun < 100000
                      && fId > 0
                      )
            {
                //                runCounter++;
                continue;
            }


            //fill bin with nEvents for current run
            fHistNumberOfEvents->GetXaxis()->SetBinLabel( runCounter+1, strRunNum );
            fHistNumberOfEvents->SetBinContent( runCounter+1, nEventsInRun );


            cout << hist1D_multInWin->GetName() << endl;
            cout << hist1D_vertexZ->GetName() << ", " << nEventsInRun << endl;


            //modify bin content: set errors, scale
            for ( int bin = 0; bin < hist1D_multInWin->GetNbinsX(); bin++ )
            {
                double binContent = hist1D_multInWin->GetBinContent(bin+1);
                //            cout << "bin=" << bin+1 << ", binContent=" << binContent << endl;

                hist1D_multInWin->SetBinError( bin+1, sqrt( binContent ) );

                TString label = Form( "%.2f", -0.75+0.1*bin );
                hist1D_multInWin->GetXaxis()->SetBinLabel( bin+1, label );

            }
            hist1D_multInWin->Scale(1./nEventsInRun);


            int color = runCounter < 20 ? kOrange-9+runCounter : kPink-9+runCounter-20;
            if ( runCounter >= 40 )
            {
                if ( runCounter < 48 )
                    color = kRed-4+runCounter-40;
                else
                    color = kGray;
            }


            hist1D_multInWin->SetLineColor(color);

            hist1D_multInWin->GetYaxis()->SetRangeUser(0, 40 );
            //        hist1D_multInWin->DrawNormalized( i==0 ? "" : "same" );
            hist1D_multInWin->DrawCopy( firstHist ? "" : "same" );
            firstHist = false;

            runCounter++;
        } // end of iteration over histos in file
        nRunsInFiles[fId] = runCounter;


        //prepare canvas nEvents in runs
        TString str_canv_nEvents_inRunsName = Form( "canv_nEvents_inRuns_file%d", fId );
        TCanvas *canv_nEvents_inRuns = new TCanvas( str_canv_nEvents_inRunsName ,str_canv_nEvents_inRunsName
                                                    , 300+100*fId, 50+100*fId, 700,600 );
        tuneCanvas(canv_nEvents_inRuns);

        tuneHist1D(fHistNumberOfEvents);
        fHistNumberOfEvents->GetXaxis()->SetLabelSize(0.03);
        fHistNumberOfEvents->DrawCopy();
    } //end of taking histos from files


    //now we have arrays of histos for runs from several files!

    // compare nEvents run-by-run
    TH1D* histCompareNevents_MC = (TH1D*)histNeventsInRuns[0]->Clone( "histCompare_nEvents_runByRun_MC" );
    TH1D* histCompareNevents_Data = (TH1D*)histNeventsInRuns[0]->Clone( "histCompare_nEvents_runByRun_Data" );

    histCompareNevents_MC->Reset();
    histCompareNevents_Data->Reset();

    TCanvas *canv_Ratio_Of_Mults_MC_Data_run_by_run = new TCanvas( "canv_Ratio_Of_Mults_MC_Data_run_by_run" , "canv_Ratio_Of_Mults_MC_Data_run_by_run"
                                                            , 0, 0, 700,600 );
    tuneCanvas(canv_Ratio_Of_Mults_MC_Data_run_by_run);

    bool firstRatioOfMultsDrawing = true;

    //loop over bins with MC data, then match them with runs from real data:
    for ( int bin = 0; bin < histNeventsInRuns[0]->GetNbinsX(); bin++ )
    {
        double binContent = histNeventsInRuns[0]->GetBinContent(bin+1);

        int runNumberMC = getRunNumberFromBinLabel( histNeventsInRuns[0], bin+1 );
        cout << "bin = " << bin+1 << ", binContent = " << binContent << ", runNumberMC = " << runNumberMC;// << endl;

        //data runs loop:
        for ( int k = 0; k < 2; k++ )
        {
            int dataBlock = k+1; //to take histos with data
            for ( int binData = 0; binData < histNeventsInRuns[dataBlock]->GetNbinsX(); binData++ )
            {
                double binContentData = histNeventsInRuns[dataBlock]->GetBinContent(binData+1);
                int runNumberData = getRunNumberFromBinLabel( histNeventsInRuns[dataBlock], binData+1 );
                if ( runNumberData == runNumberMC ) //found corresponding runs in MC and Data
                {
                    cout << " >>> matched with data: binData =" << binData+1
                         << ", binContentData = " << binContentData << ", runNumberData = " << runNumberData;// << endl;

                    histCompareNevents_MC->SetBinContent( bin+1, binContent );
                    histCompareNevents_Data->SetBinContent( bin+1, binContentData );

                    // ######### take mult in eta wins for THIS run in MC and Data
                    if(1)for ( int rMC = 0; rMC < nRunsInFiles[0]; rMC++ )
                    {
                        //cout << " #####  mult hist: MC run =" << stripRunNumber( histMultFileRun[0][rMC]->GetName() ).Atoi() << " ";
                        if ( runNumberMC != stripRunNumber( histMultFileRun[0][rMC]->GetName() ).Atoi() )
                            continue;
                        //... mult hist for MC run found, now search for correspoding Data run and hist
                        for ( int rData = 0; rData < nRunsInFiles[dataBlock]; rData++ )
                        {
                            if ( runNumberData == stripRunNumber( histMultFileRun[dataBlock][rData]->GetName() ).Atoi() )
                            {
                                cout << " #####  matched mult hist: MC run =" << runNumberMC << ", Data run = " << runNumberData; // << endl;
                                // URA!!! Have mult histos from MC and Data for the same run!
                                histMultFileRun[dataBlock][rData]->Divide( histMultFileRun[0][rMC] );
                                histMultFileRun[dataBlock][rData]->GetYaxis()->SetTitle("ratio MC to Data");
                                histMultFileRun[dataBlock][rData]->GetYaxis()->SetRangeUser(0.5, 1.5);
                                histMultFileRun[dataBlock][rData]->DrawCopy( firstRatioOfMultsDrawing ? "" : "same" );
                                firstRatioOfMultsDrawing = false;


                                break;
                            }
                        }
                    }
                    break;
                }
            }
        }
        cout << endl;
    }

    // draw comparison nEvents run-by-run in MC and Data
    TCanvas *canv_compare_nEvents_run_by_run = new TCanvas( "canv_compare_nEvents_run_by_run" , "canv_compare_nEvents_run_by_run"
                                                            , 100, 250, 700,600 );
    tuneCanvas(canv_compare_nEvents_run_by_run);

    histCompareNevents_MC->SetLineColor( kBlue );
    histCompareNevents_Data->SetLineColor( kRed );

    histCompareNevents_MC->DrawCopy();
    histCompareNevents_Data->DrawCopy("same");

    //RATIO nEvents run-by-run
    TCanvas *canv_compare_nEvents_run_by_run_RATIO = new TCanvas( "canv_compare_nEvents_run_by_run_RATIO" , "canv_compare_nEvents_run_by_run_RATIO"
                                                                  , 200, 250, 700,600 );
    tuneCanvas(canv_compare_nEvents_run_by_run_RATIO);

    histCompareNevents_MC->GetYaxis()->SetTitle("ratio of n events MC to Data");
    histCompareNevents_MC->Divide( histCompareNevents_Data );
    histCompareNevents_MC->DrawCopy();

    return;
}











