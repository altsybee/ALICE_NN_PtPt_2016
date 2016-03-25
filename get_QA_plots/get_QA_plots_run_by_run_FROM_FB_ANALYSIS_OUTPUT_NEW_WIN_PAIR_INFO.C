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

const int nCentrBins = 8;

struct RunInfo
{
    int runId;
    //run-by-run histos
    TH1D *hist1D_multDistr_RunByRun_F[nCentrBins];
    TH1D *hist1D_multDistr_RunByRun_B[nCentrBins];
    TH1D *hist1D_avPtDistr_RunByRun_F[nCentrBins];
    TH1D *hist1D_avPtDistr_RunByRun_B[nCentrBins];
};

void fillBinWithRunData( const int bin, TH1D *histTarget, const TH1D* histSource, int runId )
{
    TString runName = Form( "%d", runId );
    histTarget->GetXaxis()->SetBinLabel( bin+1, runName );

    double mean = histSource->GetMean();
    double rms = histSource->GetRMS();
    histTarget->SetBinContent( bin+1, mean );
    histTarget->SetBinError( bin+1, rms );
}

void fixTitlesLabels( TH1D *hist )
{
    hist->GetYaxis()->SetLabelSize(0.055);
    hist->GetYaxis()->SetTitleSize(0.065);
    hist->GetYaxis()->SetTitleOffset(0.6);
    hist->GetYaxis()->CenterTitle();
}

void get_QA_plots_run_by_run_FROM_FB_ANALYSIS_OUTPUT_NEW_WIN_PAIR_INFO()
{
    gStyle->SetOptStat(0);

    TFile *file = new TFile( "output_FB_ANALYSIS_hist_mult_in_wins_run_by_run/output_histos_graphs_RUN_BY_RUN_NEW.root" );


    RunInfo runInfo[1000];

    TList *listKeys = file->GetListOfKeys();
    //    listKeys->ls();

    TH1D *histIter;
    TIter next(listKeys);
    int runId = -1;
    int runCounter = -1;
    int cBinCounter = 0;
    cout << "##### start iteration..." << endl;
    while ( histIter = (TH1D*)next() )
    {
        // mult hist
        TString strNameMult = histIter->GetName();
        //        hist->DrawCopy();
        //        break;

        //        cout << strNameMult << endl;

        if ( strNameMult.Contains( "multDistr_winF" ) )
        {
            TH1D *hist = (TH1D*)file->Get( strNameMult );

            //check run
            TString title = hist->GetTitle();
            if( title.Atoi() != runId )
            {
                // cout << "title.Atoi()=" << title.Atoi() << ", runCounter=" << runCounter << endl;
                runId = title.Atoi();

                for ( int i = 0; i < runCounter; i++ )
                {
                    if ( runId == runInfo[i].runId )
                    {
                        // we started to loop over runs again! (=next cBin!)
                        runCounter = -1;
                        cBinCounter++;
                        break;
                    }
                }
                runCounter++;

                runInfo[runCounter].runId = runId;
            }

            runInfo[runCounter].hist1D_multDistr_RunByRun_F[cBinCounter] = hist;
            hist->DrawCopy();
            //                    break;
        }
        else if ( strNameMult.Contains( "multDistr_winB" ) )
        {
            TH1D *hist = (TH1D*)file->Get( strNameMult );
            runInfo[runCounter].hist1D_multDistr_RunByRun_B[cBinCounter] = hist;
        }
        else if ( strNameMult.Contains( "avPtDistr_winF" ) )
        {
            TH1D *hist = (TH1D*)file->Get( strNameMult );
            runInfo[runCounter].hist1D_avPtDistr_RunByRun_F[cBinCounter] = hist;
        }
        else if ( strNameMult.Contains( "avPtDistr_winB" ) )
        {
            TH1D *hist = (TH1D*)file->Get( strNameMult );
            runInfo[runCounter].hist1D_avPtDistr_RunByRun_B[cBinCounter] = hist;

            //            cBinCounter++; // it's last hist in this cBin!
        }
    } // end of iteration over histos in file

    //    runInfo[0].hist1D_multDistr_RunByRun_F[0]->DrawCopy();


    //return;

    TH1D *hist_multF_MeanRMS[nCentrBins];
    TH1D *hist_multB_MeanRMS[nCentrBins];

    TH1D *hist_avPtF_MeanRMS[nCentrBins];
    TH1D *hist_avPtB_MeanRMS[nCentrBins];

    //additional: to see the spread of means:
    TH1D *hist_mult_variations[nCentrBins];
    TH1D *hist_avPt_variations[nCentrBins];

    for ( int cBin = 0; cBin < nCentrBins; cBin++ )
    {
        //        cout << "cBin=" << cBin << endl;

        TString strName = Form("hist_multF_MeanRMS_c%d", cBin );
        hist_multF_MeanRMS[cBin] = new TH1D( strName, ";runId;#LTn_{F}#GT #pm RMS", runCounter, -0.5, runCounter-0.5 );

        strName = Form("hist_multB_MeanRMS_c%d", cBin );
        hist_multB_MeanRMS[cBin] = new TH1D( strName, ";runId;#LTn_{B}#GT #pm RMS", runCounter, -0.5, runCounter-0.5 );

        strName = Form("hist_avPtF_MeanRMS_c%d", cBin );
        hist_avPtF_MeanRMS[cBin] = new TH1D( strName, ";runId;#LT#LTp_{T}^{F}#GT#GT #pm RMS", runCounter, -0.5, runCounter-0.5 );

        strName = Form("hist_avPtB_MeanRMS_c%d", cBin );
        hist_avPtB_MeanRMS[cBin] = new TH1D( strName, ";runId;#LT#LTp_{T}^{B}#GT#GT #pm RMS", runCounter, -0.5, runCounter-0.5 );

        //additional: to see the spread of means:
        strName = Form("hist_mult_variations_c%d", cBin );
        hist_mult_variations[cBin] = new TH1D( strName, ";#LTn#GT;n runs", 1000, -0.5, 1000-0.5 );

        strName = Form("hist_mult_variations_c%d", cBin );
        hist_avPt_variations[cBin] = new TH1D( strName, ";#LT#LTp_{T}#GT#GT;n runs", 4000, 0.2, 2.0 );

        if(1)for ( int run = 0; run < runCounter; run++ )
        {
            //            cout << "run=" << run << endl;

            if ( runInfo[run].hist1D_multDistr_RunByRun_F[cBin]->GetEntries() < 10000 )
                continue;

            //mult
            fillBinWithRunData( run, hist_multF_MeanRMS[cBin]
                                , runInfo[run].hist1D_multDistr_RunByRun_F[cBin], runInfo[run].runId );
            fillBinWithRunData( run, hist_multB_MeanRMS[cBin]
                                , runInfo[run].hist1D_multDistr_RunByRun_B[cBin], runInfo[run].runId );
            //av pT
            fillBinWithRunData( run, hist_avPtF_MeanRMS[cBin]
                                , runInfo[run].hist1D_avPtDistr_RunByRun_F[cBin], runInfo[run].runId );
            fillBinWithRunData( run, hist_avPtB_MeanRMS[cBin]
                                , runInfo[run].hist1D_avPtDistr_RunByRun_B[cBin], runInfo[run].runId );

        }

        for ( int hBin = 0; hBin < runCounter; hBin++ )
        {
            if ( hist_multF_MeanRMS[cBin]->GetBinContent(hBin+1) > 0 )
                hist_mult_variations[cBin]->Fill( hist_multF_MeanRMS[cBin]->GetBinContent(hBin+1) );
            if ( hist_multB_MeanRMS[cBin]->GetBinContent(hBin+1) > 0 )
                hist_mult_variations[cBin]->Fill( hist_multB_MeanRMS[cBin]->GetBinContent(hBin+1) );

            if ( hist_avPtF_MeanRMS[cBin]->GetBinContent(hBin+1) > 0 )
                hist_avPt_variations[cBin]->Fill( hist_avPtF_MeanRMS[cBin]->GetBinContent(hBin+1) );
            if ( hist_avPtB_MeanRMS[cBin]->GetBinContent(hBin+1) > 0 )
                hist_avPt_variations[cBin]->Fill( hist_avPtB_MeanRMS[cBin]->GetBinContent(hBin+1) );
        }

    }

    TCanvas *canv_mean_rms_run_by_run = new TCanvas( "canv_mean_rms_run_by_run" , "canv_mean_rms_run_by_run"
                                                     , 10, 10, 1000, 1200 ); //800,450 );
    //    tuneCanvas(canv_mean_rms_run_by_run);

    canv_mean_rms_run_by_run->Divide(2,4);

    //    tuneHist1D( hist_multF_MeanRMS[0] );


    for ( int cBin = 0; cBin < nCentrBins; cBin++ )
    {
        canv_mean_rms_run_by_run->cd(cBin+1);

        // !!!!! SELECT WHAT TO DRAW!
        //        TH1D *hist = hist_multF_MeanRMS[cBin];
        TH1D *hist = hist_avPtF_MeanRMS[cBin];
        fixTitlesLabels( hist );

        hist->SetMarkerStyle( 24 );
        hist->SetMarkerSize( 0.7 );
        hist->SetMarkerColor( kBlue );
        hist->SetLineColor( kBlue );

        //        hist->GetYaxis()->SetRangeUser( 0.45, 0.78 );


        TString strBin = Form( "centr bin %d-%d%%", cBin*10, (cBin+1)*10 );
        //TLatex *tex = new TLatex(0.25,0.28, strBin);
        hist->SetTitle( strBin );
        //        drawTex(tex);

        hist->DrawCopy();
    }


    // QA variations of means in cBins
    TCanvas *canv_mean_rms_run_by_run = new TCanvas( "canv_mean_rms_run_by_run" , "canv_mean_rms_run_by_run"
                                                     , 10, 10, 1000, 1200 ); //800,450 );
    //    tuneCanvas(canv_mean_rms_run_by_run);

    canv_mean_rms_run_by_run->Divide(2,4);
    for ( int cBin = 0; cBin < nCentrBins; cBin++ )
    {
        canv_mean_rms_run_by_run->cd(cBin+1)->SetLogy();

        // !!!!! SELECT WHAT TO DRAW!
                TH1D *hist = hist_mult_variations[cBin];
//        TH1D *hist = hist_avPt_variations[cBin];
        fixTitlesLabels( hist );
        hist->GetXaxis()->SetLabelSize(0.055);
        hist->GetXaxis()->SetTitleSize(0.065);
        hist->GetXaxis()->SetTitleOffset(0.7);
        hist->GetXaxis()->CenterTitle();

        hist->SetMarkerStyle( 24 );
        hist->SetMarkerSize( 0.7 );
        hist->SetMarkerColor( kBlue );
        hist->SetLineColor( kBlue );

//        hist->GetXaxis()->SetRangeUser( 0.6, 0.70 );
        TString strBin = Form( "centr bin %d-%d%%", cBin*10, (cBin+1)*10 );
        hist->SetTitle( strBin );

        hist->DrawCopy();

    }


    return;
}











