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

    TCanvas *canv_QA_hist = new TCanvas("canv_QA_hist","canv_QA_hist",20,50,700,600 );
    tuneCanvas(canv_QA_hist);
    gPad->SetGridx();

    //    TFile *file = new TFile( "output_FB_ANALYSIS_hist_mult_in_wins_run_by_run/histos_QA_mult_in_eta_wins_run-by-run_HIJING_reco_vertZ_8.root" );
    TFile *file = new TFile( "output_FB_ANALYSIS_hist_mult_in_wins_run_by_run/histos_QA_mult_in_eta_wins_run-by-run_vertZ_8.root" );
    //    TFile *file = new TFile( "output_FB_ANALYSIS_hist_mult_in_wins_run_by_run/histos_QA_mult_in_eta_wins_run-by-run_vertZ_5_8.root" );

    TList *listKeys = file->GetListOfKeys();
    listKeys->ls();

    const int nTrees = listKeys->GetSize()/2; //45;//30;
    cout << ">>> nTrees = " << nTrees << endl;


    //prepare hist with n events in runs
    TH1D *fHistNumberOfEvents = new TH1D("fHistNumberOfEvents",";;N_{events}", nTrees, -0.5, nTrees-0.5 );

    int counter = 0;
    bool firstHist = true;
    TString name;
    TH1D *hist;
    TIter next(listKeys);
    while ( hist=(TH1D*)next() )
    {
        // mult hist
        TString strNameMult = hist->GetName();
        cout << strNameMult << endl;
        TH1D *hist1D_multInWin = (TH1D*)file->Get(strNameMult);

        //vertex hist
        hist=(TH1D*)next();
        TString strNameVertex = hist->GetName();
        cout << strNameVertex << endl;
        TH1D *hist1D_vertexZ = (TH1D*)file->Get(strNameVertex);

        //        continue;

        int nEventsInRun = hist1D_vertexZ->GetEntries();
        if ( nEventsInRun < 10000 )
        {
            counter++;
            continue;
        }


        //fill bin with nEvents for current run
        fHistNumberOfEvents->GetXaxis()->SetBinLabel( counter+1, stripRunNumber(strNameMult) );
        fHistNumberOfEvents->SetBinContent( counter+1, nEventsInRun );


        cout << hist1D_multInWin->GetName() << endl;
        cout << hist1D_vertexZ->GetName() << ", " << nEventsInRun << endl;


        for ( int bin = 0; bin < hist1D_multInWin->GetNbinsX(); bin++ )
        {
            double binContent = hist1D_multInWin->GetBinContent(bin+1);
            //            cout << "bin=" << bin+1 << ", binContent=" << binContent << endl;

            hist1D_multInWin->SetBinError( bin+1, sqrt( binContent ) );

            TString label = Form( "%.2f", -0.75+0.1*bin );
            hist1D_multInWin->GetXaxis()->SetBinLabel( bin+1, label );

        }
        hist1D_multInWin->Scale(1./nEventsInRun);


        int color = counter < 20 ? kOrange-9+counter : kPink-9+counter-20;
        if ( counter >= 40 )
        {
            if ( counter < 48 )
                color = kRed-4+counter-40;
            else
                color = kGray;
        }


        hist1D_multInWin->SetLineColor(color);

        hist1D_multInWin->GetYaxis()->SetRangeUser(0, 100 );
        //        hist1D_multInWin->DrawNormalized( i==0 ? "" : "same" );
        hist1D_multInWin->DrawCopy( firstHist ? "" : "same" );
        firstHist = false;

        counter++;
    }


    //prepare canvas nEvents in runs
    TCanvas *canv_nEvents_inRuns = new TCanvas( "canv_nEvents_inRuns","canv_nEvents_inRuns",300,50,700,600 );
    tuneCanvas(canv_nEvents_inRuns);

    tuneHist1D(fHistNumberOfEvents);
    fHistNumberOfEvents->GetXaxis()->SetLabelSize(0.03);
    fHistNumberOfEvents->DrawCopy();

    return;
}
