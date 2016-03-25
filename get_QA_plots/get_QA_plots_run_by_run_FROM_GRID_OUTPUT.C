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



//TFile *
void list_files( int &nFiles, const char *dirname="", const char *begins="AnalysisResults", const char *ext=".root" )
{
    gStyle->SetOptStat(0);

    //    TFile *filePtr[1000];

    //prepare canvas
    TCanvas *canv_QA_hist = new TCanvas("canv_QA_hist","canv_QA_hist",20,50,700,600 );
    canv_QA_hist->Divide( 2, 2);
//    tuneCanvas(canv_QA_hist);

    //prepare hist with n events in runs
    TH1D *fHistNumberOfEvents = new TH1D("fHistNumberOfEvents",";runId;N_{events}", 45, -0.5, 45-0.5 );


    int colors[] = { kBlack, kBlue, kRed, kMagenta, kCyan };

    TLegend *leg = new TLegend(0.65,0.1,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);



    int counter = 0;

    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files)
    {
        TSystemFile *systFile;
        TString fname;
        TIter next(files);
        while ((systFile=(TSystemFile*)next()))
        {
            fname = systFile->GetName();
            if (!systFile->IsDirectory() && fname.BeginsWith(begins) && fname.EndsWith(ext))
            {
                //                filePtr[counter] = file;
                cout << fname.Data() << endl;
                TFile *file = new TFile( Form( "%s/%s", dirname, fname.Data() ) );

                //set run name on x axis in hist
                TString runName = fname;
                runName.Remove(0,19);
                runName.Remove(6,5);
                fHistNumberOfEvents->GetXaxis()->SetBinLabel( counter+1, runName );

//                int color = colors[counter];
                int color = counter < 20 ? kOrange-9+counter : kPink-9+counter-20;
                if ( counter >= 40)
                    color = kRed-4+counter-40;

                // get QA hist from file
                //vertex
                canv_QA_hist->cd(1);
                TH1D *histVertex = getFromFile( file, "fHistVz" );

                //fill bin with nEvents for current run
                int nEventsInRun = histVertex->GetEntries();
                if ( nEventsInRun < 200000 )
                {
                    counter++;
                    continue;
                }
                fHistNumberOfEvents->SetBinContent( counter+1, nEventsInRun );


                histVertex->SetLineColor( color );//colors[i] );
                tuneHist1D(histVertex);
                histVertex->DrawNormalized( counter == 0 ? "" : "same" );
                histVertex->GetYaxis()->SetRangeUser(0, 1000000000 );



                //pt
                canv_QA_hist->cd(2);
                hist = getFromFile( file, "fHistPt" );
                hist->SetLineColor( color );
                tuneHist1D(hist);
                hist->Scale(1./nEventsInRun);
                hist->DrawCopy( counter == 0 ? "" : "same" );

                //eta
                canv_QA_hist->cd(3);
                hist = getFromFile( file, "fHistEta" );
                hist->SetLineColor( color );
                tuneHist1D(hist);
                hist->Scale(1./nEventsInRun);
                hist->DrawCopy( counter == 0 ? "" : "same" );

                leg->AddEntry(hist, runName, "l");

                //phi
                canv_QA_hist->cd(4);
                hist = getFromFile( file, "fHistPhi" );
                hist->SetLineColor( color );
                tuneHist1D(hist);
                hist->Scale(1./nEventsInRun);
                hist->DrawCopy( counter == 0 ? "" : "same" );



                file->Close();
                counter++;

//                if ( counter >= 5 )
//                    break;
            }
        }
    }
    nFiles = counter;


    //draw legend
    canv_QA_hist->cd(3);
    leg->Draw();

    //prepare canvas
    TCanvas *canv_nEvents_inRuns = new TCanvas( "canv_nEvents_inRuns","canv_nEvents_inRuns",100,50,700,600 );
    tuneCanvas(canv_nEvents_inRuns);

    tuneHist1D(fHistNumberOfEvents);
    fHistNumberOfEvents->GetXaxis()->SetLabelSize(0.03);
    fHistNumberOfEvents->DrawCopy();


    //    return filePtr;
}



void get_QA_plots_run_by_run_FROM_GRID_OUTPUT() // TString inputFileName = "MergedOutput.root")
{
    int nFiles = 0;
    list_files( nFiles, "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFplus_all"
                               , "AnalysisResults" );

    cout << "nFiles=" << nFiles << endl;

    return;
}

TH1D* getFromFile( TFile *file, TString what )
{
    file->ls();

    int listId = 0;//1;


    TList *listKeys = file->GetListOfKeys();
    cout << "going into list: " << listKeys->At(listId)->GetName() << endl;

    //return;
    file->cd( listKeys->At(listId)->GetName() );
    TList *listKeys2 = gDirectory->GetListOfKeys();
    TList *myTask = (TList*) gDirectory->Get( listKeys2->At(0)->GetName() );


    TString fname = file->GetName();
    fname.Remove(0,19);
    fname.Remove(6,5);
    cout << fname.Data() << endl;


        TH1D *histQA = (TH1D*)  myTask->FindObject( what );
    //    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistPt" );
    //    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistEta" );
//    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistPhi" );
    histQA->SetName( Form( "%s_%d", histQA->GetName(), fname.Data() ) );


    //histQA->DrawCopy();
    return histQA;


    TTree *t1 = (TTree*)  myTask->FindObject( "treeFB" );

    int nEvents = t1->GetEntries();
    cout <<"nEvents = " << nEvents << endl;


}
