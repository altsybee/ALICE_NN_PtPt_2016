#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>


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

#include "SupplementaryClasses.h"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;

#include "utils.C"

void get_QA_plots_from_MergedOutputs() // TString inputFileName = "MergedOutput.root")
{
    int colors[] = { kBlack, kBlue, kRed, kMagenta, kCyan };

    TCanvas *canv_QA_hist = new TCanvas("canv_QA_hist","canv_QA_hist",20,50,700,600 );
    tuneCanvas(canv_QA_hist);


    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    TString strLegend[] = { "LHC10h", "LHC11h_MF--", "LHC11h_MF+++"
                            , "LHC15o_MF--", "LHC15o_MF++" };
    for ( int fId = 0; fId < 5; fId++ )
    {
        TH1D *hist = getHistByFileId( fId );
        hist->SetLineColor( colors[fId] );
        tuneHist1D(hist);
        hist->DrawNormalized( fId == 0 ? "" : "same" );

        leg->AddEntry(hist, strLegend[fId], "l");
    }
    leg->Draw();

}

TH1D* getHistByFileId(int fileId)
{
    TFile *myFile;

    //LHC10h
    //    if (fileId==0) fileId = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_09_PbPb_Data_2_76TeV_LHC10h_AOD0160_3etaWins_phi1_pt02_20_FB_TREE_blocks123/block1/AnalysisResults.139465.root" );
    if (fileId==0) myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_09_PbPb_Data_2_76TeV_LHC10h_AOD0160_3etaWins_phi1_pt02_20_FB_TREE_blocks123/MergedOutput.root" );

    //LHC11h
    if (fileId==1) myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoMinus/MergedOutput.root" );
    if (fileId==2) myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoPlus/MergedOutput.root" );

    //LHC11h Femto Plus Minus merged
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoPlus/MergedResults/MergedOutput_Femto_Plus_Minus.root" );

    //LHC15o MF ++ --
    if (fileId==3) myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_22_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldMM/MergedOutput.root" );
    if (fileId==4) myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_21_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldPP/MergedOutput.root" );

    //LHC15o MF merged
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_22_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldMM/MergedResults/MergedOutput_PP_and_MM.root" );

    if (!myFile)
    {
        cout << "No input file!" << endl;
        return;
    }
    myFile->ls();

    // !!! important for listId, also for branch V0M below!
    bool isAnalysing502 = fileId>2;

    // !!! listId can be different for 2.76 and 5.02!
    int listId = 0;//1;
    if ( isAnalysing502 )
        listId = 1;


    TList *listKeys = myFile->GetListOfKeys();
    cout << "going into list: " << listKeys->At(listId)->GetName() << endl;

    //return;
    myFile->cd( listKeys->At(listId)->GetName() );
    TList *listKeys2 = gDirectory->GetListOfKeys();
    TList *myTask = (TList*) gDirectory->Get( listKeys2->At(0)->GetName() );



//    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistVz" );
//    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistPt" );
//    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistEta" );
    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistPhi" );
    histQA->SetName( Form( "%s_%d", histQA->GetName(), fileId ) );


    //histQA->DrawCopy();
    return histQA;


    TTree *t1 = (TTree*)  myTask->FindObject( "treeFB" );

    int nEvents = t1->GetEntries();
    cout <<"nEvents = " << nEvents << endl;


}
