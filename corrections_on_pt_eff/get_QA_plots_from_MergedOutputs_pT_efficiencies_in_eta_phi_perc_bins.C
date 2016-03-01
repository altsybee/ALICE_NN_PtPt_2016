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

#include "../SupplementaryClasses.h"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;

#include "../utils.C"


TH1D* get_QA_plots_from_MergedOutputs_pT_efficiencies_in_eta_phi_perc_bins() //int fileId)
{
    TFile *myFile;

    //HIJING LHC11a10a_bis 900k events AOD162
    myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_28_PbPb_MCAOD_LHC11a10a_bis_AOD162_TEST_Efficiency_kine_vs_reco_try5/MergedOutput.root" );

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

    TCanvas *canv_QA_diff = new TCanvas("canv_QA_diff","canv_QA_diff",350,50,700,600 );
    tuneCanvas(canv_QA_diff);

    const int nEtaWins = 8;
    const int nPhiWins = 1;
    const int nPtBins = 1;

    TH2D *fHist2D_ptDistrInEtaPhiVsCentr_winF_kine[nPtBins][nEtaWins][nPhiWins];
    TH2D *fHist2D_ptDistrInEtaPhiVsCentr_winB_kine[nPtBins][nEtaWins][nPhiWins];
    TH2D *fHist2D_ptDistrInEtaPhiVsCentr_winF_reco[nPtBins][nEtaWins][nPhiWins];
    TH2D *fHist2D_ptDistrInEtaPhiVsCentr_winB_reco[nPtBins][nEtaWins][nPhiWins];

    for ( int etaW = 0; etaW < nEtaWins; etaW++ )
        for ( int phiW = 0; phiW < nPhiWins; phiW++ )
            for ( int ptW = 0; ptW < nPtBins; ptW++ )
            {
                TString namePostfix = Form("eta_%d_phi%d_pt_%d"
                                           , etaW, phiW, ptW );
                // ### for kine:
                //F:
                fHist2D_ptDistrInEtaPhiVsCentr_winF_kine[ptW][etaW][phiW] = (TH2D*)  myTask->FindObject( Form("fHist2D_ptDistrInEtaPhiVsCentr_winF_kine_%s", namePostfix.Data() ) );

                //B:
                fHist2D_ptDistrInEtaPhiVsCentr_winB_kine[ptW][etaW][phiW] = (TH2D*)  myTask->FindObject( Form("fHist2D_ptDistrInEtaPhiVsCentr_winB_kine_%s", namePostfix.Data() ) );

                // ### for reco:
                //F:
                fHist2D_ptDistrInEtaPhiVsCentr_winF_reco[ptW][etaW][phiW] = (TH2D*)  myTask->FindObject( Form("fHist2D_ptDistrInEtaPhiVsCentr_winF_reco_%s", namePostfix.Data() ) );

                //B:
                fHist2D_ptDistrInEtaPhiVsCentr_winB_reco[ptW][etaW][phiW] = (TH2D*)  myTask->FindObject( Form("fHist2D_ptDistrInEtaPhiVsCentr_winB_reco_%s", namePostfix.Data() ) );


                // ### prepare slices:
                TH1D *projBin_kine_1 = fHist2D_ptDistrInEtaPhiVsCentr_winB_kine[ptW][etaW][phiW]->ProjectionY( Form("%s_proj%d", fHist2D_ptDistrInEtaPhiVsCentr_winF_kine[ptW][etaW][phiW]->GetName(), 0 ),
                        1, 2);

//                projBin_kine_1->DrawCopy();

                TH1D *projBin_reco_1 = fHist2D_ptDistrInEtaPhiVsCentr_winB_reco[ptW][etaW][phiW]->ProjectionY( Form("%s_proj%d", fHist2D_ptDistrInEtaPhiVsCentr_winF_reco[ptW][etaW][phiW]->GetName(), 0 ),
                        1, 2);
                projBin_reco_1->SetLineColor(kOrange-5+etaW);
//                projBin_reco_1->DrawCopy( "same" );

                // ### eff:
                canv_QA_eff->cd();
                TH1D *histEff = projBin_reco_1->Clone( Form( "%s_eff_clone", projBin_reco_1->GetName() ) );
                histEff->Divide( projBin_kine_1 );
                histEff->DrawCopy( (etaW==0 && phiW==0 && ptW==0) ? "" : "same" );

                // ### diff
                canv_QA_diff->cd();
                TH1D *histDiff_reco = projBin_reco_1->Clone( Form( "%s_diff_clone", projBin_reco_1->GetName() ) );
                TH1D *histDiff_kine = projBin_kine_1->Clone( Form( "%s_diff_clone", projBin_kine_1->GetName() ) );
                histDiff_kine->Add( histDiff_reco, -1 );
//                histDiff_kine->Divide( projBin_kine_1 );//, -1 );
                histDiff_kine->DrawCopy( (etaW==0 && phiW==0 && ptW==0) ? "" : "same" );
            }



    //    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistVz" );
    //    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistPt" );
    //    TH1D *histQA = (TH1D*)  myTask->FindObject( "fHistEta" );
    //TH1D *histQA = (TH1D*)  myTask->FindObject( "fHist2D_ptDistrInEtaPhiVsCentr_winF_kine_" );
    //histQA->SetName( Form( "%s_%d", histQA->GetName(), fileId ) );


    //histQA->DrawCopy();


    return;

    int colors[] = { kBlack, kBlue, kRed, kMagenta, kCyan };



    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    TString strLegend[] = { "LHC10h", "LHC11h_MF--", "LHC11h_MF+++"
                            , "LHC15o_MF--", "LHC15o_MF++" };
    for ( int fId = 0; fId < 1; fId++ )
    {
        //        TH1D *hist = getHistByFileId( fId );
        //        hist->SetLineColor( colors[fId] );
        //        tuneHist1D(hist);
        //        hist->DrawNormalized( fId == 0 ? "" : "same" );

        leg->AddEntry(hist, strLegend[fId], "l");
    }
    leg->Draw();







    TTree *t1 = (TTree*)  myTask->FindObject( "treeFB" );

    int nEvents = t1->GetEntries();
    cout <<"nEvents = " << nEvents << endl;


}
