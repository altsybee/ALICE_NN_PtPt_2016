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

#include "../utils.C"



void draw_graphs()
{
    //        gROOT->ProcessLine(".L ../utils.C");
    //        gROOT->ProcessLine(".L AliLRCFit.cxx");
    TFile *f[8];

    //    f[0] = new TFile( "output_classesByV0M_LHC10h.root" );
    //    f[0] = new TFile( "output_classesByV0M_LHC10h_c10_5_1.root" );
    f[0] = new TFile( "output_classesByV0M_LHC10h_c10_5_25_1_05.root" );
//    f[0] = new TFile( "output_classesByZDCZEM_LHC10h_c10_5_25_1_05.root" );



    //    f[0]->ls();
    //    const int nCW = 2; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35 }; //n centrality bins
    //    const int nCW = 3; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 1.0 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5, 1.0 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35, 90 }; //n centrality bins

//    const int nCW = 4; //nCentrWidths
//    const double cWidths[nCW] = { 10, 5, 1.0, 0.5 }; //width of the centrality bins
//    const double cStep[nCW] = { 5, 2.5, 1.0, 1.0 }; //centrality bins step
//    const int nCentrBins[nCW] = { 17, 35, 90, 90 }; //n centrality bins

    const int nCW = 4; //nCentrWidths
    const double cWidths[] = { 10, 5, 2.5, 1.0, 0.5 }; //width of the centrality bins
    const double cStep[] = { 5, 2.5,  2.5, 1.0, 1.0 }; //centrality bins step
    const int nCentrBins[] = { 17, 35, 36,  90, 90 }; //n centrality bins


    int colors[] = { kBlack, kBlack, kBlue, kRed, kMagenta };
    int markers[] = { 20, 24, 5, 2, 3 };


    const int nEtaWins = 3;

    const int kCorrType = 1; //0-NN, 1-PtPt, 2-PtN

    TGraphErrors *graphCorr[nCW][nEtaWins];

    TGraphErrors *grNN[nCW][nEtaWins];
    TGraphErrors *grPtPt[nCW][nEtaWins];
    TGraphErrors *grPtN[nCW][nEtaWins];

    for ( int cW = 0; cW < nCW; cW++ )
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
        {
            grNN[cW][etaW] = (TGraphErrors*)f[0]->Get( Form( "grNN_c%d_eta%d", cW, etaW ) );
            grPtPt[cW][etaW] = (TGraphErrors*)f[0]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
            grPtN[cW][etaW] = (TGraphErrors*)f[0]->Get( Form( "grPtN_c%d_eta%d", cW, etaW ) );


//            cout << grNN[cW][etaW] << endl;

            if ( kCorrType == 0) graphCorr[cW][etaW] = grNN[cW][etaW];
            else if ( kCorrType == 1) graphCorr[cW][etaW] = grPtPt[cW][etaW];
            else if ( kCorrType == 2) graphCorr[cW][etaW] = grPtN[cW][etaW];
        }
    const int etaId = 0;

    // ########## Graphs
    TCanvas *canv_graphCorr = new TCanvas("canv_graphCorr","canv_graphCorr",20,50,700,600 );
    tuneCanvas(canv_graphCorr);
    graphCorr[0][etaId]->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels(graphCorr[0][etaId]);


    //remove some points by hand:
//    graphCorr[0][etaId]->RemovePoint( nCentrBins[0]-1 );
//    if (nCW>1) graphCorr[1][etaId]->RemovePoint( nCentrBins[1]-1 );

    //draw for several centralities:
    for ( int cW = 0; cW < nCW; cW++ )
        drawGraph(graphCorr[cW][etaId], markers[cW], colors[cW], cW == 0 ? "AP" : "P");

    graphCorr[0][etaId]->SetMinimum( 0 );

    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(graphCorr[cW][etaId], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();

//    graphCorr[0][etaId]->GetYaxis()->SetRangeUser( 0, 0.92 );

//    TLatex *tex = new TLatex(0.4,0.89, "#eta_{gap}=0.8, #delta#eta=0.4");
    TLatex *tex = new TLatex(0.5,0.3, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);

//    tex = new TLatex(0.4,0.89, "NN");
//    drawTex(tex, 0.045);


    //    TString strPostfix;

    //    if (flag_V0M_ZDC==0)
    //        strPostfix = Form("V0M.eps");
    //    else
    //        strPostfix = Form("ZDCZEM.eps");

    //    canv_graphCorr->SaveAs( Form("NN_%s", strPostfix.Data() ) );




    return;


    //    gROOT->ProcessLine( ".q");
}
