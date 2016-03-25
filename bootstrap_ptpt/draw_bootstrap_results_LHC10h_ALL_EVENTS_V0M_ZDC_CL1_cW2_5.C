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

#include "../utils.C"


void draw_bootstrap_results_LHC10h_ALL_EVENTS_V0M_ZDC_CL1_cW2_5()
{
    //    gROOT->ProcessLine(".L ../utils.C");

    TFile *f[9];

//    f[0] = new TFile( "DATA_ALL_EVENTS_with_BS_cW_10_5_2_5/output_histos_graphs_cW10_V0M_noPileUp_rejection.root" );
    f[0] = new TFile( "DATA_ALL_EVENTS_with_BS_cW_10_5_2_5/output_histos_graphs_cW2_5_V0M.root" );
    f[1] = new TFile( "DATA_ALL_EVENTS_with_BS_cW_10_5_2_5/output_histos_graphs_cW2_5_ZDC.root" );
    f[2] = new TFile( "DATA_ALL_EVENTS_with_BS_cW_10_5_2_5/output_histos_graphs_cW2_5_CL1.root" );

//    const int nFiles = 1;//3;
    const int nFiles = 3;


    const int nCW = 1; //nCentrWidths
    const double cWidths[nCW] = { 2.5 }; //width of the centrality bins
    const double cStep[nCW] = { 2.5 }; //centrality bins step
    const int nCentrBins[nCW] = { 34 };//18 };//34 }; //n centrality bins

//    const int nCW = 1; //nCentrWidths
//    const double cWidths[nCW] = { 10 }; //width of the centrality bins
//    const double cStep[nCW] = { 10 }; //centrality bins step
//    const int nCentrBins[nCW] = { 8 }; //n centrality bins


    const int nEtaWins = 1;

    TGraphErrors *graphs[10][10][10];
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int etaW = 0; etaW < 1; etaW++ )
            for ( int fId = 0; fId < nFiles; fId++ )
                graphs[fId][cW][etaW] = (TGraphErrors*)f[fId]->Get( Form( "gr_BS_PtPt_cW%d_etaW%d_phiW0", cW, etaW ) );


    int etaId = 0;

    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",150,50,700,600 );
    tuneCanvas(canv_graphs);


    graphs[0][0][etaId]->SetTitle( ";centrality percentile;b_{corr}" );

    TLegend *leg = new TLegend(0.7,0.75,0.95,0.95);
    //    leg->SetHeader( "Number of events:" );
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    TString strNumEventsInFiles[] = { "V0M", "ZDC", "CL1" };

    int colors[] = { kBlue, kMagenta, kRed };
    int markers[] = { 21, 24, 20 };

    //LHC10h
    for ( int fId = 0; fId < nFiles; fId++ )
    {
        TGraphErrors *gr = graphs[fId][0][etaId];

        tuneGraphAxisLabels( gr );
        gr->GetXaxis()->CenterTitle();
        gr->GetYaxis()->CenterTitle();
        drawGraph( gr, markers[fId], colors[fId] /*kOrange-5+fId*/, fId==0 ? "AP" : "P", 0.8 );
        leg->AddEntry( gr, Form( "%s", strNumEventsInFiles[fId].Data() ), "p");
    }

    leg->Draw();

    TLatex *tex = new TLatex(0.25,0.28, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex);

    // ####### for QA: graph with bootstrap abs error as func of centrality bin
    TCanvas *canv_BS_abs_error = new TCanvas("canv_BS_abs_error","canv_BS_abs_error",450,150,700,600 );
    tuneCanvas(canv_BS_abs_error);

    TGraphErrors *graphsHighStat_BS_abs_error[nFiles][nCW][nEtaWins];
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
        {
            for ( int fId = 0; fId < nFiles; fId++ )
            {
                graphsHighStat_BS_abs_error[fId][cW][etaW] = new TGraphErrors;
                graphsHighStat_BS_abs_error[fId][cW][etaW]->SetName( Form( "grHighStat_BS_abs_error_cW%d_etaW%d", cW, etaW ) );
                graphsHighStat_BS_abs_error[fId][cW][etaW]->SetTitle( ";centrality percentile;abs stat. error of b_{corr}" );

                for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
                {
                    float cBinMin = cStep[cW] * cBin;
                    float cBinMax = cWidths[cW] + cStep[cW] * cBin;

                    float cBinCentre = cBinMin + (cBinMax - cBinMin)/2;

                    double pointError = graphs[fId][cW][etaId]->GetErrorY( cBin );
                    if ( fId == 1 /*ZDC*/ && cBinCentre > 44 )
                        continue;

                    graphsHighStat_BS_abs_error[fId][cW][etaId]->SetPoint( cBin, cBinCentre, pointError );
                    graphsHighStat_BS_abs_error[fId][cW][etaId]->SetPointError( cBin, 0, 0 );

                    double x,y;
                    graphs[fId][cW][etaId]->GetPoint( cBin, x, y );
                    cout << "QA: pointError/bCorr = " << pointError/y << endl;

                }
                tuneGraphAxisLabels( graphsHighStat_BS_abs_error[fId][cW][etaId] );
                drawGraph( graphsHighStat_BS_abs_error[fId][cW][etaId], 20, colors[fId], fId == 0 ? "APL" : "PL" );
            }
        }
    leg->Draw();


    return;
}
