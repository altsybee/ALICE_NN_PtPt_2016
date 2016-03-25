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


void draw_bootstrap_results_LHC10h_ALL_EVENTS_V0M_cW10_5_2_5()
{
    //    gROOT->ProcessLine(".L ../utils.C");

    TFile *f[9];

    f[0] = new TFile( "DATA_ALL_EVENTS_with_BS_cW_10_5_2_5/output_histos_graphs_cW10_V0M.root" );
    f[1] = new TFile( "DATA_ALL_EVENTS_with_BS_cW_10_5_2_5/output_histos_graphs_cW5_V0M.root" );
    f[2] = new TFile( "DATA_ALL_EVENTS_with_BS_cW_10_5_2_5/output_histos_graphs_cW2_5_V0M.root" );


    const int nFiles = 3;


    //    const int nCW = 1; //nCentrWidths
//    const double cWidths[] = { 10, 5, 2.5, 1.0, 0.5 }; //width of the centrality bins
//    const double cStep[] = { 5, 2.5,  2.5, 1.0, 1.0 }; //centrality bins step
//    const int nCentrBins[] = { 17, 35, 36,  90, 90 }; //n centrality bins


    TGraphErrors *graphs[10][10][10];
    for ( int cW = 0; cW < 1; cW++ )
        for ( int etaW = 0; etaW < 1; etaW++ )
            for ( int fId = 0; fId < nFiles; fId++ )
                graphs[fId][cW][etaW] = (TGraphErrors*)f[fId]->Get( Form( "gr_BS_PtPt_cW%d_etaW%d_phiW0", cW, etaW ) );


    int etaId = 0;

    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",150,50,700,600 );
    tuneCanvas(canv_graphs);


    graphs[0][0][etaId]->SetTitle( ";centrality percentile;b_{corr}" );

    TLegend *leg = new TLegend(0.6,0.75,0.995,0.95);
    //    leg->SetHeader( "Number of events:" );
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    TString strNumEventsInFiles[] = { "class width 10%", "class width 5%", "class width 2.5%" };

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


    //TLatex *tex = new TLatex(0.2,0.9, "#delta#eta=0.4");
    TLatex *tex = new TLatex(0.25,0.28, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex);

    tex = new TLatex(0.75,0.7, "(V0M)");
    drawTex(tex);


    return;
}
