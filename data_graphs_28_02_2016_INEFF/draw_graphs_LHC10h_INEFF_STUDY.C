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


void draw_graphs_LHC10h_INEFF_STUDY()
{
    //    gROOT->ProcessLine(".L ../utils.C");

    TFile *f[8];

    f[0] = new TFile( "output_histos_graphs_INEFF_0.root" );
    f[1] = new TFile( "output_histos_graphs_INEFF_05.root" );
    f[2] = new TFile( "output_histos_graphs_INEFF_10.root" );
    f[3] = new TFile( "output_histos_graphs_INEFF_20.root" );
    f[4] = new TFile( "output_histos_graphs_INEFF_30.root" );
    f[5] = new TFile( "output_histos_graphs_INEFF_40.root" );
    f[6] = new TFile( "output_histos_graphs_INEFF_50.root" );


    const int nFiles = 7;

    TGraphErrors *graphs[10][10][10];
    for ( int cW = 0; cW < 2; cW++ )
        for ( int etaW = 0; etaW < 3; etaW++ )
            for ( int fId = 0; fId < nFiles; fId++ )
                graphs[fId][cW][etaW] = (TGraphErrors*)f[fId]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );


    int etaId = 2;

    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",150,50,700,600 );
    tuneCanvas(canv_graphs);

    tuneGraphAxisLabels(graphs[0][0][etaId]);

    graphs[0][0][etaId]->SetTitle( ";centrality percentile;b_{corr}" );

    //LHC10h MF plus
    for ( int fId = 0; fId < nFiles; fId++ )
        drawGraph( graphs[fId][0][etaId], 20, kBlack, fId==0 ? "APL" : "PL" );


    //prepare graphs with dependence on ineff
    const float trueEff = 0.8;
    double ineffValues[] = { 0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 };
//    for ( int fId = 0; fId < nFiles; fId++ )
//        ineffValues[fId] *= trueEff;

    const int nPoints = graphs[0][0][etaId]->GetN();
    TGraphErrors *grDepOnIneff[nPoints];

    double x, y;
    for ( int i = 0; i < nPoints; i++ )
    {
        grDepOnIneff[i] = new TGraphErrors;
        for ( int fId = 0; fId < nFiles; fId++ )
        {
            graphs[fId][0][etaId]->GetPoint( i, x, y );
            grDepOnIneff[i]->SetPoint( fId, 100*trueEff*(1-ineffValues[fId]), y );
        }
    }


    //draw dependence on ineff
    TCanvas *canv_dep_on_ineff = new TCanvas("canv_dep_on_ineff","canv_dep_on_ineff",250,150,700,600 );
    tuneCanvas(canv_dep_on_ineff);

    TGraphErrors *grFrame = new TGraphErrors;
    grFrame->SetPoint(0,0,0);
    grFrame->SetPoint(1,0,0.15);
    grFrame->SetPoint(2,100,0.15);
    grFrame->SetPoint(3,100,0);
    grFrame->SetMarkerColor(kWhite);
    grFrame->SetTitle( ";efficiency, %;b_{corr}" );

    tuneGraphAxisLabels( grFrame );
    grFrame->Draw("AP");

    TF1 *fFitFunc = new TF1("fitFunc","[0]*x/(1+[0]*x)",0,110);
    fFitFunc->SetParameter(0,0.1);

    int nPointsToDraw = 12;
    for ( int i = 0; i < nPointsToDraw; i++ )
    {
        drawGraph( grDepOnIneff[i], 20, kOrange-5+i, "P" ); //i==0 ? "APL" : "PL" );

        grDepOnIneff[i]->Fit( fFitFunc );//,"Q");

        fFitFunc->SetLineColor(kOrange-5+i);
        fFitFunc->DrawCopy("same");
    }







    //    drawGraph( gr[0][1][etaId], 24, kBlack, "PL" );

    //    //LHC10h MF minus
    //    drawGraph( gr[1][0][etaId], 20, kBlue, "PL" );
    //    drawGraph( gr[1][1][etaId], 24, kBlue, "PL" );

    //    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    //    leg->SetFillColor(kWhite);
    //    leg->SetFillStyle(0);
    //    leg->SetBorderSize(0);

    //    leg->AddEntry( gr[0][0][etaId], "LHC10h MF++, class 10", "p");
    //    leg->AddEntry( gr[0][1][etaId], "LHC10h MF++, class 5", "p");
    //    leg->AddEntry( gr[1][0][etaId], "LHC10h MF--, class 10", "p");
    //    leg->AddEntry( gr[1][1][etaId], "LHC10h MF--, class 5", "p");

    //    leg->Draw();


    //    gROOT->ProcessLine( ".q");
}
