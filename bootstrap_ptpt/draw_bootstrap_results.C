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


void draw_bootstrap_results()
{
    //    gROOT->ProcessLine(".L ../utils.C");

    TFile *f[9];

    f[0] = new TFile( "output_histos_graphs_20k.root" );
    f[1] = new TFile( "output_histos_graphs_50k.root" );
    f[2] = new TFile( "output_histos_graphs_100k.root" );
    f[3] = new TFile( "output_histos_graphs_200k.root" );
    f[4] = new TFile( "output_histos_graphs_500k.root" );
    f[5] = new TFile( "output_histos_graphs_1mln.root" );
    f[6] = new TFile( "output_histos_graphs_2mln.root" );
    f[7] = new TFile( "output_histos_graphs_5mln.root" );
    f[8] = new TFile( "output_histos_graphs_ALL_14mln.root" );


    const int nFiles = 9;


//    const int nCW = 1; //nCentrWidths
    const double cWidths[] = { 10, 5, 2.5, 1.0, 0.5 }; //width of the centrality bins
    const double cStep[] = { 5, 2.5,  2.5, 1.0, 1.0 }; //centrality bins step
    const int nCentrBins[] = { 17, 35, 36,  90, 90 }; //n centrality bins


    TGraphErrors *graphs[10][10][10];
    for ( int cW = 0; cW < 1; cW++ )
        for ( int etaW = 0; etaW < 1; etaW++ )
            for ( int fId = 0; fId < nFiles; fId++ )
                graphs[fId][cW][etaW] = (TGraphErrors*)f[fId]->Get( Form( "gr_BS_PtPt_cW%d_etaW%d_phiW0", cW, etaW ) );


    int etaId = 0;

    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",150,50,700,600 );
    tuneCanvas(canv_graphs);


    graphs[0][0][etaId]->SetTitle( ";centrality percentile;b_{corr}" );

    TLegend *leg = new TLegend(0.7,0.55,0.999,0.95);
    leg->SetHeader( "Number of events:" );
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    TString strNumEventsInFiles[nFiles] = { "20k", "50k", "100k", "200k", "500k", "1mln", "2mln", "5mln", "14mln" };

    //LHC10h
    for ( int fId = 0; fId < nFiles; fId++ )
    {
        shiftPointX(graphs[fId][0][etaId], fId*0.25);

        tuneGraphAxisLabels( graphs[fId][0][etaId] );
        graphs[fId][0][etaId]->GetXaxis()->CenterTitle();
        drawGraph( graphs[fId][0][etaId], 20, kOrange-5+fId, fId==0 ? "AP" : "P" );
        leg->AddEntry( graphs[fId][0][etaId], Form( "%s", strNumEventsInFiles[fId].Data() ), "p");
    }

    leg->Draw();

//return;
    //prepare graphs with dependence on ineff
//    const float trueEff = 0.8;
//    double ineffValues[] = { 0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 };
//    for ( int fId = 0; fId < nFiles; fId++ )
//        ineffValues[fId] *= trueEff;

    //const int nPoints = graphs[0][0][etaId]->GetN();
    int cW = 0;
    TGraphErrors *grErrorFromNevents[ nCentrBins[cW] ];

    //double x, y;
    for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
    {
        float cBinMin = cStep[cW] * cBin;
        float cBinMax = cWidths[cW] + cStep[cW] * cBin;

        grErrorFromNevents[cBin] = new TGraphErrors;

        for ( int fId = 0; fId < nFiles; fId++ )
        {
            TH2D *hist2D = (TH2D*)f[fId]->Get( Form("hist2D_PtPt_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0) );

            //take point and error from graph for this centrality class
            //graphs[fId][0][etaId]->GetPoint( i, x, y );
            double error = graphs[fId][0][etaId]->GetErrorY(cBin);

            //set to spec graph (point bootstrap error VS nEvents in centr. class)
            grErrorFromNevents[cBin]->SetPoint( fId, hist2D->GetEntries(), error );
        }
    }


    //draw dependence on ineff
    TCanvas *canv_error_from_nEvents = new TCanvas("canv_error_from_nEvents","canv_error_from_nEvents",250,150,700,600 );
    tuneCanvas(canv_error_from_nEvents);

//    TGraphErrors *grFrame = new TGraphErrors;
//    grFrame->SetPoint(0,0,0);
//    grFrame->SetPoint(1,0,0.15);
//    grFrame->SetPoint(2,100,0.15);
//    grFrame->SetPoint(3,100,0);
//    grFrame->SetMarkerColor(kWhite);
//    grFrame->SetTitle( ";efficiency, %;b_{corr}" );

//    tuneGraphAxisLabels( grFrame );
//    grFrame->Draw("AP");

    tuneGraphAxisLabels( grErrorFromNevents[0] );
    grErrorFromNevents[0]->SetTitle( ";number of events;bootstrap error of b_{corr}" );
    grErrorFromNevents[0]->GetXaxis()->CenterTitle();
    grErrorFromNevents[0]->GetYaxis()->CenterTitle();

//    TF1 *fFitFunc = new TF1("fitFunc","[0]/([1]+[2]*sqrt(x))",0,10000000);
//    fFitFunc->SetParameters( 0.1, 0., 1.);
    TF1 *fFitFunc = new TF1("fitFunc","[0]/sqrt(x)",0,10000000);
    fFitFunc->SetParameter( 0, 0.1 );

    int nPointsToDraw = 1;
    for ( int i = 0; i < nPointsToDraw; i++ )
    {
        drawGraph( grErrorFromNevents[i], 20, kOrange-5+i, i==0 ? "AP" : "P" );

        grErrorFromNevents[i]->Fit( fFitFunc );//,"Q");

        fFitFunc->SetLineColor(kOrange-5+i);
        fFitFunc->DrawCopy("same");
    }

    canv_error_from_nEvents->SetLogx();
    canv_error_from_nEvents->SetLogy();




}
