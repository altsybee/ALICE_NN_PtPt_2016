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

void draw_a_graphs()
{
    TFile *f[8];

    f[0] = new TFile( "graph_A_LHC10h_eta0.root" );
    f[1] = new TFile( "graph_A_LHC15o_eta0.root" );

    f[2] = new TFile( "graph_A_LHC10h_eta2.root" );
    f[3] = new TFile( "graph_A_LHC15o_eta2.root" );


    TGraphErrors *gr[10];
    gr[0] = (TGraphErrors*)f[0]->Get("grA");
    gr[1] = (TGraphErrors*)f[1]->Get("grA");

    gr[2] = (TGraphErrors*)f[2]->Get("grA");
    gr[3] = (TGraphErrors*)f[3]->Get("grA");


    int colors[] = { kBlack, kBlack, kBlue, kRed, kMagenta };
    int markers[] = { 20, 24, 5, 2, 3 };

    TCanvas *canv_par_a = new TCanvas("canv_par_a","canv_par_a",200,300,700,600 );
    tuneCanvas(canv_par_a);
    tuneGraphAxisLabels( gr[0] );

    gr[0]->SetTitle( ";centrality percentile;coeff a" );

    for ( int i = 0; i < 2; i++ )
    {
        drawGraph( gr[i], markers[i], colors[i], i == 0 ? "APL" : "PL");
    }

    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry( gr[0], "LHC10h #eta_{gap}=0.8", "p");
    leg->AddEntry( gr[1], "LHC15o #eta_{gap}=0.8", "p");
//    leg->AddEntry( gr[2], "LHC10h #eta_{gap}=0", "p");
//    leg->AddEntry( gr[3], "LHC15o #eta_{gap}=0", "p");

    leg->Draw();

}
