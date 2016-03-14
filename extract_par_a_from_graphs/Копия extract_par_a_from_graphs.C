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



void extract_par_a_from_graphs()
{
    TFile *f[8];

//    f[0] = new TFile( "../data_graphs_23_02_2016/new_eW_binning_27_02_2016/output_histos_graphs_LHC10h_eW04_MFplus.root" );
    f[0] = new TFile( "../data_graphs_23_02_2016/output_classesByV0M_LHC15o_fieldCOMBINED_c10_5_CUT_OUTLIERS.root" );

    const int nCW = 2; //nCentrWidths
    const double cWidths[] = { 10, 5 /*, 2.5, 1.0, 0.5*/ }; //width of the centrality bins
    const double cStep[] = { 5, 2.5/*,  2.5, 1.0, 1.0*/ }; //centrality bins step
    const int nCentrBins[] = { 17, 35/*, 36,  90, 90*/ }; //n centrality bins

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


    // ###### extract a:
    TCanvas *canv_histMults = new TCanvas("canv_histMults","canv_histMults",20,50,700,600 );
    tuneCanvas(canv_histMults);
    canv_histMults->SetGrid(1,1);

    TGraphErrors *grFrame = new TGraphErrors;
    grFrame->SetPoint(0,0,0);
    grFrame->SetPoint(1,0,1);
    grFrame->SetPoint(2,1000,1);
    grFrame->SetPoint(3,1000,0);
    grFrame->SetMarkerColor(kWhite);

    tuneGraphAxisLabels( grFrame );
    grFrame->SetTitle(";N_{ch};b_{corr}");
    grFrame->Draw("AP");


    const int cW = 1;
    TGraphErrors *gr = graphCorr[cW][etaId];

    TGraphErrors *grA = new TGraphErrors;
    grA->SetTitle(";centrality percentile;coeff a");
    grA->SetName("grA");
    double x, y;
    for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
    {
        //mult hist: need <n> for cBin
        float cBinMin = cStep[cW] * cBin;
        float cBinMax = cWidths[cW] + cStep[cW] * cBin;

        TString multNameF = Form( "hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
        TH1D *histMultF = (TH1D*)f[0]->Get( multNameF );

        TString multNameB = Form( "hist1D_multDistrB_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
        TH1D *histMultB = (TH1D*)f[0]->Get( multNameB );

        //histMult->DrawCopy( cBin == 0 ? "" : "same" );
        double meanNinCentrBin = histMultF->GetMean() + histMultB->GetMean();


        //get bCorr in this cBin
        gr->GetPoint( cBin, x, y );
        TGraphErrors *grForFit = new TGraphErrors;

        grForFit->SetPoint(0, 0, 0);
        grForFit->SetPoint(1, meanNinCentrBin, y );

        grForFit->SetMarkerStyle(30);
        grForFit->DrawClone( cBin==0? "P" : "P" );

        //fit func
        TF1 *fFitFunc = new TF1("fitFunc","[0]*x/(1+[0]*x)",0,1200);
        fFitFunc->SetParameter(0,0.01);


        for ( int i = 0; i < 2; i++ )
            grForFit->Fit( fFitFunc );//,"Q");

        fFitFunc->SetLineColor(kOrange-5+cBin);
        fFitFunc->DrawCopy("same");

        //get parameter a
        double a = fFitFunc->GetParameter(0);
        grA->SetPoint(cBin, x, a);

    }

    TCanvas *canv_par_a = new TCanvas("canv_par_a","canv_par_a",200,300,700,600 );
    tuneCanvas(canv_par_a);
    tuneGraphAxisLabels( grA );
    grA->SetMarkerStyle(20);
    grA->Draw("APL");

//    grA->SaveAs( "graph_A_LHC10h_eta0.root" );
    grA->SaveAs( "graph_A_LHC15o_eta0.root" );


    return;
}
