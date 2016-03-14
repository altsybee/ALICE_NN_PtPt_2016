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


void draw_graphs_INEFF_NEW_AMPLIFIED_HIJING()
{
    //    gROOT->ProcessLine(".L ../utils.C");

    const int nFiles = 9*2;
    TFile *f[nFiles];

//    f[0] = new TFile( "nEv_900k/output_histos_graphs_ineff_0.root" );
//    f[1] = new TFile( "nEv_900k/output_histos_graphs_ineff_0_5.root" );
//    f[2] = new TFile( "nEv_900k/output_histos_graphs_ineff_1_0.root" );
//    f[3] = new TFile( "nEv_900k/output_histos_graphs_ineff_1_5.root" );
//    f[4] = new TFile( "nEv_900k/output_histos_graphs_ineff_2_0.root" );
//    f[5] = new TFile( "nEv_900k/output_histos_graphs_ineff_2_5.root" );
//    f[6] = new TFile( "nEv_900k/output_histos_graphs_ineff_3_0.root" );
//    f[7] = new TFile( "nEv_900k/output_histos_graphs_ineff_3_5.root" );
//    f[8] = new TFile( "nEv_900k/output_histos_graphs_ineff_4_0.root" );

    for ( int fId = 0; fId < nFiles/2; fId++ )
//        f[fId] = new TFile( Form("nEv_3mln/output_histos_graphs_ineff_file%d.root", fId) );
    f[fId] = new TFile( Form("MFplus_nEv_3mln/output_histos_graphs_ineff_file%d.root", fId) );

    for ( int fId = nFiles/2; fId < nFiles; fId++ )
        f[fId] = new TFile( Form("MFminus_nEv_3mln/output_histos_graphs_ineff_file%d.root", fId-nFiles/2) );


    const int nCW = 1; //nCentrWidths
    const double cWidths[] = { 10, 5 /*, 2.5, 1.0, 0.5*/ }; //width of the centrality bins
    const double cStep[] = { 5, 2.5/*,  2.5, 1.0, 1.0*/ }; //centrality bins step
    const int nCentrBins[] = { 17, 35/*, 36,  90, 90*/ }; //n centrality bins

//    int colors[] = { kBlack, kBlack, kBlue, kRed, kMagenta };
//    int markers[] = { 20, 24, 5, 2, 3 };



    TGraphErrors *graphs[100][100][100];
    for ( int cW = 0; cW < 1; cW++ )
        for ( int etaW = 0; etaW < 1; etaW++ )
            for ( int fId = 0; fId < nFiles; fId++ )
                if ( fId < nFiles/2)
                    graphs[fId][cW][etaW] = (TGraphErrors*)f[fId]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
                else
                    graphs[fId][cW][etaW] = (TGraphErrors*)f[fId]->Get( Form( "gr_BS_PtPt_cW%d_etaW%d_phiW0", cW, etaW ) );

    int etaId = 0;

    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",150,50,700,600 );
    tuneCanvas(canv_graphs);

    tuneGraphAxisLabels(graphs[0][0][etaId]);

    graphs[0][0][etaId]->SetTitle( ";centrality percentile;b_{corr}" );

    //draw graphs
    for ( int fId = 0; fId < nFiles; fId++ )
        drawGraph( graphs[fId][0][etaId], 20, kBlack, fId==0 ? "APL" : "PL" );


//    for ( int fId = 0; fId < nFiles; fId++ )
//    {
//        cout << "f[fId]=" << f[fId] << endl;
//        f[fId]->ls();
//    }
//return;


    // ##### get mult in cBins in files
    TCanvas *canv_graphs_ineff_from_mult = new TCanvas("canv_graphs_ineff_from_mult","canv_graphs_ineff_from_mult",250,150,700,600 );
    tuneCanvas(canv_graphs_ineff_from_mult);


    TGraphErrors *grFrame = new TGraphErrors;
    grFrame->SetPoint(0,0,0);
    grFrame->SetPoint(1,0,0.2);
    grFrame->SetPoint(2,1000,0.2);
    grFrame->SetPoint(3,1000,0);
    grFrame->SetMarkerColor(kWhite);

    tuneGraphAxisLabels( grFrame );
    grFrame->SetTitle(";N_{ch};b_{corr}");
        grFrame->Draw("AP");


    TGraphErrors *grIneffFromMult[1000];
    int cW = 0;
    double x, y;
    for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
    {
        //mult hist: need <n> for cBin
        float cBinMin = cStep[cW] * cBin;
        float cBinMax = cWidths[cW] + cStep[cW] * cBin;

        grIneffFromMult[cBin] = new TGraphErrors;
        TString strGrName = Form( "grIneffFromMult_cBin%d", cBin );
        grIneffFromMult[cBin]->SetName( strGrName.Data() );
        tuneGraphAxisLabels( grIneffFromMult[cBin] );

        for ( int fId = 0; fId < nFiles; fId++ )
        {
            TString multNameF = Form( "hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
            TH1D *histMultF = (TH1D*)f[fId]->Get( multNameF );

            TString multNameB = Form( "hist1D_multDistrB_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
            TH1D *histMultB = (TH1D*)f[fId]->Get( multNameB );

            //histMult->DrawCopy( cBin == 0 ? "" : "same" );
            double meanNinCentrBin = histMultF->GetMean() + histMultB->GetMean();

            //get bCorr in this cBin
            graphs[fId][0][etaId]->GetPoint( cBin, x, y );
            double error = graphs[fId][0][etaId]->GetErrorY(cBin);

            //set point on ineff dep graph
            grIneffFromMult[cBin]->SetPoint( fId, meanNinCentrBin, y );
            grIneffFromMult[cBin]->SetPointError( fId, 0, error );
        }
        drawGraph( grIneffFromMult[cBin], 20, kOrange-9+cBin, "PL");//cBin==0 ? "APL" : "PL" );

        //fit func
        TF1 *fFitFunc = 0x0;
        if (1)
        {
            fFitFunc = new TF1("fitFunc","[0]*x/(1+[0]*x)",0,1200);
    //            TF1 *fFitFunc = new TF1("fitFunc","x/([0]+x)",0,1200);
    //            TF1 *fFitFunc = new TF1("fitFunc", "x/([0]+x)",0,1200);
            fFitFunc->SetParameter(0, 0.01);
        }
        if (0)
        {
            fFitFunc = new TF1("fFitFunc","[0]*(x+[1])/(1+[0]*(x+[1]))+[2]",0,1000);
            fFitFunc->SetParameter(0, 0.01);
            fFitFunc->SetParameter(1, 5);
            fFitFunc->SetParameter(2, 0.01);
        }
        if (0)
        {
            fFitFunc = new TF1("fFitFunc","[0]*(x+[1])/([2]+[0]*(x+[1]))",0,1000);
            fFitFunc->SetParameter(0, 0.01);
            fFitFunc->SetParameter(1, 5);
            fFitFunc->SetParameter(2, 1);
        }
        if (0)
        {
            fFitFunc = new TF1("fFitFunc","[0]*sqrt(x+[1])",0,1000);
            fFitFunc->SetParameter(0, 0.004);
            fFitFunc->SetParameter(1, 5);
        }


        for ( int i = 0; i < 2; i++ )
            grIneffFromMult[cBin]->Fit( fFitFunc,"Q");

        fFitFunc->SetLineColor(kOrange-9+cBin);
        //fFitFunc->SetLineStyle( fileId==0 ? 1 : 2 );
        fFitFunc->DrawCopy("same");

    }







    return;


}
