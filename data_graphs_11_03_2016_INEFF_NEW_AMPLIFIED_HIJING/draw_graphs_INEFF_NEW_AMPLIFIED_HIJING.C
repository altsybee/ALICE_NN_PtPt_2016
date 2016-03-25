#include <TFile.h>


#include "TObject.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1.h"
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
    gSystem->cd( "HIJING_RECO_KINE" );
    gROOT->ProcessLine(".L extract_correction_factors_from_HIJING_RECO_KINE.C");
    double *corrFactors = extract_correction_factors_from_HIJING_RECO_KINE();
    gSystem->cd( ".." );

    //return;

    const int nFiles = 8;//9;//9*2;
    TFile *f[nFiles];

    //INEFF data
    for ( int fId = 0; fId < nFiles; fId++ )
        f[fId] = new TFile( Form("MF_all_combined_6mln/output_histos_graphs_ineff_file%d.root", fId) );


    TFile *fileHighStat = new TFile( "../bootstrap_ptpt/DATA_ALL_EVENTS_with_BS_cW_10_5_2_5/output_histos_graphs_cW10_V0M.root" );  //output_histos_graphs_ALL_14mln.root" );



    const int nCW = 1; //nCentrWidths
    const double cWidths[] = { 10, 5 /*, 2.5, 1.0, 0.5*/ }; //width of the centrality bins
    const double cStep[] = { 5, 2.5/*,  2.5, 1.0, 1.0*/ }; //centrality bins step
    const int nCentrBins[] = { 15, 35/*, 36,  90, 90*/ }; //n centrality bins

    const int nEtaWins = 1;
    //    int colors[] = { kBlack, kBlack, kBlue, kRed, kMagenta };
    //    int markers[] = { 20, 24, 5, 2, 3 };

    double ineffValues[] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 };

    TGraphErrors *graphs_ineff[nFiles][nCW][nEtaWins]; //ineff
    TGraphErrors *graphsHighStat_cBin[nCW][nEtaWins]; //main high stat graphs to be corrected
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
        {
            //high stat graphs
            graphsHighStat_cBin[cW][etaW] = (TGraphErrors*)fileHighStat->Get( Form( "gr_BS_PtPt_cW%d_etaW%d_phiW0", cW, etaW ) );
            //ineff graphs
            for ( int fId = 0; fId < nFiles; fId++ )
            {
                if ( fId == 0 )// nFiles/2)
                    graphs_ineff[fId][cW][etaW] = (TGraphErrors*)f[fId]->Get( Form( "gr_BS_PtPt_cW%d_etaW%d_phiW0", cW, etaW ) );
                else
                {
                    graphs_ineff[fId][cW][etaW] = (TGraphErrors*)f[fId]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
                    // !!! take error from fId=0 !!! (=without amplified ineff!)
                    for ( int p = 0; p < graphs_ineff[fId][cW][etaW]->GetN(); p++ )
                    {
                        double error = graphs_ineff[0][cW][etaW]->GetErrorY(p);
                        graphs_ineff[fId][cW][etaW]->SetPointError( p, 0, error );
                    }
                }
            }
        }

    int etaId = 0;



    //draw ineff graphs (for QA)
    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",50,50,700,600 );
    tuneCanvas(canv_graphs);

    tuneGraphAxisLabels(graphs_ineff[0][0][etaId]);
    graphs_ineff[0][0][etaId]->SetTitle( ";centrality percentile;b_{corr}" );

    for ( int fId = 0; fId < nFiles; fId++ )
        drawGraph( graphs_ineff[fId][0][etaId], 20, kBlack, fId==0 ? "APL" : "PL" );



    // ##### get mult in cBins in files
    TCanvas *canv_graphs_ineff_from_mult = new TCanvas("canv_graphs_ineff_from_mult","canv_graphs_ineff_from_mult",250,150,700,600 );
    tuneCanvas(canv_graphs_ineff_from_mult);


    TGraphErrors *grFrame = new TGraphErrors;
    grFrame->SetPoint(0,0,0);
    grFrame->SetPoint(1,0,0.2);
    grFrame->SetPoint(2,1000,0.2);
    grFrame->SetPoint(3,1000,0);
    //    grFrame->SetPoint(2,1,0.2);
    //    grFrame->SetPoint(3,1,0);
    grFrame->SetMarkerColor(kWhite);

    tuneGraphAxisLabels( grFrame );
    grFrame->SetTitle(";N_{ch};b_{corr}");
    grFrame->Draw("AP");


    //prepare graphs with mult on x-axis
    TGraphErrors *graphsHighStat_mult[nCW][nEtaWins];

    // for QA: graph with bootstrap abs error as func of centrality bin
    TGraphErrors *graphsHighStat_BS_abs_error[nCW][nEtaWins];

    //prepare graphs with corrected points
    TGraphErrors *graphsCorrected_mult[nCW][nEtaWins];
    TGraphErrors *graphsCorrected_cBin[nCW][nEtaWins];
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
        {
            graphsHighStat_mult[cW][etaW] = (TGraphErrors*)graphsHighStat_cBin[cW][etaW]->Clone( Form( "grHighStat_multOnX_cW%d_etaW%d", cW, etaW ));

            //graphsHighStat_BS_abs_error[cW][etaW] = (TGraphErrors*)graphsHighStat_cBin[cW][etaW]->Clone( Form( "grHighStat_BS_abs_error_cW%d_etaW%d", cW, etaW ));
            graphsHighStat_BS_abs_error[cW][etaW] = new TGraphErrors;
            graphsHighStat_BS_abs_error[cW][etaW]->SetName( Form( "grHighStat_BS_abs_error_cW%d_etaW%d", cW, etaW ) );

            graphsCorrected_mult[cW][etaW] = new TGraphErrors;
            graphsCorrected_mult[cW][etaW]->SetName( Form( "corrected_points_cW%d_etaW%d", cW, etaW ) );

            graphsCorrected_cBin[cW][etaW] = new TGraphErrors;
            graphsCorrected_cBin[cW][etaW]->SetName( Form( "corrected_points_xAxis_cBin_cW%d_etaW%d", cW, etaW ) );
        }



    TGraphErrors *grIneffFromMult[1000];
    int cW = 0;
    double x, y;
    for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
    {
        //mult hist: need <n> for cBin
        float cBinMin = cStep[cW] * cBin;
        float cBinMax = cWidths[cW] + cStep[cW] * cBin;

        float cBinCentre = cBinMin + (cBinMax - cBinMin)/2;


        grIneffFromMult[cBin] = new TGraphErrors;
        TString strGrName = Form( "grIneffFromMult_cBin%d", cBin );
        grIneffFromMult[cBin]->SetName( strGrName.Data() );
        tuneGraphAxisLabels( grIneffFromMult[cBin] );

        //        double correctedMeanNf = 0;

        for ( int fId = 0; fId < nFiles; fId++ )
        {
            TString multNameF = Form( "hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
            TH1D *histMultF = (TH1D*)f[fId]->Get( multNameF );

            TString multNameB = Form( "hist1D_multDistrB_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
            TH1D *histMultB = (TH1D*)f[fId]->Get( multNameB );

            //histMult->DrawCopy( cBin == 0 ? "" : "same" );
            double meanNinCentrBin = histMultF->GetMean();// + histMultB->GetMean();

            if ( fId == 0 ) //remember mean nF for correcitons (below)
            {
                //                correctedMeanNf = meanNinCentrBin*corrFactors[cBin];
            }

            //get bCorr in this cBin
            graphs_ineff[fId][cW][etaId]->GetPoint( cBin, x, y );
            double error = graphs_ineff[fId][cW][etaId]->GetErrorY(cBin);

            //set point on ineff dep graph
            //MODIFY X: mult->INEFF
            if (0)
                grIneffFromMult[cBin]->SetPoint( fId, 0.85-0.15*ineffValues[fId], y );
            else
                grIneffFromMult[cBin]->SetPoint( fId, meanNinCentrBin, y );

            grIneffFromMult[cBin]->SetPointError( fId, 0, error );
        }
        //        drawGraph( grIneffFromMult[cBin], 20, kOrange-9+cBin, "PL");//cBin==0 ? "APL" : "PL" );
        drawGraph( grIneffFromMult[cBin], 20, kWhite, "");//cBin==0 ? "APL" : "PL" );

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

        // GRAPH WITH HIGH-STAT POINTS
        double meanNinCentrBin = 0;
        double meanNinCentrBinError = 0;
        double correctedMeanNf = 0;
        double pointError = 0;
        if (1)
        {
            //high-stat graph modification to have multF on x-axis:
            TString multNameF_highStatGraph = Form( "hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
            TH1D *histMultF_highStatGraph = (TH1D*)fileHighStat->Get( multNameF_highStatGraph );

            //CORRECT meanNf
            meanNinCentrBin = histMultF_highStatGraph->GetMean();// + histMultB->GetMean();
            meanNinCentrBinError = histMultF_highStatGraph->GetMeanError();// + histMultB->GetMean();
            correctedMeanNf = meanNinCentrBin*corrFactors[cBin];

            double x1, y1;
            graphsHighStat_cBin[cW][etaId]->GetPoint( cBin, x1, y1 );
            graphsHighStat_mult[cW][etaId]->SetPoint( cBin, meanNinCentrBin, y1 );

            //re-do fit using high-stat point
            TGraphErrors grForFit;
            grForFit.SetPoint(0, meanNinCentrBin, y1 );
            pointError = graphsHighStat_mult[cW][etaId]->GetErrorY( cBin );
            grForFit.SetPointError(0, 0, pointError );
            if(1)for ( int i = 0; i < 2; i++ )
                grForFit.Fit( fFitFunc,"Q");

        }



        fFitFunc->SetLineColor(kOrange-9+cBin);
        //fFitFunc->SetLineStyle( fileId==0 ? 1 : 2 );
        fFitFunc->DrawCopy("same");

        // draw variation of the fit
        if (0)
        {
            double a = fFitFunc->GetParameter(0);
            double aErr = fFitFunc->GetParError(0);

            fFitFunc->SetLineColor(kGray);//-9+cBin);

            //upper bound
            fFitFunc->SetParameter(0,a+aErr);
            fFitFunc->DrawCopy("same");

            //lower bound
            fFitFunc->SetParameter(0,a-aErr);
            fFitFunc->DrawCopy("same");

            //set back the a-value
            fFitFunc->SetParameter(0,a);
        }

        // GRAPH WITH CORRECTED POINTS
        if (1)
        {
            //corrected point - from correspondingFitLineValue
            double corrected_bCor = fFitFunc->Eval( correctedMeanNf );
            //double propagatedPointError = pointError*correctedMeanNf/meanNinCentrBin;
            //DO NOT TRY TO MULTIPLY ERROR BY FACTOR - SIMPLY COPY RAW STAT ERROR!
            double propagatedPointError = pointError; //*correctedMeanNf/meanNinCentrBin;
            graphsCorrected_mult[cW][etaId]->SetPoint(cBin, correctedMeanNf, corrected_bCor);
            graphsCorrected_mult[cW][etaId]->SetPointError( cBin, 0, propagatedPointError );


            graphsCorrected_cBin[cW][etaId]->SetPoint(cBin, cBinCentre, corrected_bCor);
            graphsCorrected_cBin[cW][etaId]->SetPointError( cBin, meanNinCentrBinError, propagatedPointError );
        }

        //QA plot: abs error of bootstrap VS centrality bins
        if ( 1 ) //cBin < nCentrBins[cW] )
        {
            graphsHighStat_BS_abs_error[cW][etaId]->SetPoint( cBin, cBinCentre, pointError );
            graphsHighStat_BS_abs_error[cW][etaId]->SetPointError( cBin, 0, 0 );
        }

    }

    //TMP for the plot: remove last two points from graphsHighStat_mult, because these are for 80-90 cBins
    if (1)
    {
        graphsHighStat_mult[cW][etaId]->RemovePoint( graphsHighStat_mult[cW][etaId]->GetN()-1 );
        graphsHighStat_mult[cW][etaId]->RemovePoint( graphsHighStat_mult[cW][etaId]->GetN()-1 );

        graphsHighStat_cBin[cW][etaId]->RemovePoint( graphsHighStat_cBin[cW][etaId]->GetN()-1 );
        graphsHighStat_cBin[cW][etaId]->RemovePoint( graphsHighStat_cBin[cW][etaId]->GetN()-1 );
    }

    // GRAPH WITH HIGH-STATISTICS POINTS
    drawGraph( graphsHighStat_mult[cW][etaId], 30, kBlue+1, "P" );

    // DRAW GRAPH WITH CORRECTED POINTS
    drawGraph( graphsCorrected_mult[cW][etaId], 29, kRed, "PL" );

    TLegend *leg = new TLegend(0.7,0.75,0.95,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry( graphsHighStat_mult[cW][etaId], "uncorrected", "p");
    leg->AddEntry( graphsCorrected_mult[cW][etaId], "corrected", "p");
    leg->Draw();

    gPad->SetGridy();

    // ##### Draw raw and corrected points AS FUNC OF CENTRALITY PERC
    TCanvas *canv_raw_corrected = new TCanvas("canv_raw_corrected","canv_raw_corrected",350,180,700,600 );
    tuneCanvas(canv_raw_corrected);
    tuneGraphAxisLabels( graphsCorrected_cBin[cW][etaId] );
    graphsCorrected_cBin[cW][etaId]->SetTitle( ";centrality percentile;b_{corr}" );
    graphsCorrected_cBin[cW][etaId]->GetXaxis()->CenterTitle();
    graphsCorrected_cBin[cW][etaId]->GetYaxis()->CenterTitle();

    drawGraph( graphsCorrected_cBin[cW][etaId], 21, kRed, "APL" );
    drawGraph( graphsHighStat_cBin[cW][etaId], 20, kBlue, "PL" );

    leg = new TLegend(0.7,0.75,0.95,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry( graphsHighStat_cBin[cW][etaId], "uncorrected", "p");
    leg->AddEntry( graphsCorrected_cBin[cW][etaId], "corrected", "p");
    leg->Draw();

    gPad->SetGridy();


    // draw bootstrap abs error as func of centrality bin
    TCanvas *canv_BS_abs_error = new TCanvas("canv_BS_abs_error","canv_BS_abs_error",450,150,700,600 );
    tuneCanvas(canv_BS_abs_error);

    for ( int cW = 0; cW < nCW; cW++ )
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
            drawGraph( graphsHighStat_BS_abs_error[cW][etaId], 20, kBlue, cW == 0 ? "APL" : "PL" );


    return;
}
