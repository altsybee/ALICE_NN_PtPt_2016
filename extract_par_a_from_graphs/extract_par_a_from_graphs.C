#include <TFile.h>
#include "TGraph.h"
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
    TCanvas *canv_graphCorr = new TCanvas("canv_graphCorr","canv_graphCorr",20,50,700,600 );
    tuneCanvas(canv_graphCorr);

    TCanvas *canv_fit_graph_with_func = new TCanvas("canv_fit_graph_with_func","canv_fit_graph_with_func",200,50,700,600 );
    tuneCanvas(canv_fit_graph_with_func);


    TGraphErrors *grA[2];


    const int nCW = 1; //nCentrWidths
    const double cWidths[] = { 10, 5 /*, 2.5, 1.0, 0.5*/ }; //width of the centrality bins
    const double cStep[] = { 5, 2.5/*,  2.5, 1.0, 1.0*/ }; //centrality bins step
    const int nCentrBins[] = { 17, 35/*, 36,  90, 90*/ }; //n centrality bins

    int colors[] = { kBlack, kBlack, kBlue, kRed, kMagenta };
    int markers[] = { 20, 24, 5, 2, 3 };

    const int nEtaWins = 1;

    const int kCorrType = 1; //0-NN, 1-PtPt, 2-PtN

    TFile *f[8];
    f[0] = new TFile( "../bootstrap_ptpt/output_histos_graphs_ALL_14mln.root" );
//    f[0] = new TFile( "../data_graphs_23_02_2016/new_eW_binning_27_02_2016/output_histos_graphs_LHC10h_eW04_MFplus.root" );
    f[1] = new TFile( "../data_graphs_23_02_2016/output_classesByV0M_LHC15o_fieldCOMBINED_c10_5_CUT_OUTLIERS.root" );

    TGraphErrors *graphCorr[nCW][nEtaWins];

    TGraphErrors *grNN[nCW][nEtaWins];
    TGraphErrors *grPtPt[nCW][nEtaWins];
    TGraphErrors *grPtN[nCW][nEtaWins];


    TLegend *legFits = new TLegend(0.65,0.2,0.999,0.45);



    //    int fileId = 0;
    for ( int fileId = 0; fileId < 2; fileId++ )
    {
        canv_graphCorr->cd();
        canv_graphCorr->Update();

        for ( int cW = 0; cW < nCW; cW++ )
            for ( int etaW = 0; etaW < nEtaWins; etaW++ )
            {
                grNN[cW][etaW] = (TGraphErrors*)f[fileId]->Get( Form( "grNN_c%d_eta%d", cW, etaW ) );
                if (fileId==1) grPtPt[cW][etaW] = (TGraphErrors*)f[fileId]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
                if (fileId==0) grPtPt[cW][etaW] = (TGraphErrors*)f[fileId]->Get( Form( "gr_BS_PtPt_cW%d_etaW%d_phiW0", cW, etaW ) );
                grPtN[cW][etaW] = (TGraphErrors*)f[fileId]->Get( Form( "grPtN_c%d_eta%d", cW, etaW ) );


                //            cout << grNN[cW][etaW] << endl;

                if ( kCorrType == 0) graphCorr[cW][etaW] = grNN[cW][etaW];
                else if ( kCorrType == 1) graphCorr[cW][etaW] = grPtPt[cW][etaW];
                else if ( kCorrType == 2) graphCorr[cW][etaW] = grPtN[cW][etaW];

                cout << graphCorr[cW][etaW]->GetName() << endl;
                TString grName = Form( "%s_file%d", graphCorr[cW][etaW]->GetName(), fileId );
                graphCorr[cW][etaW]->SetName( grName.Data() );
            }
        const int etaId = 0;

        // ########## Graphs
        canv_graphCorr->cd();
        graphCorr[0][etaId]->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
        tuneGraphAxisLabels(graphCorr[0][etaId]);


        //remove some points by hand:
        //    graphCorr[0][etaId]->RemovePoint( nCentrBins[0]-1 );
        //    if (nCW>1) graphCorr[1][etaId]->RemovePoint( nCentrBins[1]-1 );

        //draw for several centralities:
        for ( int cW = 0; cW < nCW; cW++ )
            drawGraph(graphCorr[cW][etaId], markers[cW], colors[cW], (cW == 0 && fileId == 0) ? "AP" : "P" );

        graphCorr[0][etaId]->SetMinimum( 0 );

        TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
        leg->SetFillColor(kWhite);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        for ( int cW = 0; cW < nCW; cW++ )
            leg->AddEntry(graphCorr[cW][etaId], Form("class width %.1f", cWidths[cW]), "p");
        if ( fileId == 0 )
            leg->Draw();

        //    graphCorr[0][etaId]->GetYaxis()->SetRangeUser( 0, 0.92 );

        //    TLatex *tex = new TLatex(0.4,0.89, "#eta_{gap}=0.8, #delta#eta=0.4");
        TLatex *tex = new TLatex(0.5,0.3, "#eta_{gap}=0.8, #delta#eta=0.4");
        if ( fileId == 0 )
            drawTex(tex, 0.045);


        // ###### extract a:
        canv_fit_graph_with_func->cd();
        canv_fit_graph_with_func->Update();
        canv_fit_graph_with_func->SetGrid(1,1);

        TGraphErrors *grFrame = new TGraphErrors;
        grFrame->SetPoint(0,0,0);
        grFrame->SetPoint(1,0,0.2);
        grFrame->SetPoint(2,1000,0.2);
        grFrame->SetPoint(3,1000,0);
        grFrame->SetMarkerColor(kWhite);

        tuneGraphAxisLabels( grFrame );
        grFrame->SetTitle(";N_{ch};b_{corr}");
        if ( fileId == 0 )
            grFrame->Draw("AP");


        const int cW = 0;
        TGraphErrors *gr = graphCorr[cW][etaId];

        grA[fileId] = new TGraphErrors;
        grA[fileId]->SetTitle(";centrality percentile;coeff a");
        grA[fileId]->SetName( Form("grA_file%d", fileId ) );
        double x, y;
        for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
        {
            //mult hist: need <n> for cBin
            float cBinMin = cStep[cW] * cBin;
            float cBinMax = cWidths[cW] + cStep[cW] * cBin;

            TString multNameF = Form( "hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
            TH1D *histMultF = (TH1D*)f[fileId]->Get( multNameF );

            TString multNameB = Form( "hist1D_multDistrB_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaId, 0);
            TH1D *histMultB = (TH1D*)f[fileId]->Get( multNameB );

            //histMult->DrawCopy( cBin == 0 ? "" : "same" );
            double meanNinCentrBin = histMultF->GetMean() + histMultB->GetMean();


            //get bCorr in this cBin
            gr->GetPoint( cBin, x, y );
            double error = gr->GetErrorY(cBin);

            TGraphErrors *grForFit = new TGraphErrors;
            TString strGrFit = Form( "grFit_cBin%d_file%d", cBin, fileId);
            grForFit->SetName( strGrFit.Data() );

            grForFit->SetPoint(0, 0, 0);
            grForFit->SetPoint(1, meanNinCentrBin, y );
            grForFit->SetPointError(1, 0, error );

            grForFit->SetMarkerStyle( fileId==0 ? 24 : 25 );//29 : 30 );
//            grForFit->SetMarkerSize( fileId==0 ? 1.5 : 1.4 );
            grForFit->SetMarkerSize( 0.8 );
            grForFit->SetMarkerColor( fileId==0 ? kRed : kBlue );
            grForFit->SetLineColor( fileId==0 ? kRed : kBlue );
            grForFit->DrawClone( cBin==0? "P" : "P" );

            if (cBin==0 && fileId==0)
                legFits->AddEntry(grForFit, "2.76 TeV", "p");
            else if (cBin==0 && fileId==1)
                legFits->AddEntry(grForFit, "5.02 TeV", "p");

            //fit func
            TF1 *fFitFunc = new TF1("fitFunc","[0]*x/(1+[0]*x)",0,1200);
//            TF1 *fFitFunc = new TF1("fitFunc","x/([0]+x)",0,1200);
//            TF1 *fFitFunc = new TF1("fitFunc", "x/([0]+x)",0,1200);
            fFitFunc->SetParameter(0,0.01);


            for ( int i = 0; i < 2; i++ )
                grForFit->Fit( fFitFunc,"Q");

            fFitFunc->SetLineColor(kOrange-5+cBin);
            fFitFunc->SetLineStyle( fileId==0 ? 1 : 2 );
            fFitFunc->DrawCopy("same");


            //get parameter a
            double a = fFitFunc->GetParameter(0);
            double aErr = fFitFunc->GetParError(0);
//            grA[fileId]->SetPoint(cBin, x, a/ meanNinCentrBin );
//            grA[fileId]->SetPointError(cBin, 0, aErr/ meanNinCentrBin );
            grA[fileId]->SetPoint(cBin, x, a );
            grA[fileId]->SetPointError(cBin, 0, aErr );

        } // end of cBin loop
    } //end of files loop


    legFits->SetFillColor(kWhite);
    legFits->SetFillStyle(0);
    legFits->SetBorderSize(0);
    legFits->Draw();


    TCanvas *canv_par_a = new TCanvas("canv_par_a","canv_par_a",280,300,700,600 );
    tuneCanvas(canv_par_a);
    tuneGraphAxisLabels( grA[0] );
    //        grA[0]->SetMarkerStyle(20);
    //        grA[fileId]->Draw( (fileId == 0) ? "APL" : "PL" );

    drawGraph( grA[0], markers[0], colors[0], "APL" );
    drawGraph( grA[1], markers[1], colors[1], "PL" );

    //    grA[0]->SaveAs( "graph_A_LHC10h_eta0.root" );
    //    grA[1]->SaveAs( "graph_A_LHC15o_eta0.root" );


    //RATIO A:
    TCanvas *canv_ratio_a = new TCanvas("canv_ratio_a","canv_ratio_a",320,300,700,600 );
    tuneCanvas(canv_ratio_a);

    calcPointsRatio( grA[0], grA[1] );
    tuneGraphAxisLabels( grA[0] );
    drawGraph( grA[0], markers[0], colors[0], "APL" );


    return;
}
