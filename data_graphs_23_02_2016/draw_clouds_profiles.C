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

void removeBinsWithFewEntries( TH2D *hist2D )
{
    for ( int i = 1; i < hist2D->GetNbinsX()+1; i++ )
        for ( int j = 1; j < hist2D->GetNbinsY()+1; j++ )
        {
            if ( hist2D->GetBinContent( i, j ) <= 1 )
                hist2D->SetBinContent( i, j, 0 );
            //            cout << i << " " << j << endl;
        }
}

void deleteProfileEmptyBinErrors( TProfile *pr )
{
    for ( int i = 1; i < pr->GetNbinsX()+1; i++ )
    {
        if ( /*pr->GetBinError(i) < 0.0001 && */
             pr->GetBinEntries(i) < 50 )
            //                pr->SetBinError(i, 0.001);//0.01);
            pr->SetBinEntries(i,0);// Error(i, 0.1);
        //            cout << i << ": " << pr->GetBinEntries(i) << " " << pr->GetBinError(i) << endl;
    }
    //    cout << endl;
}



//TH1D* getHist1DWithProperErrors(TH2D* source, const char *name)
//{
//    TH1D* fPrAbs = source->ProjectionX();
//    fPrAbs->Reset();

//    // Calculate errors for NN
//    TProfile* profX = (TProfile*) source->ProfileX(name, 1, source->GetNbinsY());
//    for(int i = 0; i < profX->GetNbinsX(); i++)
//    {
//        fPrAbs->SetBinContent(i, profX->GetBinContent(i));
//        if(fPrf->GetBinContent(i)!=0)
//        {
//            fPrAbs->SetBinError(i,sqrt(profX->GetBinContent(i)/fPrf->GetBinContent(i)));
//        }
//    }
//    return fPrAbs;
//}


void tunePad(TVirtualPad *pad)
{
    pad->SetLeftMargin(0.18);
    pad->SetRightMargin(0.05);
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.11);
    //    pad->SetGridy();
}

void tuneHist2D_onPad( TH2D *hist )
{
    hist->GetYaxis()->SetTitleOffset( 1.45 );
    hist->GetYaxis()->SetTitleSize( 0.048 );
    hist->GetYaxis()->SetLabelSize( 0.042 );

    hist->GetXaxis()->SetTitleOffset( 1.05 );
    hist->GetXaxis()->SetTitleSize( 0.048 );
    hist->GetXaxis()->SetLabelSize( 0.042 );
}
void tuneProfile_onPad( TProfile *prof )
{
    prof->GetYaxis()->SetTitleOffset( 1.45 );
    prof->GetYaxis()->SetTitleSize( 0.048 );
    prof->GetYaxis()->SetLabelSize( 0.042 );

    prof->GetXaxis()->SetTitleOffset( 1.05 );
    prof->GetXaxis()->SetTitleSize( 0.048 );
    prof->GetXaxis()->SetLabelSize( 0.042 );
}


void draw_clouds_profiles()
{
    //        gROOT->ProcessLine(".L ../utils.C");
    //        gROOT->ProcessLine(".L AliLRCFit.cxx");
    TFile *f[8];

    //    f[0] = new TFile( "output_classesByV0M_LHC10h.root" );
    //    f[0] = new TFile( "output_classesByV0M_LHC10h_c10_5_1.root" );
        f[0] = new TFile( "output_classesByV0M_LHC10h_c10_5_25_1_05.root" );
//    f[0] = new TFile( "output_classesByV0M_LHC11h_FemtoPlus_c10_5_CUT_OUTLIERS.root" );

//    f[0] = new TFile( "output_classesByV0M_LHC15o_fieldMM_c10_5_CUT_OUTLIERS.root" );
//    f[0] = new TFile( "output_classesByV0M_LHC15o_fieldPP_c10_5_CUT_OUTLIERS.root" );

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

        const int nCW = 5; //nCentrWidths
        const double cWidths[nCW] = { 10, 5, 2.5, 1.0, 0.5 }; //width of the centrality bins
        const double cStep[nCW] = { 5, 2.5,  2.5, 1.0, 1.0 }; //centrality bins step
        const int nCentrBins[nCW] = { 17, 35, 36,  90, 90 }; //n centrality bins


    TH2D *hist2D;//[200][3];
    TProfile *profile;//[200][3];
    int cW = 2;
    int etaW = 1;
    int phiW = 0;

    const int kCorrType = 1; //0-NN, 1-PtPt, 2-PtN

    TCanvas *canv_tmp_for_fit = new TCanvas("canv_tmp_for_fit","canv_tmp_for_fit",50,50,300,300 );
    TCanvas *canv_2D_clouds = new TCanvas("canv_2D_clouds","canv_2D_clouds",150,250,1400,600 );
    tuneCanvas(canv_2D_clouds);
    canv_2D_clouds->Divide(2,1);

    gStyle->SetOptStat( kFALSE );

    TGraph *grFromFit2D = new TGraph;

    bool firstDraw = true;
    //    for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
    for ( int cBin = nCentrBins[cW]-1; cBin >= 0; cBin-- )
    {
                if (cBin%2!=0)
                    continue;

        //        cout << "cBin=" << cBin << endl;

        float cBinMin = cStep[cW] * cBin;
        float cBinMax = cWidths[cW] + cStep[cW] * cBin;

        // ##### pad 1 - clouds
        tunePad( canv_2D_clouds->cd(1) );
        if ( kCorrType == 0 )
        {
            hist2D = (TH2D*)f[0]->Get( Form("hist2D_NN_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW) );
            hist2D->SetTitle( "");
            hist2D->GetXaxis()->SetTitle( "N_{ch} Forward");
            hist2D->GetYaxis()->SetTitle( "N_{ch} Backward");
            hist2D->GetXaxis()->SetRangeUser(0,650);
            hist2D->GetYaxis()->SetRangeUser(0,650);
        }
        else if ( kCorrType == 1 )
        {
            hist2D = (TH2D*)f[0]->Get( Form("hist2D_PtPt_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW) );
            hist2D->SetTitle( "");
            hist2D->GetXaxis()->SetTitle( "#LTp_{T}#GT Forward");
            hist2D->GetYaxis()->SetTitle( "#LTp_{T}#GT Backward");
        }     hist2D->SetMarkerColor(kOrange-9+cBin);
        tuneHist2D_onPad(hist2D);
        hist2D->GetXaxis()->CenterTitle();
        hist2D->GetYaxis()->CenterTitle();

        //        removeBinsWithFewEntries(hist2D);

        hist2D->DrawCopy( firstDraw ? "" : "same" );

        // ##### pad 2 - profiles
        tunePad( canv_2D_clouds->cd(2) );
        profile = hist2D->ProfileX();    //(TProfile*)f[0]->Get( Form("hist2D_c%.1f-%.1f_etaW_%d_phiW_%d_pfx", cBinMin, cBinMax, etaW, phiW) );
        if ( kCorrType == 0 )
        {
            profile->SetTitle( "");
            profile->GetYaxis()->SetTitle( "#LTN_{ch}#GT Backward");
            profile->GetXaxis()->SetRangeUser(0,650);
            profile->GetYaxis()->SetRangeUser(0,650);
        }
        else if ( kCorrType == 1 )
        {
            profile->SetTitle( "");
            profile->GetYaxis()->SetTitle( "#LT#LTp_{T}#GT#GT Backward");
        }
        profile->SetLineColor(kOrange-9+cBin);
        profile->SetMarkerStyle(7);
        tuneProfile_onPad( profile );
        profile->GetYaxis()->CenterTitle();

        deleteProfileEmptyBinErrors(profile);


        canv_tmp_for_fit->cd();
        profile->Fit("pol1","Q");//,"",0.25,1.2);//,"N");
        canv_2D_clouds->cd(2);

        TF1 *fit = profile->GetFunction("pol1");
        Double_t p0 = fit->GetParameter(0);
        Double_t p1 = fit->GetParameter(1);
        grFromFit2D->SetPoint(grFromFit2D->GetN(), (cBinMax+cBinMin)/2, p1);

        double meanX = hist2D->ProjectionX()->GetMean();
        double rmsX = hist2D->ProjectionX()->GetRMS();
        fit->SetRange( meanX-3*rmsX, meanX+3*rmsX );
        fit->SetLineColorAlpha( kRed, 0.6 );

        fit->Draw("same");

        profile->DrawCopy( firstDraw ? "" : "same" );



        firstDraw = false;
    }

    TGraphErrors *grByFormula; /*[cW][etaW]*/
    if ( kCorrType == 0 )
    {
        grByFormula = (TGraphErrors*)f[0]->Get( Form( "grNN_c%d_eta%d", cW, etaW ) );
    }
    else if ( kCorrType == 1 )
    {
        grByFormula = (TGraphErrors*)f[0]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
    }



    TCanvas *canv_GrCoeff = new TCanvas("canv_GrCoeff","canv_GrCoeff",20,150,900,700 );
    grByFormula->Draw("APL");

    grFromFit2D->SetLineColor(kRed);
    grFromFit2D->DrawClone("PL");



    return;

    TGraphErrors *gr[10][10];
    for ( int cW = 0; cW < 2; cW++ )
        for ( int etaW = 0; etaW < 3; etaW++ )
            gr[cW][etaW] = (TGraphErrors*)f[0]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );



    drawGraph( gr[1][0], 24, kBlack, "APL" );
    drawGraph( gr[0][0], 20, kBlack, "PL" );

    //    drawGraph( gr[1][1], 24, kBlue, "PL" );
    //    drawGraph( gr[0][1], 20, kBlue, "PL" );

    //    drawGraph( gr[1][2], 24, kGreen, "PL" );
    //    drawGraph( gr[0][2], 20, kGreen, "PL" );

    ////    drawGraph( gr[6], 24, kRed, "PL" );
    //    drawGraph( gr[7], 20, kRed, "PL" );



    f[1] = new TFile( "output_histos_graphs_LHC15o_fieldMM.root" );
    TGraphErrors *grMM[10][10];
    for ( int cW = 0; cW < 2; cW++ )
        for ( int etaW = 0; etaW < 3; etaW++ )
            grMM[cW][etaW] = (TGraphErrors*)f[1]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );

    drawGraph( grMM[1][0], 24, kGreen, "PL" );
    drawGraph( grMM[0][0], 20, kGreen, "PL" );

    //    drawGraph( grMM[1][1], 24, kRed, "PL" );
    //    drawGraph( grMM[0][1], 20, kRed, "PL" );

    //    drawGraph( grMM[1][2], 24, kGreen, "PL" );
    //    drawGraph( grMM[0][2], 20, kGreen, "PL" );




    //    gROOT->ProcessLine( ".q");
}
