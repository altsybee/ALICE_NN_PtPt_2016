#include <TTree.h>
#include <TBranch.h>
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
#include "TRandom3.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;




void tuneGraphAxisLabels( TGraphErrors *gr )
{
    gr->GetYaxis()->SetTitleOffset( 1.45 );
    gr->GetYaxis()->SetTitleSize( 0.048 );
    gr->GetYaxis()->SetLabelSize( 0.042 );

    gr->GetXaxis()->SetTitleOffset( 0.95 );
    gr->GetXaxis()->SetTitleSize( 0.048 );
    gr->GetXaxis()->SetLabelSize( 0.042 );

}


void tuneCanvas(TCanvas *canvas)
{
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.11);
    canvas->SetGridy();
}

void drawGraph(TGraphErrors *gr, int marker, int color, const char* drawOpt =  "P same", double markerSize = -1 )
{
    gr->SetMarkerStyle( marker );
    gr->SetMarkerColor( color );
    gr->SetLineColor( color );
    if (markerSize>0)
        gr->SetMarkerSize( markerSize );
    gr->DrawClone( drawOpt );
}

void drawTex(TLatex *tex, double fontSize = 0.045)
{
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextSize(fontSize);
    //    tex->SetLineWidth(2);
    tex->DrawClone();

}


void getQuantiles(TH1D *h, const int nq, double *yq)
{
    cout << "nq=" << nq << endl;
    //    const Int_t nq = 10;
    //    const Int_t nshots = 10;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    //    Double_t yq[nq];  // array to contain the quantiles
    for (Int_t i=0;i<nq;i++)
        xq[i] = Float_t(i+1)/nq;

    //    TGraph *gr70 = new TGraph(nshots);
    //   TH1D *h = new TH1D("h","demo quantiles",50,-3,3);

    h->GetQuantiles(nq,yq,xq);
    for (Int_t i=0;i<nq;i++)
        cout << yq[i] << ", ";
    cout << endl;

    //    TCanvas *c1 = new TCanvas("c1","demo quantiles",10,10,600,900);
    //    TGraph *gr = new TGraph(nq,xq,yq);
    //    gr->SetTitle("final quantiles");
    //    gr->SetMarkerStyle(21);
    //    gr->SetMarkerColor(kRed);
    //    gr->SetMarkerSize(0.3);
    //    gr->Draw("ap");
}


void drawCanvasWithClasses(TH1D *hist1D, TString label, const int nCentrBins, double *centrMultBounds)
{
    TCanvas *canv_mult_withClasses = new TCanvas(Form("canv_mult_withClasses_%s",label.Data())
                                                 ,Form("canv_mult_withClasses_%s",label.Data()),300,100,800,600 );
    tuneCanvas(canv_mult_withClasses);

    hist1D->SetLineColor(kRed);
    hist1D->DrawCopy();

    gPad->SetLogy();

    //draw centrality classes on mult hist
    //    int nBinsForClassHists = fHistParticlesInCutConditionInEvent->GetNbinsX();
    TH1D **fHistCentrClass = new TH1D*[nCentrBins];
    for ( int iCentrClass = 0; iCentrClass < nCentrBins; iCentrClass++ )
    {
        fHistCentrClass[iCentrClass] = (TH1D*)hist1D->Clone( Form("fHistCentrClass%d", iCentrClass) );   //new TH1D( Form("fHistCentrClass%d", iCentrClass), Form("iCentrClass%d", iCentrClass)
        //                                                 , hist1D->GetNbinsX(), hist1D->GetBinLowEdge(1), hist1D->GetBin );
        fHistCentrClass[iCentrClass]->Reset();
    }

    int centrClassId = 0;
    for ( int iBin = 0; iBin < hist1D->GetNbinsX(); iBin++ )
    {
        //        cout << fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) << " " << centrMultBounds[centrClassId] << endl;
        if ( hist1D->GetBinCenter(iBin+1) > centrMultBounds[centrClassId] )
            centrClassId++;
        if ( centrClassId >= nCentrBins )
            break;
        double binContent = hist1D->GetBinContent( iBin+1 );
        fHistCentrClass[centrClassId]->SetBinContent( iBin+1, binContent );
    }
    for ( int iCentrClass = 0; iCentrClass < nCentrBins; iCentrClass++ )
    {
        fHistCentrClass[iCentrClass]->SetFillColor( kOrange - 5 + iCentrClass );
        fHistCentrClass[iCentrClass]->SetLineColor(kBlue);
        fHistCentrClass[iCentrClass]->DrawCopy( "same" );

        //                double meanNchInCentrBin = fHistCentrClass[iCentrClass]->GetMean();
        //                fHistMultClassMeanNch->SetBinContent( iCentrClass+1, meanNchInCentrBin );
        //                cout << "meanNch in centrality bin " << iCentrClass << ": " << meanNchInCentrBin << endl;
    }

    canv_mult_withClasses->SaveAs( Form("output/canv_%s_classes_%d.eps", label.Data(), nCentrBins));

}




class WinPair
{
public:
    float cBinMin;
    float cBinMax;

    double NN_nF;
    double NN_nB;
    double NN_nF_nB;
    double NN_nF2;
    int NN_Nevents;
    TH2D *hist2D_NN;

    double PtPt_PtF;
    double PtPt_PtB;
    double PtPt_PtF_PtB;
    double PtPt_PtF2;
    int PtPt_Nevents;
    TH2D *hist2D_PtPt;

    double PtN_nF;
    double PtN_PtB;
    double PtN_nF_PtB;
    double PtN_nF2;
    int PtN_Nevents;
    TH2D *hist2D_PtN;

    //corr coeffs
    double NN_bCorr;
    double NN_C2;
    double PtPt_bCorr;
    double PtPt_C2;
    double PtN_bCorr;
    double PtN_C2;

    TH1D *hist1D_EstimatorEntries;
    TH1D *hist1D_multDistrF;
    TH1D *hist1D_multDistrB;

    WinPair() :
        cBinMin(0)
      , cBinMax(100)
      , NN_nF(0)
      , NN_nB(0)
      , NN_nF_nB(0)
      , NN_nF2(0)
      , NN_Nevents(0)
      , PtPt_PtF(0)
      , PtPt_PtB(0)
      , PtPt_PtF_PtB(0)
      , PtPt_PtF2(0)
      , PtPt_Nevents(0)
      , PtN_nF(0)
      , PtN_PtB(0)
      , PtN_nF_PtB(0)
      , PtN_nF2(0)
      , PtN_Nevents(0)

      , NN_bCorr(-1000)
      , NN_C2(-1000)
      , PtPt_bCorr(-1000)
      , PtPt_C2(-1000)
      , PtN_bCorr(-1000)
      , PtN_C2(-1000)
    {}
    void init(float _cMin, float _cMax, int etaW, int phiW)
    {
        cBinMin = _cMin;
        cBinMax = _cMax;

        TString strNN = Form("hist2D_NN_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
        hist2D_NN = new TH2D( strNN, strNN, 800, -0.5, 799.5, 800, -0.5, 799.5 );

        TString strPtPt = Form("hist2D_PtPt_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
        hist2D_PtPt = new TH2D( strPtPt, strPtPt, 100, 0, 2, 100, 0, 2);

        TString strPtN = Form("hist2D_PtN_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
        hist2D_PtN = new TH2D( strPtN, strPtN, 1000, -0.5, 999.5, 1000, 0, 2);

        //QA histos:
        TString strEstPerc = Form("hist1D_EstimatorEntries_c%.1f-%.1f_etaW_%d_phiW_%d;percentile;entries", cBinMin, cBinMax, etaW, phiW);
        hist1D_EstimatorEntries = new TH1D( strEstPerc, strEstPerc, 1001, -0.5, 100.5);

        TString strMultDistrF = Form("hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d;percentile;entries", cBinMin, cBinMax, etaW, phiW);
        hist1D_multDistrF = new TH1D( strMultDistrF, strMultDistrF, 601, -0.5, 600.5);

        TString strMultDistrB= Form("hist1D_multDistrB_c%.1f-%.1f_etaW_%d_phiW_%d;percentile;entries", cBinMin, cBinMax, etaW, phiW);
        hist1D_multDistrB = new TH1D( strMultDistrB, strMultDistrB, 601, -0.5, 600.5);
    }
    void fill( Float_t cPerc, UShort_t nF, UShort_t nB, Float_t PtF, Float_t PtB )
    {
        //check if we in centrality bin
        if ( cPerc < cBinMin || cPerc > cBinMax )
            return;

        //QA fills:
        hist1D_EstimatorEntries->Fill(cPerc);
        hist1D_multDistrF->Fill(nF);
        hist1D_multDistrB->Fill(nB);

        //NN
        hist2D_NN->Fill( nF, nB );

        NN_nF       += nF;
        NN_nB       += nB;
        NN_nF_nB    += nF*nB;
        NN_nF2      += nF*nF;
        NN_Nevents++;

        //PtPt
        if ( nF > 0 && nB > 0 )
        {
            hist2D_PtPt->Fill( PtF, PtB );

            PtPt_PtF += PtF;
            PtPt_PtB += PtB;
            PtPt_PtF_PtB += PtF*PtB;
            PtPt_PtF2 += PtF*PtF;
            PtPt_Nevents++;
        }

        //PtN
        if ( nB > 0 )
        {
            hist2D_PtN->Fill( nF, PtB );

            PtN_nF += nF;
            PtN_PtB += PtB;
            PtN_nF_PtB += nF*PtB;
            PtN_nF2 += nF*nF;
            PtN_Nevents++;
        }
    }

    void calcCorrCoeffs()
    {
        NN_bCorr    = _calc( 0, NN_nF, NN_nB, NN_nF_nB, NN_nF2, NN_Nevents );
        NN_C2       = _calc( 1, NN_nF, NN_nB, NN_nF_nB, NN_nF2, NN_Nevents );

        PtPt_bCorr    = _calc( 0, PtPt_PtF, PtPt_PtB, PtPt_PtF_PtB, PtPt_PtF2, PtPt_Nevents );
        PtPt_C2       = _calc( 1, PtPt_PtF, PtPt_PtB, PtPt_PtF_PtB, PtPt_PtF2, PtPt_Nevents );

        PtN_bCorr    = _calc( 0, PtN_nF, PtN_PtB, PtN_nF_PtB, PtN_nF2, PtN_Nevents, 1 );
        PtN_C2       = _calc( 1, PtN_nF, PtN_PtB, PtN_nF_PtB, PtN_nF2, PtN_Nevents );
    }

private:
    double _calc(int type, double F, double B, double FB, double F2, int nEvents, int ifRel=0)
    {
        if ( nEvents <= 0 )
            return -1000;
        double meanF     =  F   / nEvents;     //printf( "<F> = %f\n", meanF );
        double meanB     =  B   / nEvents;     //printf( "<B> = %f\n", meanB );
        double meanFB    =  FB  / nEvents;    //printf( "<FB> = %f\n", meanFB );
        double meanF2    =  F2  / nEvents;    //printf( "<FB> = %f\n", meanFB );

        //        printf( "nEvents = %d, ", nEvents );
        //        printf( "meanF = %f, ", meanF  );
        //        printf( "meanB = %f, ", meanB  );
        //        printf( "meanFB = %f, ",  meanFB );
        //        printf( "meanF2 = %f\n, ",  meanF2 );

        double numenator = meanFB - meanF * meanB;
        double denominator = 0;
        if (type == 0) // bCorr
            denominator = meanF2 - meanF*meanF;
        else if (type == 1) // C2
            denominator = meanF * meanB;

        double value = 0;
        if (denominator!=0)
            value = numenator / denominator;

        if (ifRel)
            value *= meanF/meanB;


        return value;
    }
};

struct BranchFB
{
    UShort_t nF;
    UShort_t nB;
    Float_t PtF;
    Float_t PtB;
    BranchFB() : nF(0), nB(0), PtF(0), PtB(0)
    {}
};


struct CentralityOccupancy
{
public:
    float cBinMin;
    float cBinMax;

    int nEventsV0M;
    int nEventsZDCZEM;
    int nEventsV0M_and_ZDCZEM;

    CentralityOccupancy() :
        cBinMin(0)
      , cBinMax(100)
      , nEventsV0M(0)
      , nEventsZDCZEM(0)
      , nEventsV0M_and_ZDCZEM(0)
    {}
    void fill( Float_t cPercV0M, Float_t cPercZDCvsZNA )
    {
        //check if we in centrality bin
        bool isV0M = false;
        if ( cPercV0M > cBinMin && cPercV0M < cBinMax )
        {
            nEventsV0M++;
            isV0M = true;
        }

        bool isZDCvsZNA = false;
        if ( cPercZDCvsZNA > cBinMin && cPercZDCvsZNA < cBinMax )
        {
            nEventsZDCZEM++;
            isZDCvsZNA = true;
        }

        if ( isV0M && isZDCvsZNA )
            nEventsV0M_and_ZDCZEM++;
    }

};



void get_nn_ptpt_FROM_FB_TREE() // TString inputFileName = "MergedOutput.root")
{
    //            TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_09_PbPb_Data_AOD_3etaWins_phi1_pt02_20_FB_TREE_blocks123/block1/AnalysisResults.139465.root" );
    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_21_PbPb_Data_AOD_3etaWins_phi1_pt02_20_FB_TREE_TEST_NEW_DATA_5_02TeV_tryAOD2/MergedOutput.root" );
    if (!myFile)
    {
        cout << "No input file!" << endl;
        return;
    }
    myFile->ls();
    TList *listKeys = myFile->GetListOfKeys();
    cout << "going into list: " << listKeys->At(1)->GetName() << endl;

    //return;
    myFile->cd( listKeys->At(1)->GetName() );
    TList *listKeys2 = gDirectory->GetListOfKeys();
    TList *myTask = (TList*) gDirectory->Get( listKeys2->At(0)->GetName() );
    TTree *t1 = (TTree*)  myTask->FindObject( "treeFB" );


    //    cout << t1->GetEntries() << endl;
    //t1->Print();
    //return;
    //    TFile *f = new TFile( "file_event_tree_merged_impPar_2files.root" );
    //    TTree *t1 = (TTree*) f->Get( "t1" );

    //    TFile *f = new TFile( "output_1_list.grid.root" );
    //    TList *list = (TList*) f->Get( "coutput_V0_mult_0.000000_100000.000000" );
    //    TTree *t1 = (TTree*) list->FindObject( "t1" );

    //other vars:
    //    const int nVars = 3;
    //    TString strVarName[] = {
    //        "centr_V0M" ,
    //        "centr_CL1" ,
    //        "centr_ZEMvsZDC"     ,
    //    };
    //    Float_t varBranch[nVars];

    //    for(Int_t var = 0; var < nVars; var++)
    //        t1->SetBranchAddress(strVarName[var],&varBranch[var]);

    Float_t brV0M;
    //    t1->SetBranchAddress( "centr_V0M", &brV0M );
    t1->SetBranchAddress( "centrV0M_NEW_MULT_SEL", &brV0M );

    Float_t brZEMvsZDC;
    t1->SetBranchAddress( "centr_ZEMvsZDC", &brZEMvsZDC );

    //        const int nCW = 1; //nCentrWidths
    //        const double cWidths[nCW] = { 10 }; //width of the centrality bins
    //        const double cStep[nCW] = { 5 }; //centrality bins step
    //        const int nCentrBins[nCW] = { 17 }; //n centrality bins

    const int nCW = 2; //nCentrWidths
    const double cWidths[nCW] = { 10, 5 }; //width of the centrality bins
    const double cStep[nCW] = { 5, 2.5 }; //centrality bins step
    const int nCentrBins[nCW] = { 17, 35 }; //n centrality bins

    //    const int nCW = 5; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 2.5, 1.0, 0.5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5, 2.5, 1.0, 1.0 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35, 36, 90, 90 }; //n centrality bins

    const int nEtaWins = 3;
    const int nPhiWins = 1;

    const int maxNCentrBins = 100; //TMath::MaxElement(nCW, &nCentrBins);
    WinPair wins[nCW][maxNCentrBins][nEtaWins][nPhiWins];
    CentralityOccupancy cOccupancy[nCW][maxNCentrBins];

    for ( int cW = 0; cW < nCW; cW++ )
        for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
        {
            float cBinMin = cStep[cW] * cBin;
            float cBinMax = cWidths[cW] + cStep[cW] * cBin;

            cOccupancy[cW][cBin].cBinMin = cBinMin;
            cOccupancy[cW][cBin].cBinMax = cBinMax;

            for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                    wins[cW][cBin][etaW][phiW].init(cBinMin, cBinMax, etaW, phiW);
        }


    BranchFB br[nEtaWins][nPhiWins];

    for ( int etaW = 0; etaW < nEtaWins; etaW++ )
        for ( int phiW = 0; phiW < nPhiWins; phiW++ )
        {
            int ptW = 0;
            TString brNamePostfix = Form("eta_%d_phi%d_pt_%d"
                                         , etaW, phiW, ptW );
            t1->SetBranchAddress( Form("nF_%s", brNamePostfix.Data() ),  &br[etaW][phiW].nF );
            t1->SetBranchAddress( Form("nB_%s", brNamePostfix.Data() ),  &br[etaW][phiW].nB );
            t1->SetBranchAddress( Form("PtF_%s", brNamePostfix.Data() ), &br[etaW][phiW].PtF );
            t1->SetBranchAddress( Form("PtB_%s", brNamePostfix.Data() ), &br[etaW][phiW].PtB );
        }



    TH1D *hist1D_QA_percentilesEstimator = new TH1D( "hist1D_QA_percentilesEstimator", "hist1D_QA_percentilesEstimator;percentile;entries", 3001, -0.5, 300.5);
    TH1D *hist1D_QA_multALL = new TH1D( "hist1D_QA_multALL", "hist1D_QA_multALL;mult;entries", 3001, -0.5, 3000.5);


    // ##### prepare for tree loop
    int nEvents = t1->GetEntries();
    cout <<"nEvents = " << nEvents << endl;

    //    float **BS_Nf = new float*[nPhiWins];
    //    float **BS_Nb = new float*[nPhiWins];
    //    for ( int w = 0; w < nPhiWins; w++ )
    //    {
    //        BS_Nf[w] = new float[nEvents];
    //        BS_Nb[w] = new float[nEvents];
    //    }

    // ##### main loop over events
    int flag_V0M_ZDC = 0;//1;
    int nAccepted = 0;
    for (Long64_t i=0; i < nEvents; i++)
    {
        if ( i % 100000 == 0 )
            cout << "getting " << (int)i << endl;
        //                cout <<"getting " << (int)i << "\r"; cout.flush();

        t1->GetEntry( i );
        //            t1->GetEntry( TMath::Nint( gRandom->Uniform(-0.5,nEvents-0.5) ) );

        float cEstimator = -1;
        if ( flag_V0M_ZDC==0 )
        {
            if ( brV0M > 90 ) //V0M cut
                continue;
            cEstimator = brV0M;
        }
        else
        {
            if ( brZEMvsZDC > 50 ) //ZDCvsZEM cut
                continue;
            cEstimator = brZEMvsZDC;
        }
        hist1D_QA_percentilesEstimator->Fill(cEstimator);


        //calc mult in whole TPC using wins:
        int multTPC = 0;
        for ( int phiW = 0; phiW < nPhiWins; phiW++ )
            multTPC += (br[0][phiW].nF + br[2][phiW].nF) + (br[0][phiW].nB + br[2][phiW].nB);
        hist1D_QA_multALL->Fill(multTPC);


        //        for ( int w = 0; w < nPhiWins; w++ )
        //        {
        //            BS_Nf[w][nAccepted] = Nf[w];
        //            BS_Nb[w][nAccepted] = Nb[w];
        //        }

        for ( int cW = 0; cW < nCW; cW++ )
            for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
            {
                cOccupancy[cW][cBin].fill(brV0M, brZEMvsZDC);
                for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                    for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                        wins[cW][cBin][etaW][phiW].fill( cEstimator, br[etaW][phiW].nF, br[etaW][phiW].nB, br[etaW][phiW].PtF, br[etaW][phiW].PtB );
            }

        nAccepted++;


    } // end of events
    cout << "nAccepted = " << nAccepted << endl;
    cout << "nAccepted/nAll = " << (float)nAccepted/nEvents << endl;


    // ########## QA PLOTTING:
    TCanvas *canv_estimatorPercentiles_QA_all = new TCanvas("canv_estimatorPercentiles_QA_all","canv_estimatorPercentiles_QA_all",0,0,700,600 );
    hist1D_QA_percentilesEstimator->DrawCopy();

    TCanvas *canv_hist1D_QA_multALL = new TCanvas("canv_hist1D_QA_multALL","canv_hist1D_QA_multALL",50,50,700,600 );
    hist1D_QA_multALL->DrawCopy();


    // MULT BINNING:
    int nCentrBinsMult = 10;
    cout << "nCentrBins=" << nCentrBinsMult << endl;
    double *estBounds = new double[nCentrBinsMult]; // array to contain the quantiles
    getQuantiles(hist1D_QA_multALL, nCentrBinsMult, estBounds);
    drawCanvasWithClasses( hist1D_QA_multALL, "byMultTPC", nCentrBinsMult, estBounds );




    // ##########MAIN PLOTTING FOR CORRS:

    TGraphErrors *grNN [nCW][nEtaWins];
    TGraphErrors *grPtPt[nCW][nEtaWins];
    TGraphErrors *grPtN[nCW][nEtaWins];

    TGraphErrors *grFractEstByV0M[nCW];
    TGraphErrors *grFractEstByZDC[nCW];

    for ( int cW = 0; cW < nCW; cW++ )
    {
        grFractEstByV0M[cW] = new TGraphErrors;
        grFractEstByZDC[cW] = new TGraphErrors;
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
        {
            grNN[cW][etaW] = new TGraphErrors;
            grPtPt[cW][etaW] = new TGraphErrors;
            grPtN[cW][etaW] = new TGraphErrors;
        }
    }


    //calc (1) - occupancies in centr bins, (2) - corr coeffs
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
        {
            if ( cOccupancy[cW][cBin].nEventsV0M > 0 )
            {
                CentralityOccupancy *c = &cOccupancy[cW][cBin];
                float centr = c->cBinMin + (c->cBinMax - c->cBinMin)/2;
                float cRatio = 0;
                if (c->nEventsV0M>0)
                    cRatio = (float)c->nEventsV0M_and_ZDCZEM / c->nEventsV0M;
                grFractEstByV0M[cW]->SetPoint(grFractEstByV0M[cW]->GetN(), centr, cRatio);
                if (c->nEventsZDCZEM>0)
                    cRatio = (float)c->nEventsV0M_and_ZDCZEM / c->nEventsZDCZEM;
                grFractEstByZDC[cW]->SetPoint(grFractEstByZDC[cW]->GetN(), centr, cRatio);
            }
            for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                {
                    WinPair *w = &wins[cW][cBin][etaW][phiW];
                    float centr = w->cBinMin + (w->cBinMax - w->cBinMin)/2;
                    w->calcCorrCoeffs();
                    if(0)cout << "cMin=" << w->cBinMin << ", cMax=" << w->cBinMax << ", etaW=" << etaW
                              << ", NN_bCorr= " << w->NN_bCorr
                              << ", PtPt_bCorr= " << w->PtPt_bCorr
                              << endl;

                    //fill graphs
                    TGraphErrors *gr;
                    //gr NN
                    gr = grNN[cW][etaW];
                    if ( fabs(w->NN_bCorr) < 10 )
                        gr->SetPoint(gr->GetN(), centr, w->NN_bCorr);
                    //gr PtPt
                    gr = grPtPt[cW][etaW];
                    if ( fabs(w->PtPt_bCorr) < 10 )
                        gr->SetPoint(gr->GetN(), centr, w->PtPt_bCorr);
                    //gr PtN
                    gr = grPtN[cW][etaW];
                    if ( fabs(w->PtN_bCorr) < 10 )
                        gr->SetPoint(gr->GetN(), centr, w->PtN_bCorr);

                }
        }

    int colors[] = { kBlack, kBlack, kBlue, kRed, kMagenta };
    int markers[] = { 20, 24, 5, 2, 21 };

    // NN
    TCanvas *canv_grNN = new TCanvas("canv_grNN","canv_grNN",20,50,700,600 );
    tuneCanvas(canv_grNN);
    grNN[0][0]->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels(grNN[0][0]);


    //centr 10
    drawGraph(grNN[0][0], 20, kBlack, "AP");
    //    drawGraph(grNN[0][1], 21, kBlack, "P");
    //    drawGraph(grNN[0][2], 22, kBlack, "P");

    for ( int cW = 1; cW < nCW; cW++ )
        drawGraph(grNN[cW][0], markers[cW], colors[cW], "P");

    grNN[0][0]->SetMinimum( 0 );

    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grNN[cW][0], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();

    TLatex *tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);

    tex = new TLatex(0.4,0.89, "NN");
    drawTex(tex, 0.045);


    TString strPostfix;

    if (flag_V0M_ZDC==0)
        strPostfix = Form("V0M.eps");
    else
        strPostfix = Form("ZDCZEM.eps");

    canv_grNN->SaveAs( Form("NN_%s", strPostfix.Data() ) );

    // PtPt
    TCanvas *canv_grPtPt = new TCanvas("canv_grPtPt","canv_grPtPt",250,50,700,600 );
    tuneCanvas(canv_grPtPt);
    grPtPt[0][0]->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels(grPtPt[0][0]);


    //centr 10
    drawGraph(grPtPt[0][0], 20, kBlack, "AP");
    //    drawGraph(grPtPt[0][1], 21, kBlack, "P");
    //    drawGraph(grPtPt[0][2], 22, kBlack, "P");

    for ( int cW = 1; cW < nCW; cW++ )
        drawGraph(grPtPt[cW][0], markers[cW], colors[cW], "P");

    grPtPt[0][0]->SetMinimum( 0 );

    leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grPtPt[cW][0], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();

    tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);

    tex = new TLatex(0.4,0.89, "PtPt");
    drawTex(tex, 0.045);


    canv_grNN->SaveAs( Form("PtPt_%s", strPostfix.Data() ) );

    // PtN
    TCanvas *canv_grPtN = new TCanvas("canv_grPtN","canv_grPtN",450,50,700,600 );
    tuneCanvas(canv_grPtN);
    grPtN[0][0]->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels(grPtN[0][0]);


    //centr 10
    drawGraph(grPtN[0][0], 20, kBlack, "AP");
    //    drawGraph(grPtN[0][1], 21, kBlack, "P");
    //    drawGraph(grPtN[0][2], 22, kBlack, "P");

    for ( int cW = 1; cW < nCW; cW++ )
        drawGraph(grPtN[cW][0], markers[cW], colors[cW], "P");

    grPtN[0][0]->SetMinimum( 0 );
    leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grPtN[cW][0], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();

    tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);

    tex = new TLatex(0.4,0.89, "PtN");
    drawTex(tex, 0.045);


    canv_grNN->SaveAs( Form("PtN_%s", strPostfix.Data() ) );



    TCanvas *canv_grPtN_2D = new TCanvas("canv_grPtN_2D","canv_grPtN_2D",450,50,700,600 );
    tuneCanvas(canv_grPtN_2D);


    //    wins[0][0][0][0].hist2D_PtN->DrawCopy();
    wins[0][0][0][0].hist2D_PtN->ProfileX()->DrawCopy();



    // CENTR ESTIMATOR EVENT RATIO:
    TCanvas *canv_grCentrRatio = new TCanvas("canv_grCentrRatio","canv_grCentrRatio",450,250,700,600 );
    tuneCanvas(canv_grCentrRatio);
    grFractEstByV0M[0]->SetTitle(";centrality percentile;ratio");
    tuneGraphAxisLabels(grFractEstByV0M[0]);
    //centr 10
    drawGraph(grFractEstByV0M[0], 20, kBlack, "AP");
    //    drawGraph(grFractEstByZDC[0], 20, kBlack, "L");

    for ( int cW = 1; cW < nCW; cW++ )
        drawGraph(grFractEstByV0M[cW], markers[cW], colors[cW], "P");

    leg = new TLegend(0.65,0.65,0.999,0.95,"ratio #frac{V0M-and-ZEMZDC}{V0M}");
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grFractEstByV0M[cW], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();
    //    leg->SetHeader("ratio #frac{V0M-and-ZEMZDC}/{V0M}");

    canv_grCentrRatio->SaveAs("ratio_V0M-and-ZEMZDC_by_V0M.eps");

    // CENTR ESTIMATOR PERCENTILES QA:
    TCanvas *canv_estimatorPercentiles_QA = new TCanvas("canv_estimatorPercentiles_QA","canv_estimatorPercentiles_QA",50,350,700,600 );
    for ( int cBin = 0; cBin < nCentrBins[0]; cBin++ )
    {
        TH1D *h = wins[0][cBin][0][0].hist1D_EstimatorEntries;
        h->SetLineColor(kOrange-9+cBin);
        if ( cBin == 0 )
            h->DrawCopy();
        else
            h->DrawCopy("same");
    }


    // MULT F IN CENTR CLASSES QA:
    TCanvas *canv_mult_F_in_centr_QA = new TCanvas("canv_mult_F_in_centr_QA","canv_mult_F_in_centr_QA",50,400,700,600 );
    for ( int cBin = 0; cBin < nCentrBins[0]; cBin++ )
    {
        TH1D *h = wins[0][cBin][0][0].hist1D_multDistrF;
        h->SetLineColor(kOrange-9+cBin);
        if ( cBin == 0 )
            h->DrawCopy();
        else
            h->DrawCopy("same");
    }



    return;
}



