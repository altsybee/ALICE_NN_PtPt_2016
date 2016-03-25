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

#include "TFile.h"
#include "TString.h"

#include "SupplementaryClasses.h"

#include <iostream>
#include <fstream>

using namespace std;


BootstrapHistos::BootstrapHistos():
  hist_avF(0)
, hist_avB(0)
, hist_avFB(0)
, hist_avF2(0)
, hist_DFB(0)
, hist_DF2(0)
{}

void BootstrapHistos::InitHistos( const char *strType, float _cMin, float _cMax, int _etaW, int _phiW, double histRangeFB )
{
    cout << "InitHistos for BS - " << strType << endl;
    initWithNameTitle( &hist_avF  , strType, "avF"    , _cMin, _cMax, _etaW, _phiW, histRangeFB );
    initWithNameTitle( &hist_avB  , strType, "avB"    , _cMin, _cMax, _etaW, _phiW, histRangeFB );
    initWithNameTitle( &hist_avFB , strType, "avFB"   , _cMin, _cMax, _etaW, _phiW, histRangeFB*histRangeFB );
    initWithNameTitle( &hist_avF2 , strType, "avF2"   , _cMin, _cMax, _etaW, _phiW, histRangeFB*histRangeFB );
    initWithNameTitle( &hist_DFB  , strType, "DFB"    , _cMin, _cMax, _etaW, _phiW, histRangeFB*10 );
    initWithNameTitle( &hist_DF2  , strType, "DF2"    , _cMin, _cMax, _etaW, _phiW, histRangeFB*10 );
    initWithNameTitle( &hist_bCorr, strType, "bCorr"  , _cMin, _cMax, _etaW, _phiW );
    initWithNameTitle( &hist_C2   , strType, "C2"     , _cMin, _cMax, _etaW, _phiW );
}

void BootstrapHistos::initWithNameTitle( TH1D** hist, const char *strType, const char *strName, float _cMin, float _cMax, int _etaW, int _phiW, double histRange )
{
    TString strBSname  = Form("hist1D_BS_%s_%s_c%.1f-%.1f_etaW_%d_phiW_%d", strType, strName, _cMin, _cMax, _etaW, _phiW);
    TString strBStitle = Form("hist1D_BS_%s_%s_c%.1f-%.1f_etaW_%d_phiW_%d;bCorr;entries", strType, strName, _cMin, _cMax, _etaW, _phiW);
    *hist = new TH1D( strBSname, strBStitle, 4000, -histRange, histRange );
}

void BootstrapHistos::FillHistos( const CorrCoeffInfo &corrInfo )
{
//    cout << "FillHistos for BS" << endl;
    hist_avF  ->Fill( corrInfo.avF   );
    hist_avB  ->Fill( corrInfo.avB   );
    hist_avFB ->Fill( corrInfo.avFB  );
    hist_avF2 ->Fill( corrInfo.avF2  );
    hist_DFB  ->Fill( corrInfo.DFB   );
    hist_DF2  ->Fill( corrInfo.DF2   );
    hist_bCorr->Fill( corrInfo.bCorr );
    hist_C2   ->Fill( corrInfo.C2    );
}

void BootstrapHistos::WriteHistos()
{
    TString strDirName = Form( "dir_%s" , hist_bCorr->GetName() );
    gFile->mkdir( strDirName.Data() );
    gFile->cd( strDirName.Data() );

    hist_bCorr->Write();
    hist_C2   ->Write();

    hist_DFB -> Write();
    hist_DF2 -> Write();
    hist_avFB-> Write();
    hist_avF2-> Write();
    hist_avF -> Write();
    hist_avB -> Write();
    gFile->cd();
}


// ##### WinPair
WinPair::WinPair() :
    cBinMin(0)
  , cBinMax(100)
  , etaW(-1)
  , phiW(-1)
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
  , nRuns(-1)

//  , NN_bCorr(-1000)
//  , NN_C2(-1000)
//  , PtPt_bCorr(-1000)
//  , PtPt_C2(-1000)
//  , PtN_bCorr(-1000)
//  , PtN_C2(-1000)

  , hist2D_NN                (0)
  , hist2D_PtPt              (0)
  , hist2D_PtN               (0)
  , hist1D_EstimatorEntries  (0)
  , hist1D_multDistrF        (0)
  , hist1D_multDistrB        (0)
  , hist1D_QA_PtF            (0)
  , hist1D_QA_PtB            (0)

  , doBootstrap(false)
  , BS_PtF(0)
  , BS_PtB(0)

  , hist1D_multDistr_RunByRun_F(0x0)
  , hist1D_multDistr_RunByRun_B(0x0)
  , hist1D_avPtDistr_RunByRun_F(0x0)
  , hist1D_avPtDistr_RunByRun_B(0x0)

{}
void WinPair::init(float _cMin, float _cMax, int _etaW, int _phiW, int MAX_N_EVENTS_FOR_BOOTSTRAP )
{
    cBinMin = _cMin;
    cBinMax = _cMax;
    etaW = _etaW;
    phiW = _phiW;

    TString strNN = Form("hist2D_NN_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    hist2D_NN = new TH2D( strNN, strNN, 800, -0.5, 799.5, 800, -0.5, 799.5 );

    TString strPtPt = Form("hist2D_PtPt_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    hist2D_PtPt = new TH2D( strPtPt, strPtPt, 400, 0, 2, 400, 0, 2);

    TString strPtN = Form("hist2D_PtN_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    hist2D_PtN = new TH2D( strPtN, strPtN, 1000, -0.5, 999.5, 400, 0, 2);

    //QA histos:
    TString strEstPerc = Form("hist1D_EstimatorEntries_c%.1f-%.1f_etaW_%d_phiW_%d;percentile;entries", cBinMin, cBinMax, etaW, phiW);
    hist1D_EstimatorEntries = new TH1D( strEstPerc, strEstPerc, 20001, -0.5, 1000.5);

    TString strMultDistrF_name  = Form("hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    TString strMultDistrF_title = Form("hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d;n tracks;n events", cBinMin, cBinMax, etaW, phiW);
    hist1D_multDistrF = new TH1D( strMultDistrF_name, strMultDistrF_title, 3001, -0.5, 3000.5);

    TString strMultDistrB_name  = Form("hist1D_multDistrB_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    TString strMultDistrB_title = Form("hist1D_multDistrB_c%.1f-%.1f_etaW_%d_phiW_%d;n tracks;n events", cBinMin, cBinMax, etaW, phiW);
    hist1D_multDistrB = new TH1D( strMultDistrB_name, strMultDistrB_title, 3001, -0.5, 3000.5);

    TString strPtF_name = Form("hist1D_PtF_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    TString strPtF_title = Form("hist1D_PtF_c%.1f-%.1f_etaW_%d_phiW_%d;#LTp_{T}#GT Forward;n events", cBinMin, cBinMax, etaW, phiW);
    hist1D_QA_PtF = new TH1D( strPtF_name, strPtF_title, 2002, -2, 2 );

    TString strPtB_name = Form("hist1D_PtB_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    TString strPtB_title = Form("hist1D_PtB_c%.1f-%.1f_etaW_%d_phiW_%d;#LTp_{T}#GT Backward;n events", cBinMin, cBinMax, etaW, phiW);
    hist1D_QA_PtB = new TH1D( strPtB_name, strPtB_title, 2002, -2, 2 );

    if ( MAX_N_EVENTS_FOR_BOOTSTRAP > 0 )
    {
        doBootstrap = true;
        BS_nF = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        BS_nB = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        BS_PtF = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        BS_PtB = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];

        histos_BS_NN.InitHistos( "NN", cBinMin, cBinMax, etaW, phiW, 1000 );
        histos_BS_PtPt.InitHistos( "PtPt", cBinMin, cBinMax, etaW, phiW );

//        TString strBS_NN_name  = Form("hist1D_bCorr_BS_NN_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
//        TString strBS_NN_title = Form("hist1D_bCorr_BS_NN_c%.1f-%.1f_etaW_%d_phiW_%d;bCorr;entries", cBinMin, cBinMax, etaW, phiW);
//        hist1D_bCorr_BS_NN = new TH1D( strBS_NN_name, strBS_NN_title, 2000, -1, 1 );

//        TString strBS_PtPt_name  = Form("hist1D_bCorr_BS_PtPt_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
//        TString strBS_PtPt_title = Form("hist1D_bCorr_BS_PtPt_c%.1f-%.1f_etaW_%d_phiW_%d;bCorr;entries", cBinMin, cBinMax, etaW, phiW);
//        hist1D_bCorr_BS_PtPt = new TH1D( strBS_PtPt_name, strBS_PtPt_title, 2000, -1, 1 );

    }

}

void WinPair::initRunByRunHistos( int _nRuns, const int *runListNumbers )
{
    nRuns = _nRuns;
    hist1D_multDistr_RunByRun_F = new TH1D*[nRuns];
    hist1D_multDistr_RunByRun_B = new TH1D*[nRuns];
    hist1D_avPtDistr_RunByRun_F = new TH1D*[nRuns];
    hist1D_avPtDistr_RunByRun_B = new TH1D*[nRuns];

    TString str_hist_name;
    for ( int r = 0; r < nRuns; r++ )
    {
        int runNumber = runListNumbers[r];
        TString str_hist_title  = Form("%d", runNumber);

        //mult
        str_hist_name  = Form("hist1D_multDistr_winF_run_%d_c%.1f-%.1f_etaW_%d_phiW_%d", runNumber, cBinMin, cBinMax, etaW, phiW);
        hist1D_multDistr_RunByRun_F[r] = new TH1D( str_hist_name, ";n tracks;n events", 1000, -0.5, 1000-0.5 );

        str_hist_name  = Form("hist1D_multDistr_winB_run_%d_c%.1f-%.1f_etaW_%d_phiW_%d", runNumber, cBinMin, cBinMax, etaW, phiW);
        hist1D_multDistr_RunByRun_B[r] = new TH1D( str_hist_name, ";n tracks;n events", 1000, -0.5, 1000-0.5 );

        hist1D_multDistr_RunByRun_F[r]->SetTitle( str_hist_title );
        hist1D_multDistr_RunByRun_B[r]->SetTitle( str_hist_title );

        //av pT
        str_hist_name  = Form("hist1D_avPtDistr_winF_run_%d_c%.1f-%.1f_etaW_%d_phiW_%d", runNumber, cBinMin, cBinMax, etaW, phiW);
        hist1D_avPtDistr_RunByRun_F[r] = new TH1D( str_hist_name, ";#LTp_{T}#GT;n events", 1000, 0, 2 );

        str_hist_name  = Form("hist1D_avPtDistr_winB_run_%d_c%.1f-%.1f_etaW_%d_phiW_%d", runNumber, cBinMin, cBinMax, etaW, phiW);
        hist1D_avPtDistr_RunByRun_B[r] = new TH1D( str_hist_name, ";#LTp_{T}#GT;n events", 1000, 0, 2 );

        hist1D_avPtDistr_RunByRun_F[r]->SetTitle( str_hist_title );
        hist1D_avPtDistr_RunByRun_B[r]->SetTitle( str_hist_title );
    }

}


void WinPair::fill(Float_t cPerc, UShort_t nF, UShort_t nB, Float_t PtF, Float_t PtB, int treeId )
{
    //check if we in centrality bin
    if ( cPerc < cBinMin || cPerc > cBinMax )
        return;

    //QA fills:
    hist1D_EstimatorEntries->Fill(cPerc);
    hist1D_multDistrF->Fill(nF);
    hist1D_multDistrB->Fill(nB);

    //QA run-by-run
    if ( hist1D_multDistr_RunByRun_F != 0x0 )
    {
        //F
        hist1D_multDistr_RunByRun_F[treeId]->Fill( nF);
        if ( nF > 0 )
            hist1D_avPtDistr_RunByRun_F[treeId]->Fill( PtF );

        //B
        hist1D_multDistr_RunByRun_B[treeId]->Fill( nB);
        if ( nB > 0 )
            hist1D_avPtDistr_RunByRun_B[treeId]->Fill( PtB );
    }


    //NN
    hist2D_NN->Fill( nF, nB );

    NN_nF       += nF;
    NN_nB       += nB;
    NN_nF_nB    += nF*nB;
    NN_nF2      += nF*nF;

    if ( doBootstrap )
    {
        BS_nF[NN_Nevents] = nF;
        BS_nB[NN_Nevents] = nB;
    }

    NN_Nevents++;

    hist1D_QA_PtF->Fill( PtF );
    hist1D_QA_PtB->Fill( PtB );

    //PtPt
    if ( nF > 0 && nB > 0 )
    {
        hist2D_PtPt->Fill( PtF, PtB );

        PtPt_PtF += PtF;
        PtPt_PtB += PtB;
        PtPt_PtF_PtB += PtF*PtB;
        PtPt_PtF2 += PtF*PtF;

        if ( doBootstrap )
        {
            BS_PtF[PtPt_Nevents] = PtF;
            BS_PtB[PtPt_Nevents] = PtB;
        }

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

void WinPair::calcCorrCoeffs()
{
    corrInfo_NN    = _calc( NN_nF, NN_nB, NN_nF_nB, NN_nF2, NN_Nevents );
//    NN_C2       = _calc( NN_nF, NN_nB, NN_nF_nB, NN_nF2, NN_Nevents );

    corrInfo_PtPt  = _calc( PtPt_PtF, PtPt_PtB, PtPt_PtF_PtB, PtPt_PtF2, PtPt_Nevents );
//    PtPt_C2       = _calc( PtPt_PtF, PtPt_PtB, PtPt_PtF_PtB, PtPt_PtF2, PtPt_Nevents );

    corrInfo_PtN   = _calc( PtN_nF, PtN_PtB, PtN_nF_PtB, PtN_nF2, PtN_Nevents, 1 );
//    PtN_C2       = _calc( PtN_nF, PtN_PtB, PtN_nF_PtB, PtN_nF2, PtN_Nevents );
}

void WinPair::writeHistos()
{
    hist2D_NN->Write();
    hist2D_PtPt->Write();
    hist2D_PtN->Write();

    hist1D_multDistrF->Write();
    hist1D_multDistrB->Write();

    hist1D_QA_PtF->Write();
    hist1D_QA_PtB->Write();

    if ( nRuns > 0 ) //we have run-by-run histos
    {
        for ( int treeId = 0; treeId < nRuns; treeId++ )
        {
            hist1D_multDistr_RunByRun_F[treeId]->Write();
            hist1D_multDistr_RunByRun_B[treeId]->Write();
            hist1D_avPtDistr_RunByRun_F[treeId]->Write();
            hist1D_avPtDistr_RunByRun_B[treeId]->Write();
        }
    }
    //                    wins[cW][cBin][etaW][phiW].hist2D_NN->ProfileX()->Write();
    //                    wins[cW][cBin][etaW][phiW].hist2D_PtPt->ProfileX()->Write();
    //                    wins[cW][cBin][etaW][phiW].hist2D_PtN->ProfileX()->Write();

}

//double WinPair::_calc(int type, double F, double B, double FB, double F2, int nEvents, int ifRel )
CorrCoeffInfo WinPair::_calc( const double &F, const double &B, const double &FB
                       , const double &F2, const int &nEvents, const int &ifRel )
{
    CorrCoeffInfo corrInfo;

    if ( nEvents <= 0 )
        return corrInfo;
    double meanF     =  F   / nEvents;     //printf( "<F> = %f\n", meanF );
    double meanB     =  B   / nEvents;     //printf( "<B> = %f\n", meanB );
    double meanFB    =  FB  / nEvents;    //printf( "<FB> = %f\n", meanFB );
    double meanF2    =  F2  / nEvents;    //printf( "<FB> = %f\n", meanFB );

    //        printf( "nEvents = %d, ", nEvents );
    //        printf( "meanF = %f, ", meanF  );
    //        printf( "meanB = %f, ", meanB  );
    //        printf( "meanFB = %f, ",  meanFB );
    //        printf( "meanF2 = %f\n, ",  meanF2 );

    double numerator = meanFB - meanF * meanB;
    double denominator_bCorr = meanF2 - meanF*meanF;
    double denominator_C2 = meanF * meanB;
//    if (type == 0) // bCorr
//        denominator = meanF2 - meanF*meanF;
//    else if (type == 1) // C2
//        denominator = meanF * meanB;

    //bcorr
    double bcorr = -1000;
    if ( denominator_bCorr != 0 )
    {
        bcorr = numerator / denominator_bCorr;
        if (ifRel)
            bcorr *= meanF/meanB;
    }

    //C2
    double C2 = -1000;
    if ( denominator_C2 !=0 )
        C2 = numerator / denominator_C2;

    //fill corrInfo data:
    corrInfo.avF = meanF;
    corrInfo.avB = meanB;
    corrInfo.avFB = meanFB;
    corrInfo.avF2 = meanF2;

    corrInfo.DFB = numerator;
    corrInfo.DF2 = denominator_bCorr;

    corrInfo.bCorr = bcorr;
    corrInfo.C2 = C2;

    return corrInfo;
}



void WinPair::performBootstrapping(int corrType)
{
    for ( int t = 0; t < 1000; t++ )
    {
        int _nDataEvents = 0;

        //corr type for bootstrap:
        if ( corrType == 0 )
            _nDataEvents = NN_Nevents;
        else if ( corrType == 1 )
            _nDataEvents = PtPt_Nevents;
        else
            cout << ">>> WARNING: up to now - bootstrapping for NN and PtPt only!" << endl;
//        else if ( corrType == 2 )
//            _nDataEvents = PtN_Nevents;

        double BS_F = 0;
        double BS_B = 0;
        double BS_FB = 0;
        double BS_F2 = 0;

        Long64_t BS_Nevents = 0;

        double *F;
        double *B;

//        cout << ">>> start bootstrap event loop..." << endl;
        for ( int ev = 0; ev < _nDataEvents; ev++ )
        {
            int id = TMath::Nint( gRandom->Uniform(-0.5, _nDataEvents-0.5) );

            if ( corrType == 0 )
            {
                F = &BS_nF[id];
                B = &BS_nB[id];
            }
            else if ( corrType == 1 )
            {
                F = &BS_PtF[id];
                B = &BS_PtB[id];
//                        cout << ">>> *F = " << *F << ", *B = " << *B << ", ";
            }

            BS_F += *F;
            BS_B += *B;
            BS_FB += (*F)*(*B);
            BS_F2 += (*F)*(*F);

            BS_Nevents++;

//            cout << ">>> BS_F = " << BS_F << ", BS_B = " << BS_B
//                << ">>> BS_FB = " << BS_FB << ", BS_F2 = " << BS_F2
//                << ", BS_Nevents = " << BS_Nevents << endl;
//            int a;
//            cin >> a;
        }
//        cout << ">>> BS_Nevents = " << BS_Nevents << endl;

        CorrCoeffInfo BS_corrInfo = _calc( BS_F, BS_B, BS_FB, BS_F2, BS_Nevents );
//        double BS_bCorr = BS_corrInfo.bCorr;
//                cout << ">>> BS_bCorr = " << BS_bCorr << endl;
        //double BS_C2    = _calc( 1, BS_F, BS_B, BS_FB, BS_F2, BS_Nevents );


        if ( corrType == 0 )
            //hist1D_bCorr_BS_NN->Fill( BS_bCorr );
            histos_BS_NN.FillHistos( BS_corrInfo );
        else if ( corrType == 1 )
//            hist1D_bCorr_BS_PtPt->Fill( BS_bCorr );
            histos_BS_PtPt.FillHistos( BS_corrInfo );

    }
//    cout << ">>> END BS " << endl;
}






// ##### CentralityOccupancy

CentralityOccupancy::CentralityOccupancy() :
    cBinMin(0)
  , cBinMax(100)
  , nEventsV0M(0)
  , nEventsZDCZEM(0)
  , nEventsV0M_and_ZDCZEM(0)
{}
void CentralityOccupancy::fill( Float_t cPercV0M, Float_t cPercZDCvsZNA )
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




void GraphsCorrInfo::WriteGraphs()
{
    TString strDirName = Form( "dir_%s" , gr_bCorr->GetName() );
    gFile->mkdir( strDirName.Data() );
    gFile->cd( strDirName.Data() );

    gr_bCorr->Write();
    gr_C2   ->Write();

    gr_DFB -> Write();
    gr_DF2 -> Write();
    gr_avFB-> Write();
    gr_avF2-> Write();
    gr_avF -> Write();
    gr_avB -> Write();

    gFile->cd();

}
