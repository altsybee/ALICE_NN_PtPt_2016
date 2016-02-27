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

#include "SupplementaryClasses.h"

#include <iostream>
#include <fstream>

using namespace std;



// ##### WinPair


WinPair::WinPair() :
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

{}
void WinPair::init(float _cMin, float _cMax, int etaW, int phiW, int MAX_N_EVENTS_FOR_BOOTSTRAP )
{
    cBinMin = _cMin;
    cBinMax = _cMax;

    TString strNN = Form("hist2D_NN_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    hist2D_NN = new TH2D( strNN, strNN, 800, -0.5, 799.5, 800, -0.5, 799.5 );

    TString strPtPt = Form("hist2D_PtPt_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    hist2D_PtPt = new TH2D( strPtPt, strPtPt, 400, 0, 2, 400, 0, 2);

    TString strPtN = Form("hist2D_PtN_c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    hist2D_PtN = new TH2D( strPtN, strPtN, 1000, -0.5, 999.5, 400, 0, 2);

    //QA histos:
    TString strEstPerc = Form("hist1D_EstimatorEntries_c%.1f-%.1f_etaW_%d_phiW_%d;percentile;entries", cBinMin, cBinMax, etaW, phiW);
    hist1D_EstimatorEntries = new TH1D( strEstPerc, strEstPerc, 20001, -0.5, 1000.5);

    TString strMultDistrF = Form("hist1D_multDistrF_c%.1f-%.1f_etaW_%d_phiW_%d;n tracks;n events", cBinMin, cBinMax, etaW, phiW);
    hist1D_multDistrF = new TH1D( strMultDistrF, strMultDistrF, 3001, -0.5, 3000.5);

    TString strMultDistrB= Form("hist1D_multDistrB_c%.1f-%.1f_etaW_%d_phiW_%d;n tracks;n events", cBinMin, cBinMax, etaW, phiW);
    hist1D_multDistrB = new TH1D( strMultDistrB, strMultDistrB, 3001, -0.5, 3000.5);

    TString strPtF = Form("hist1D_PtF_c%.1f-%.1f_etaW_%d_phiW_%d;#LTp_{T}#GT Forward;n events", cBinMin, cBinMax, etaW, phiW);
    hist1D_QA_PtF = new TH1D( strPtF, strPtF, 2002, -2, 2 );

    TString strPtB= Form("hist1D_PtB_c%.1f-%.1f_etaW_%d_phiW_%d;#LTp_{T}#GT Backward;n events", cBinMin, cBinMax, etaW, phiW);
    hist1D_QA_PtB = new TH1D( strPtB, strPtB, 2002, -2, 2 );

    if ( MAX_N_EVENTS_FOR_BOOTSTRAP > 0 )
    {
        doBootstrap = true;
        BS_nF = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        BS_nB = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        BS_PtF = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        BS_PtB = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];

        TString strBS_NN= Form("hist1D_bCorr_BS_NN_c%.1f-%.1f_etaW_%d_phiW_%d;bCorr;entries", cBinMin, cBinMax, etaW, phiW);
        hist1D_bCorr_BS_NN = new TH1D( strBS_NN, strBS_NN, 2000, -1, 1 );

        TString strBS_PtPt= Form("hist1D_bCorr_BS_PtPt_c%.1f-%.1f_etaW_%d_phiW_%d;bCorr;entries", cBinMin, cBinMax, etaW, phiW);
        hist1D_bCorr_BS_PtPt = new TH1D( strBS_PtPt, strBS_PtPt, 2000, -1, 1 );

    }

}
void WinPair::fill( Float_t cPerc, UShort_t nF, UShort_t nB, Float_t PtF, Float_t PtB )
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
    NN_bCorr    = _calc( 0, NN_nF, NN_nB, NN_nF_nB, NN_nF2, NN_Nevents );
    NN_C2       = _calc( 1, NN_nF, NN_nB, NN_nF_nB, NN_nF2, NN_Nevents );

    PtPt_bCorr    = _calc( 0, PtPt_PtF, PtPt_PtB, PtPt_PtF_PtB, PtPt_PtF2, PtPt_Nevents );
    PtPt_C2       = _calc( 1, PtPt_PtF, PtPt_PtB, PtPt_PtF_PtB, PtPt_PtF2, PtPt_Nevents );

    PtN_bCorr    = _calc( 0, PtN_nF, PtN_PtB, PtN_nF_PtB, PtN_nF2, PtN_Nevents, 1 );
    PtN_C2       = _calc( 1, PtN_nF, PtN_PtB, PtN_nF_PtB, PtN_nF2, PtN_Nevents );
}

//double WinPair::_calc(int type, double F, double B, double FB, double F2, int nEvents, int ifRel )
double WinPair::_calc( int type, const double &F, const double &B, const double &FB
                       , const double &F2, const int &nEvents, const int &ifRel )
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

    double numerator = meanFB - meanF * meanB;
    double denominator = 0;
    if (type == 0) // bCorr
        denominator = meanF2 - meanF*meanF;
    else if (type == 1) // C2
        denominator = meanF * meanB;

    double value = 0;
    if (denominator!=0)
        value = numerator / denominator;

    if (ifRel)
        value *= meanF/meanB;


    return value;
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

        cout << ">>> start bootstrap event loop..." << endl;
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

        double BS_bCorr = _calc( 0, BS_F, BS_B, BS_FB, BS_F2, BS_Nevents );
                cout << ">>> BS_bCorr = " << BS_bCorr << endl;
        //double BS_C2    = _calc( 1, BS_F, BS_B, BS_FB, BS_F2, BS_Nevents );


        if ( corrType == 0 )
            hist1D_bCorr_BS_NN->Fill( BS_bCorr );
        else if ( corrType == 1 )
            hist1D_bCorr_BS_PtPt->Fill( BS_bCorr );
    }
    cout << ">>> END BS " << endl;
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

