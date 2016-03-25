//#include <TBranch.h>
//#include <TFile.h>
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"


#include <iostream>
#include <fstream>

using namespace std;


struct BranchFB
{
    UShort_t nF;
    UShort_t nB;
    Float_t PtF;
    Float_t PtB;
    BranchFB() : nF(0), nB(0), PtF(0), PtB(0)
    {}
};

struct CorrCoeffInfo
{
    double avF;
    double avB;
    double avFB;
    double avF2;

    double DFB;
    double DF2;

    double bCorr;
    double C2;

    CorrCoeffInfo():
      avF(-1000)
    , avB(-1000)
    , avFB(-1000)
    , avF2(-1000)
    , DFB(-1000)
    , DF2(-1000)
    {}
};

struct BootstrapHistos
{
    TH1D *hist_avF;
    TH1D *hist_avB;
    TH1D *hist_avFB;
    TH1D *hist_avF2;

    TH1D *hist_DFB;
    TH1D *hist_DF2;

    TH1D *hist_bCorr;
    TH1D *hist_C2;

    BootstrapHistos();
    void InitHistos(const char *strType, float _cMin, float _cMax, int _etaW, int _phiW , double histRangeFB = 1. );
    void FillHistos(const CorrCoeffInfo &corrInfo );
    void WriteHistos();
private:
    void initWithNameTitle( TH1D** hist, const char *strType, const char *strName, float _cMin, float _cMax, int _etaW, int _phiW, double histRange = 1. );
};

class WinPair
{
public:
    float cBinMin;
    float cBinMax;
    int etaW;
    int phiW;

    double NN_nF;
    double NN_nB;
    double NN_nF_nB;
    double NN_nF2;
    Long64_t NN_Nevents;
    TH2D *hist2D_NN;

    double PtPt_PtF;
    double PtPt_PtB;
    double PtPt_PtF_PtB;
    double PtPt_PtF2;
    Long64_t PtPt_Nevents;
    TH2D *hist2D_PtPt;

    double PtN_nF;
    double PtN_PtB;
    double PtN_nF_PtB;
    double PtN_nF2;
    Long64_t PtN_Nevents;
    TH2D *hist2D_PtN;

    //corr coeffs
//    double NN_bCorr;
//    double NN_C2;
//    double PtPt_bCorr;
//    double PtPt_C2;
//    double PtN_bCorr;
//    double PtN_C2;
    CorrCoeffInfo corrInfo_NN;
    CorrCoeffInfo corrInfo_PtPt;
    CorrCoeffInfo corrInfo_PtN;

    TH1D *hist1D_EstimatorEntries;
    TH1D *hist1D_multDistrF;
    TH1D *hist1D_multDistrB;

    TH1D *hist1D_QA_PtF;
    TH1D *hist1D_QA_PtB;

    // bootstrap data:
    bool doBootstrap;
    double *BS_nF;
    double *BS_nB;
    double *BS_PtF;
    double *BS_PtB;
//    TH1D *hist1D_bCorr_BS_NN;
//    TH1D *hist1D_bCorr_BS_PtPt;
    BootstrapHistos histos_BS_NN;
    BootstrapHistos histos_BS_PtPt;

    //run-by-run histos
    int nRuns;
    TH1D **hist1D_multDistr_RunByRun_F;//[nTrees][nEtaBr];
    TH1D **hist1D_multDistr_RunByRun_B;//[nTrees][nEtaBr];
    TH1D **hist1D_avPtDistr_RunByRun_F;//[nTrees][nEtaBr];
    TH1D **hist1D_avPtDistr_RunByRun_B;//[nTrees][nEtaBr];

    WinPair();

    void init(float _cMin, float _cMax, int _etaW, int _phiW, int MAX_N_EVENTS_FOR_BOOTSTRAP = 0 );
    void initRunByRunHistos(int _nRuns, const int *runListNumbers);
    void writeHistos();
    void fill(Float_t cPerc, UShort_t nF, UShort_t nB, Float_t PtF, Float_t PtB, int treeId = -1 );
    void calcCorrCoeffs();
    void performBootstrapping(int corrType); //0-NN, 1-PtPt, 2-PtN

private:
    CorrCoeffInfo _calc( const double &F, const double &B, const double &FB
                 , const double &F2, const int &nEvents, const int &ifRel=0);
};

//class WinQA
//{
//    TH1D *hist1D_multInWin;

//};

struct CentralityOccupancy
{
public:
    float cBinMin;
    float cBinMax;

    int nEventsV0M;
    int nEventsZDCZEM;
    int nEventsV0M_and_ZDCZEM;

    CentralityOccupancy();
    void fill( Float_t cPercV0M, Float_t cPercZDCvsZNA );
};


struct GraphsCorrInfo
{
    // separate class, don't need others except for CorrCoeffInfo (which is needed to SetPoints)
    TGraphErrors *gr_bCorr;
    TGraphErrors *gr_C2;

    TGraphErrors *gr_DFB;
    TGraphErrors *gr_DF2;
    TGraphErrors *gr_avFB;
    TGraphErrors *gr_avF2;
    TGraphErrors *gr_avF;
    TGraphErrors *gr_avB;

    GraphsCorrInfo()
    {
        gr_bCorr = new TGraphErrors;
        gr_C2 = new TGraphErrors;

        gr_DFB = new TGraphErrors;
        gr_DF2 = new TGraphErrors;
        gr_avFB = new TGraphErrors;
        gr_avF2 = new TGraphErrors;
        gr_avF = new TGraphErrors;
        gr_avB = new TGraphErrors;
    }

    void SetNames( const char *strType, int cW, int etaW, int phiW )
    {
        gr_bCorr->SetName( Form("gr_%s_bCorr_cW%d_etaW%d_phiW%d", strType, cW, etaW, phiW) );
        gr_C2->SetName( Form("gr_%s_C2_cW%d_etaW%d_phiW%d", strType, cW, etaW, phiW) );

        gr_DFB ->SetName( Form("gr_%s_DFB_cW%d_etaW%d_phiW%d", strType, cW, etaW, phiW) );
        gr_DF2 ->SetName( Form("gr_%s_DF2_cW%d_etaW%d_phiW%d", strType, cW, etaW, phiW) );
        gr_avFB->SetName( Form("gr_%s_avFB_cW%d_etaW%d_phiW%d", strType, cW, etaW, phiW) );
        gr_avF2->SetName( Form("gr_%s_avF2_cW%d_etaW%d_phiW%d", strType, cW, etaW, phiW) );
        gr_avF ->SetName( Form("gr_%s_avF_cW%d_etaW%d_phiW%d", strType, cW, etaW, phiW) );
        gr_avB ->SetName( Form("gr_%s_avB_cW%d_etaW%d_phiW%d", strType, cW, etaW, phiW) );
    }
    void SetPoints( double xValue, const CorrCoeffInfo *corrInfo )
    {
        int point = gr_bCorr->GetN();

        gr_bCorr->SetPoint( point, xValue, corrInfo->bCorr );
        gr_C2   ->SetPoint( point, xValue, corrInfo->C2 );

        gr_DFB -> SetPoint( point, xValue, corrInfo->DFB );
        gr_DF2 -> SetPoint( point, xValue, corrInfo->DF2 );
        gr_avFB-> SetPoint( point, xValue, corrInfo->avFB );
        gr_avF2-> SetPoint( point, xValue, corrInfo->avF2 );
        gr_avF -> SetPoint( point, xValue, corrInfo->avF );
        gr_avB -> SetPoint( point, xValue, corrInfo->avB );
    }

    void SetPointsWithErrorsFromHistos( double xValue, const BootstrapHistos &BS_histos, bool copyFromFirst_cBin = false )
    {
        int point = gr_bCorr->GetN();

        //points
        gr_bCorr->SetPoint( point, xValue, BS_histos.hist_bCorr->GetMean() );
        gr_C2   ->SetPoint( point, xValue, BS_histos.hist_C2->GetMean() );

        gr_DFB -> SetPoint( point, xValue, BS_histos.hist_DFB ->GetMean() );
        gr_DF2 -> SetPoint( point, xValue, BS_histos.hist_DF2 ->GetMean() );
        gr_avFB-> SetPoint( point, xValue, BS_histos.hist_avFB->GetMean() );
        gr_avF2-> SetPoint( point, xValue, BS_histos.hist_avF2->GetMean() );
        gr_avF -> SetPoint( point, xValue, BS_histos.hist_avF ->GetMean() );
        gr_avB -> SetPoint( point, xValue, BS_histos.hist_avB ->GetMean() );

        //errors
        if ( !copyFromFirst_cBin || ( copyFromFirst_cBin && point == 0 ) )
        {
            // in case when asked for boostrap in all cBins OR for the first cBin
            gr_bCorr->SetPointError( point, 0, BS_histos.hist_bCorr->GetRMS() );
            gr_C2   ->SetPointError( point, 0, BS_histos.hist_C2->GetRMS() );

            gr_DFB -> SetPointError( point, 0, BS_histos.hist_DFB ->GetRMS() );
            gr_DF2 -> SetPointError( point, 0, BS_histos.hist_DF2 ->GetRMS() );
            gr_avFB-> SetPointError( point, 0, BS_histos.hist_avFB->GetRMS() );
            gr_avF2-> SetPointError( point, 0, BS_histos.hist_avF2->GetRMS() );
            gr_avF -> SetPointError( point, 0, BS_histos.hist_avF ->GetRMS() );
            gr_avB -> SetPointError( point, 0, BS_histos.hist_avB ->GetRMS() );
        }
        else //copy errors from cBin=0 point
        {
            gr_bCorr->SetPointError( point, 0, gr_bCorr->GetErrorY( 0 ) );
            gr_C2   ->SetPointError( point, 0, gr_C2   ->GetErrorY( 0 ) );

            gr_DFB -> SetPointError( point, 0, gr_DFB -> GetErrorY( 0 ) );
            gr_DF2 -> SetPointError( point, 0, gr_DF2 -> GetErrorY( 0 ) );
            gr_avFB-> SetPointError( point, 0, gr_avFB-> GetErrorY( 0 ) );
            gr_avF2-> SetPointError( point, 0, gr_avF2-> GetErrorY( 0 ) );
            gr_avF -> SetPointError( point, 0, gr_avF -> GetErrorY( 0 ) );
            gr_avB -> SetPointError( point, 0, gr_avB -> GetErrorY( 0 ) );
        }
    }

    void WriteGraphs();
};
