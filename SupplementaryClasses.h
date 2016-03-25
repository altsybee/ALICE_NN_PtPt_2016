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
    double NN_bCorr;
    double NN_C2;
    double PtPt_bCorr;
    double PtPt_C2;
    double PtN_bCorr;
    double PtN_C2;

    TH1D *hist1D_EstimatorEntries;
    TH1D *hist1D_multDistrF;
    TH1D *hist1D_multDistrB;

    TH1D *hist1D_QA_PtF;
    TH1D *hist1D_QA_PtB;

    bool doBootstrap;
    double *BS_nF;
    double *BS_nB;
    double *BS_PtF;
    double *BS_PtB;
    TH1D *hist1D_bCorr_BS_NN;
    TH1D *hist1D_bCorr_BS_PtPt;

    //run-by-run histos
    TH1D **hist1D_multDistr_RunByRun_F;//[nTrees][nEtaBr];
    TH1D **hist1D_multDistr_RunByRun_B;//[nTrees][nEtaBr];
    TH1D **hist1D_avPtDistr_RunByRun_F;//[nTrees][nEtaBr];
    TH1D **hist1D_avPtDistr_RunByRun_B;//[nTrees][nEtaBr];

    WinPair();

    void init(float _cMin, float _cMax, int _etaW, int _phiW, int MAX_N_EVENTS_FOR_BOOTSTRAP = 0 );
    void initRunByRunHistos(int nRuns, const int *runListNumbers);
    void fill(Float_t cPerc, UShort_t nF, UShort_t nB, Float_t PtF, Float_t PtB, int treeId = -1 );
    void calcCorrCoeffs();
    void performBootstrapping(int corrType); //0-NN, 1-PtPt, 2-PtN

private:
    double _calc(int type, const double &F, const double &B, const double &FB
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
