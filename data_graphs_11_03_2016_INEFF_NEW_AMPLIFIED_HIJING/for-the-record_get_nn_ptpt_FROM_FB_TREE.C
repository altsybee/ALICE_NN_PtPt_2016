#include <TTree.h>
#include <TBranch.h>
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

//#include "TChain.h"
//#include "TFileCollection.h"

#include "SupplementaryClasses.h"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;

#include "utils.C"


// ########## event selection parameters
const float kVertexZcut = 8;
//    if ( isAnalysing502 )
//        kVertexZcut = 7;
const bool USE_VERTEXZ_CUT = 1;
const int flag_V0M_ZDC_CL1 = 0;

bool checkEventSelection( float &cEstimator
                          , const Bool_t &brIsPileupSPD
                          , const Float_t &brZEMvsZDC
                          , const Float_t &brV0M
                          , const Float_t &brCL1
                          , const Float_t &br_vertexZ
                          , int &nRejectedByPileup
                          )
{
    cEstimator = -1;

    // check isPileupSPD
    if ( brIsPileupSPD == 1 )
    {
        nRejectedByPileup++;
        return false;
    }

    // vertex Z cut!!!
    if ( USE_VERTEXZ_CUT && fabs(br_vertexZ) > kVertexZcut
         //         && fabs(br_vertexZ) > 8
         //         && br_vertexZ < 0 //reject vertZ<0 ! -> keep vertZ>0
         //         && br_vertexZ > 0 //reject vertZ>0 ! -> keep vertZ<0
         )
        return false;


    if ( flag_V0M_ZDC_CL1 == 0 )
    {
        if ( brV0M > 90 ) //V0M cut
            return false;
        cEstimator = brV0M;
    }
    else if ( flag_V0M_ZDC_CL1 == 1 )
    {
        if ( brZEMvsZDC > 50 ) //ZDCvsZEM cut
            return false;
        cEstimator = brZEMvsZDC;
    }
    else if ( flag_V0M_ZDC_CL1 == 2 )
    {
        if ( brCL1 > 90 ) //CL1 cut
            return false;
        cEstimator = brCL1;
    }

    return true;
}



TTree *getTreeFromFile( TFile *myFile )
{
    int listId = 0;//1;
    // !!! listId can be different for 2.76 and 5.02!
    //    if ( isAnalysing502 )
    //        listId = 1;

    TList *listKeys = myFile->GetListOfKeys();
    TString listName = Form("%s",listKeys->At(listId)->GetName() );
    while ( 1 && !listName.Contains("PWGCFLRC") )
    {
        listId++;
        listName = Form("%s",listKeys->At(listId)->GetName() );
    }
    //return;
    //SET listId BY HAND:
    //    listId = 0;
    cout << "going into list: " << listKeys->At(listId)->GetName() << endl;

    //return;
    myFile->cd( listKeys->At(listId)->GetName() );
    TList *listKeys2 = gDirectory->GetListOfKeys();
    TList *myTask = (TList*) gDirectory->Get( listKeys2->At(0)->GetName() );
    TTree *t1 = (TTree*)  myTask->FindObject( "treeFB" );

    return t1;
}

void get_nn_ptpt_FROM_FB_TREE() // TString inputFileName = "MergedOutput.root")
{
    // !!! important for branch V0M below!
    bool isAnalysing502 = 0;
    bool isAMPT_kine = 0;


    const int nFiles = 10;//24+22;
    TFile *myFiles[nFiles];

    // FULL STATISTICS
    // LHC10h
    //    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFplus_all/MergedOutput.root" );
    //    myFiles[1] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFminus_all/MergedOutput.root" );

    //HIJING LHC11a10a_bis RECO level
    //    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_28_PbPb_MCAOD_LHC11a10a_bis_AOD162_TEST_Efficiency_kine_vs_reco_try5/MergedOutput.root" );
//    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_03_06_PbPb_MCAOD_LHC11a10a_bis_AOD162_TEST_Efficiency_kine_vs_reco/MergedOutput.root" );

    //AMPT AMPT_LHC13f3c_StrMelt_ON_rescatON_8_eta_wins KINE level
//    isAMPT_kine = 1;
//    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_29_PbPb_AMPT_LHC13f3c_StrMelt_ON_rescatON_8_eta_wins/merged_results_with_scaled_factor.root" );

    // 11.03.2016: INEFF LHC10h
    Int_t runlist[]=
    {
        //24 runs:
        138275, 138225, 138201, 138197,
        138192, 138190, 137848, 137844,
        137752, 137751, 137724, 137722,
        137718, 137704, 137693, 137692,
        137691, 137686, 137685, 137639,
        137638, 137608, 137595, 137549,

        //22 runs:   //137366 excluded (no such a run?) => 22 runs
        137546, 137544, 137541, 137539,
        137531, 137530, 137443, 137441,
        137440, 137439, 137434, 137432,
        137431, 137430, /*137366,*/ 137243,
        137236, 137235, 137232, 137231,
        137230, 137162, 137161
    };
    for (Int_t i=0; i<nFiles; i++)
        myFiles[i] = new TFile( Form("/Users/macbook/alice/aliceAnalysis/results/task_2016_03_10_PbPb_Data_LHC10h_AOD160_3etaWins_phi1_pt02_20_FB_TREE_INEFF_BY_HAND_HIST3D/AnalysisResults_000%d.root", runlist[i] ) );

    //    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFplus_all/MergedOutput.root" );
    //    myFiles[1] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFminus_all/MergedOutput.root" );



    // INEFF BY HAND
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_25_PbPb_Data_LHC10h_AOD160_3etaWins_phi1_pt02_20_FB_TREE_INEFF_BY_HAND/MergedOutput.root" );



    //LHC10h
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_09_PbPb_Data_2_76TeV_LHC10h_AOD0160_3etaWins_phi1_pt02_20_FB_TREE_blocks123/block1/AnalysisResults.139465.root" );
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_09_PbPb_Data_2_76TeV_LHC10h_AOD0160_3etaWins_phi1_pt02_20_FB_TREE_blocks123/MergedOutput.root" );

    //LHC11h
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoMinus/MergedOutput.root" );
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoPlus/MergedOutput.root" );

    //LHC11h Femto Plus Minus merged
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoPlus/MergedResults/MergedOutput_Femto_Plus_Minus.root" );

    //LHC15o MF ++ --
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_22_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldMM/MergedOutput.root" );
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_21_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldPP/MergedOutput.root" );

    //LHC15o MF merged
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_22_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldMM/MergedResults/MergedOutput_PP_and_MM.root" );



    const int nTrees = nFiles;
    TTree *trees[nTrees];
    int nEvInTrees[nTrees];

    //get trees from files
    int nEvents = 0;
    for ( int i = 0; i < nTrees; i++ )
    {
        if ( !myFiles[i] )
        {
            cout << "No input file! i=" << i << endl;
            return;
        }
        myFiles[i]->ls();

        if (!isAMPT_kine)
            trees[i] = (TTree*) getTreeFromFile( myFiles[i] );
        else
            trees[i] = (TTree*) myFiles[i]->Get( "t1" );


        nEvInTrees[i] = trees[i]->GetEntries();
        nEvents += nEvInTrees[i];
    }

    //    nEvents = 5000000;
    //    nEvents = trees[0]->GetEntries() + trees[1]->GetEntries();
//    nEvents = trees[0]->GetEntries();
    cout <<"nEvents = " << nEvents << endl;




    // ########## Branches for event-info
    Bool_t brIsPileupSPD = 0;
    Float_t brZEMvsZDC = 0;
    Float_t brV0M = 0;
    Float_t brCL1 = 0;
    Float_t br_vertexZ = 0;




    // ########## Number of eta-phi wins
    //    const int nEtaBr = 3;
    const int nEtaBr = 8;
    const int nPhiWins = 1;
    bool FLAG_PERCOLATING_WINS = false;//true;

    BranchFB br[nEtaBr][nPhiWins];

    for ( int iTree = 0; iTree < nTrees; iTree++ )
    {
        TTree *tr = trees[iTree];

        tr->SetBranchAddress( "isPileupSPD", &brIsPileupSPD );
        tr->SetBranchAddress( "vertexZ", &br_vertexZ );
        tr->SetBranchAddress( "centr_ZEMvsZDC", &brZEMvsZDC );
        tr->SetBranchAddress( "centr_CL1", &brCL1 );

        if ( !isAnalysing502 )
            tr->SetBranchAddress( "centr_V0M", &brV0M );
        else
            tr->SetBranchAddress( "centrV0M_NEW_MULT_SEL", &brV0M );

        // binding of wins branches:
        for ( int etaW = 0; etaW < nEtaBr; etaW++ )
            for ( int phiW = 0; phiW < nPhiWins; phiW++ )
            {
                int ptW = 0;
                TString brNamePostfix = Form("eta_%d_phi%d_pt_%d"
                                             , etaW, phiW, ptW );
                tr->SetBranchAddress( Form("nF_%s", brNamePostfix.Data() ),  &br[etaW][phiW].nF );
                tr->SetBranchAddress( Form("nB_%s", brNamePostfix.Data() ),  &br[etaW][phiW].nB );
                tr->SetBranchAddress( Form("PtF_%s", brNamePostfix.Data() ), &br[etaW][phiW].PtF );
                tr->SetBranchAddress( Form("PtB_%s", brNamePostfix.Data() ), &br[etaW][phiW].PtB );
            }
    }

    // ##### QA pre-loop (for mult binning)

    TH1D *hist1D_QA_percentilesEstimator = new TH1D( "hist1D_QA_percentilesEstimator", "hist1D_QA_percentilesEstimator;percentile;entries", 3001, -0.5, 300.5);
    TH1D *hist1D_QA_multALL = new TH1D( "hist1D_QA_multALL", "hist1D_QA_multALL;mult;entries", 3001, -0.5, 3000.5);

    TH2D *hist2D_ESTIMATOR_VS_multTPC = new TH2D( "hist2D_ESTIMATOR_VS_multTPC", "hist2D_ESTIMATOR_VS_multTPC;estimator;mult in TPC", 4080, -2, 100, 301, -0.5, 3000.5);


    // ##### pre-loop over events for mult bins
    int nAccepted_PRE_LOOP = 0;
    int nRejectedByPileup_PRE_LOOP = 0;
    int treeId = 0;
    int sumEvPrevTrees = 0;//nEvInTrees[0];
    for (Long64_t i=0; i < nEvents; i++)
    {
        if ( i % 100000 == 0 )
            cout << "pre-loop: getting " << (int)i << endl;


        if ( i >= sumEvPrevTrees + nEvInTrees[treeId] ) // it's time to go to next tree...
        {
            sumEvPrevTrees += nEvInTrees[treeId]; // sum of events in prev trees
            treeId++;
            cout << "check treeId = " << treeId << " brV0M = " << brV0M << endl;
        }
        trees[treeId]->GetEntry( i - sumEvPrevTrees );

//        cout << "check i = " << i - (treeId>0 ? nEvInTrees[treeId-1] : 0) << endl;



//        if ( i < nEvInTrees[0] )
//        {
//            trees[0]->GetEntry( i );
//            if ( i < 5 )
//                cout << "check i = " << i << " brV0M = " << brV0M << endl;
//        }
//        else
//        {
//            trees[1]->GetEntry( i-nEvInTrees[0] );
//            if ( i < nEvInTrees[0]+5 )
//                cout << "check i = " << i << " brV0M = " << brV0M << endl;
//        }

        // ### event selection
        float cEstimator;
        bool isEventSelected = checkEventSelection(
                    cEstimator
                    , brIsPileupSPD
                    , brZEMvsZDC
                    , brV0M
                    , brCL1
                    , br_vertexZ
                    , nRejectedByPileup_PRE_LOOP
                    );

        if ( !isEventSelected )
            continue;

        hist1D_QA_percentilesEstimator->Fill(cEstimator);

        //calc mult in whole TPC using wins:
        int multTPC = 0;
        for ( int phiW = 0; phiW < nPhiWins; phiW++ )
        {
            if ( !FLAG_PERCOLATING_WINS )
                for ( int etaW = 0; etaW < nEtaBr; etaW++ )
                    multTPC += br[etaW][phiW].nF + br[etaW][phiW].nB;
            else
                multTPC += (br[0][phiW].nF + br[2][phiW].nF) + (br[0][phiW].nB + br[2][phiW].nB);
        }

        //        cout << "multTPC=" << multTPC << endl;
        hist1D_QA_multALL->Fill(multTPC);

        hist2D_ESTIMATOR_VS_multTPC->Fill( cEstimator, multTPC );

        nAccepted_PRE_LOOP++;
    } // end of pre-loop


    cout << ">>> nRejectedByPileup_PRE_LOOP = " << nRejectedByPileup_PRE_LOOP << endl;



    // ########## QA PLOTTING:

    //percentiles QA hist
    TCanvas *canv_estimatorPercentiles_QA_all = new TCanvas("canv_estimatorPercentiles_QA_all","canv_estimatorPercentiles_QA_all",0,0,700,600 );
    tuneCanvas(canv_estimatorPercentiles_QA_all);
    hist1D_QA_percentilesEstimator->DrawCopy();

    TCanvas *canv_hist1D_QA_multALL = new TCanvas("canv_hist1D_QA_multALL","canv_hist1D_QA_multALL",50,50,700,600 );
    tuneCanvas(canv_hist1D_QA_multALL);
    canv_hist1D_QA_multALL->SetLogy();
    hist1D_QA_multALL->DrawCopy();

    TCanvas *canv_ESTIMATOR_VS_multTPC = new TCanvas("canv_ESTIMATOR_VS_multTPC","canv_ESTIMATOR_VS_multTPC",10,10,800,800 );
    tuneCanvas(canv_ESTIMATOR_VS_multTPC);
    canv_ESTIMATOR_VS_multTPC->SetLogz();
    tuneHist2D(hist2D_ESTIMATOR_VS_multTPC);
    hist2D_ESTIMATOR_VS_multTPC->DrawCopy("colz");

    TProfile *prof_ESTIMATOR_VS_multTPC = hist2D_ESTIMATOR_VS_multTPC->ProfileX();
    prof_ESTIMATOR_VS_multTPC->SetLineColor(kBlue+1);
    //    prof_ESTIMATOR_VS_multTPC->DrawCopy("same");

    // ##### FUNCTIONAL CUT FOR OUTLIERS:
    TF1 *fBorderToCutOutliersLower = 0x0;
    TF1 *fBorderToCutOutliersUpper = 0x0;

    if ( isAnalysing502 == 0 )
    {
        fBorderToCutOutliersLower = new TF1("fBorderToCutOutliersLower","-120+1800*exp(-0.042*x)",0,90);
        fBorderToCutOutliersUpper = new TF1("fBorderToCutOutliersUpper","50+2450*exp(-0.038*x)",0,90);
    }
    else
    {
        //        fBorderToCutOutliersLower = new TF1("fBorderToCutOutliersLower","-150+2050*exp(-0.042*x)",0,90);
        fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","-100+2150*exp(-0.042*x)",0,90);
        fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersLower","120+2800*exp(-0.042*x)",0,90);
    }

    if (!isAMPT_kine)
    {
        fBorderToCutOutliersLower->SetLineColor(kRed+1);
        fBorderToCutOutliersLower->DrawCopy("same");

        fBorderToCutOutliersUpper->SetLineColor(kRed+2);
        fBorderToCutOutliersUpper->DrawCopy("same");
    }



//    return;

    // ########## Centrality bins:
    //    const int nCW = 2; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 10, 5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 9, 18 }; //n centrality bins

    //    const int nCW = 3; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 2.5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 10, 5, 2.5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 9, 17, 35 }; //n centrality bins


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

    //    const int nCW = 5; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 2.5, 1.0, 0.5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5,  2.5, 1.0, 1.0 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35, 36,  90, 90 }; //n centrality bins

    const int nCW = 1; //nCentrWidths
    const double cWidths[nCW] = { 10 }; //width of the centrality bins
    const double cStep[nCW] = { 5 }; //centrality bins step
    const int nCentrBins[nCW] = { 17 }; //n centrality bins

    //        const int nCW = 2; //nCentrWidths
    //        const double cWidths[nCW] = { 10, 5.001 }; //width of the centrality bins
    //        const double cStep[nCW] = { 5, 2.5 }; //centrality bins step
    //        const int nCentrBins[nCW] = { 17, 35 }; //n centrality bins





    // Split QA mult hist into quantiles: FOR MULT CLASSES
    double **multBounds = new double*[nCW];
    double **multBinCenters = new double*[nCW]; //by mean of 1D-histograms in bins!

    double **boundsMin = new double*[nCW];
    double **boundsMax = new double*[nCW];

    for ( int cW = 0; cW < nCW; cW++ )
    {
        int nCBinsForQuant = nCentrBins[cW]+1;
        cout << "###### Quantiles for estimator: n bins = " << nCBinsForQuant << endl;

        multBounds[cW] = new double[nCBinsForQuant]; // array to contain the quantiles
        getQuantiles(hist1D_QA_multALL, nCBinsForQuant, multBounds[cW]);

        multBinCenters[cW] = new double[nCBinsForQuant]; // array to contain means of 1D-histograms in bins
        drawCanvasWithClasses(hist1D_QA_multALL, Form("byMultTPC_%d_bins", nCBinsForQuant)
                              , nCBinsForQuant, multBounds[cW], multBinCenters[cW] );

        int nCBins = nCentrBins[cW];
        boundsMin[cW] = new double[nCBins];
        boundsMax[cW] = new double[nCBins];
        rearrangeBoundaries(nCBins, multBounds[cW], boundsMin[cW], boundsMax[cW] );

        for ( int bin = 0; bin < nCBins; bin++ )
            multBinCenters[cW][bin] = (multBinCenters[cW][bin]+multBinCenters[cW][bin+1])/2;

    }

    //    int nCentrBinsMult = 10;
    //    cout << "nCentrBins=" << nCentrBinsMult << endl;
    //    double *estBounds = new double[nCentrBinsMult]; // array to contain the quantiles
    //    getQuantiles(hist1D_QA_multALL, nCentrBinsMult, estBounds);
    //    drawCanvasWithClasses( hist1D_QA_multALL, "byMultTPC", nCentrBinsMult, estBounds );



    //        return;





    // ########## Select and initiate windows:
    //    const int nEtaWins = 3;
    //    const int howMany = 1; //how many windows to merge from branches
    //    const int nEtaWins = 2;
    //    const int howMany = 4; //how many windows to merge from branches
    const int nEtaWins = 1;
    const int howMany = 4; //how many windows to merge from branches

    const int maxNCentrBins = 100; //TMath::MaxElement(nCW, &nCentrBins);
    WinPair wins[nCW][maxNCentrBins][nEtaWins][nPhiWins];
    CentralityOccupancy cOccupancy[nCW][maxNCentrBins];

    bool useCentrPercOrMult = 0; // 0 - perc bins, 1-mult bins
    bool doBootstrap = 1;

    for ( int cW = 0; cW < nCW; cW++ )
        for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
        {

            float cBinMin = 0;
            float cBinMax = 0;

            if ( useCentrPercOrMult==0 ) //bins according to centrality percentiles
            {
                cBinMin = cStep[cW] * cBin;
                cBinMax = cWidths[cW] + cStep[cW] * cBin;
            }
            else //bins according to mult bins
            {
                //                cBinMin = ( cBin==0 ? 0 : multBounds[cW][cBin-1] );
                //                cBinMax = multBounds[cW][cBin];
                cBinMin = boundsMin[cW][cBin];
                cBinMax = boundsMax[cW][cBin];
            }

            cOccupancy[cW][cBin].cBinMin = cBinMin;
            cOccupancy[cW][cBin].cBinMax = cBinMax;

            for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                {
                    if ( !doBootstrap ) // CHECK!!!
                        wins[cW][cBin][etaW][phiW].init(cBinMin, cBinMax, etaW, phiW);
                    else
                        wins[cW][cBin][etaW][phiW].init(cBinMin, cBinMax, etaW, phiW, nAccepted_PRE_LOOP);
                }
        }





    //    return;



    // ##### main loop over events
    int nAccepted = 0;
    int nRejectedByPileup = 0;
    treeId = 0;
    sumEvPrevTrees = 0; //nEvInTrees[0];
    for (Long64_t i=0; i < nEvents; i++)
    {
        if ( i % 100000 == 0 )
            cout << "getting " << (int)i << endl;
        //                cout <<"getting " << (int)i << "\r"; cout.flush();

        if ( i >= sumEvPrevTrees + nEvInTrees[treeId] ) // it's time to go to next tree...
        {
            sumEvPrevTrees += nEvInTrees[treeId]; // sum of events in prev trees
            treeId++;
            cout << "check treeId = " << treeId << " brV0M = " << brV0M << endl;
        }
        trees[treeId]->GetEntry( i - sumEvPrevTrees );

//        cout << "check i = " << i - (treeId>0 ? nEvInTrees[treeId-1] : 0) << endl;


//        if ( i >= nEvInTrees[treeId] )
//        {
//            treeId++;
////            if ( i < 5 )
//                cout << "check treeId = " << treeId << " brV0M = " << brV0M << endl;
//        }
//        trees[treeId]->GetEntry( i - (treeId>0 ? nEvInTrees[treeId-1] : 0) );

//        cout << "check i = " << i - (treeId>0 ? nEvInTrees[treeId-1] : 0) << endl;

//        if ( i < nEvInTrees[0] )
//        {
//            trees[0]->GetEntry( i );
//            if ( i < 5 )
//                cout << "check i = " << i << " brV0M = " << brV0M << endl;
//        }
//        else
//        {
//            trees[1]->GetEntry( i-nEvInTrees[0] );
//            if ( i < nEvInTrees[0]+5 )
//                cout << "check i = " << i << " brV0M = " << brV0M << endl;
//        }


        // ### event selection
        float cEstimator;
        bool isEventSelected = checkEventSelection(
                    cEstimator
                    , brIsPileupSPD
                    , brZEMvsZDC
                    , brV0M
                    , brCL1
                    , br_vertexZ
                    , nRejectedByPileup
                    );

        if ( !isEventSelected )
            continue;

        //calc mult in whole TPC using wins:
        int multTPC = 0;
        for ( int phiW = 0; phiW < nPhiWins; phiW++ )
        {
            if ( !FLAG_PERCOLATING_WINS )
                for ( int etaW = 0; etaW < nEtaBr; etaW++ )
                    multTPC += br[etaW][phiW].nF + br[etaW][phiW].nB;
            else
                multTPC += (br[0][phiW].nF + br[2][phiW].nF) + (br[0][phiW].nB + br[2][phiW].nB);
        }



        // !!!! test cut by line on Perc_vs_mult plot:
        if ( !isAMPT_kine )
            if ( multTPC < fBorderToCutOutliersLower->Eval(cEstimator)
                 || multTPC > fBorderToCutOutliersUpper->Eval(cEstimator) )
                continue;


        //assign "centrality" for this event: either cPerc or multTPC
        float centrValue = ( useCentrPercOrMult==0 ? cEstimator : multTPC );

        //fill wins
        for ( int cW = 0; cW < nCW; cW++ )
            for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
            {
                cOccupancy[cW][cBin].fill(brV0M, brZEMvsZDC);
                for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                {
                    // 27.02.16: fill wins "by hand" using info from branches:
                    UShort_t _nF = 0;
                    UShort_t _nB = 0;
                    Float_t _nF_PtF = 0;
                    Float_t _nB_PtB = 0;
                    Float_t _PtF = -1; // important to put -1!
                    Float_t _PtB = -1; // important to put -1!
                    for ( int _eW = howMany*etaW; _eW < howMany*(etaW+1); _eW++ )
                    {
                        //                        cout << "etaW=" << etaW << ", _eW=" << _eW << endl;
                        //                        int a;
                        //                        cin >> a;
                        _nF  += br[_eW][0].nF;
                        _nB  += br[_eW][0].nB;
                        _nF_PtF += br[_eW][0].nF * br[_eW][0].PtF;
                        _nB_PtB += br[_eW][0].nB * br[_eW][0].PtB;
                    }
                    if ( _nF > 0 )
                        _PtF = _nF_PtF / _nF;
                    if ( _nB > 0 )
                        _PtB = _nB_PtB / _nB;
                    wins[cW][cBin][etaW][0].fill( centrValue, _nF, _nB, _PtF, _PtB );
                }

                //... IA: po-normalnomu bylo tak:
                //                for ( int etaW = 0; etaW < nEtaBr; etaW++ )
                //                    for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                //                        wins[cW][cBin][etaW][phiW].fill( centrValue, br[etaW][phiW].nF, br[etaW][phiW].nB, br[etaW][phiW].PtF, br[etaW][phiW].PtB );
            }
        nAccepted++;
    } // end of events
    cout << "nAccepted = " << nAccepted << endl;
    cout << "nAccepted/nAll = " << (float)nAccepted/nEvents << endl;

    cout << ">>> main loop: nRejectedByPileup = " << nRejectedByPileup << endl;





    // ########## PREPARE OUTPUT ROOT-FILE:
    TFile *fileOutput = new TFile( "output_histos_graphs.root", "RECREATE" );

    // ########## BOOTSTRAPING
    TGraphErrors *gr_BS_PtPt[nCW][nEtaWins][nPhiWins];// = new TGraphErrors;
    if ( doBootstrap )
    {
        cout << "Start bootstrapping..." << endl;

        TCanvas *canv_bootStrapPhiWins = new TCanvas("canv_bootStrapPhiWins","canv_bootStrapPhiWins",350,150,900,700 );

        for ( int cW = 0; cW < nCW; cW++ )
        {
            //                for ( int etaW = 0; etaW < nEtaWins; etaW++ )
            for ( int phiW = 0; phiW < nPhiWins; phiW++ )
            {
                int etaW = 0; // TMP, FIXED FOR TESTS!
                gr_BS_PtPt[cW][etaW][phiW] = new TGraphErrors;
                gr_BS_PtPt[cW][etaW][phiW]->SetName( Form("gr_BS_PtPt_cW%d_etaW%d_phiW%d", cW, etaW, phiW) );
                for ( int cBin = /*nCentrBins[cW]-3*/ 0; cBin < nCentrBins[cW]; cBin++ )
                {
                    cout << "  processing cBin " << cBin << "..." << endl;

                    WinPair *w = &wins[cW][cBin][etaW][phiW];

                    // DO BOOTSTRAP
                    w->performBootstrapping(1);

                    double BS_bCorr_mean = w->hist1D_bCorr_BS_PtPt->GetMean();
                    double BS_bCorr_sigma = w->hist1D_bCorr_BS_PtPt->GetRMS();

                    float centr = -1;
                    if ( useCentrPercOrMult == 0 )
                        centr = w->cBinMin + (w->cBinMax - w->cBinMin)/2;
                    else
                        centr = multBinCenters[cW][cBin];

                    gr_BS_PtPt[cW][etaW][phiW]->SetPoint( cBin, centr, BS_bCorr_mean );
                    gr_BS_PtPt[cW][etaW][phiW]->SetPointError( cBin, 0, BS_bCorr_sigma );


                    w->hist1D_bCorr_BS_PtPt->SetLineColor(kOrange - 5 + cW);

                    if (cBin==0)
                        w->hist1D_bCorr_BS_PtPt->DrawCopy();
                    else
                        w->hist1D_bCorr_BS_PtPt->DrawCopy("same");

                    w->hist1D_bCorr_BS_PtPt->Write();
                } // end of cBin
                gr_BS_PtPt[cW][etaW][phiW]->SetLineColor(kMagenta);
                gr_BS_PtPt[cW][etaW][phiW]->SetMarkerColor(kMagenta);
                gr_BS_PtPt[cW][etaW][phiW]->SetMarkerStyle(24);
                gr_BS_PtPt[cW][etaW][phiW]->Write();
            } // end of phiW
        } // end of cW

        TCanvas *canv_GrCoeff = new TCanvas("canv_GrCoeff","canv_GrCoeff",20,150,900,700 );
        tuneCanvas(canv_GrCoeff);
        //        grC2->Draw("APL");

        //        grFromFit2D->SetLineColor(kRed);
        //        grFromFit2D->DrawClone("PL");

        tuneGraphAxisLabels(gr_BS_PtPt[0][0][0]);
        gr_BS_PtPt[0][0][0]->DrawClone("APL");

    } // end of bootstrap


    //    return;



    // ########## SAVE HISTOS TO ROOT-FILE:
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int cBin = /*nCentrBins[cW]-3*/ 0; cBin < nCentrBins[cW]; cBin++ )
            for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                {
                    wins[cW][cBin][etaW][phiW].hist2D_NN->Write();
                    wins[cW][cBin][etaW][phiW].hist2D_PtPt->Write();
                    wins[cW][cBin][etaW][phiW].hist2D_PtN->Write();

                    wins[cW][cBin][etaW][phiW].hist1D_multDistrF->Write();
                    wins[cW][cBin][etaW][phiW].hist1D_multDistrB->Write();

                    wins[cW][cBin][etaW][phiW].hist1D_QA_PtF->Write();
                    wins[cW][cBin][etaW][phiW].hist1D_QA_PtB->Write();

                    //                    wins[cW][cBin][etaW][phiW].hist2D_NN->ProfileX()->Write();
                    //                    wins[cW][cBin][etaW][phiW].hist2D_PtPt->ProfileX()->Write();
                    //                    wins[cW][cBin][etaW][phiW].hist2D_PtN->ProfileX()->Write();
                }

    hist1D_QA_percentilesEstimator->Write();
    hist1D_QA_multALL->Write();
    hist2D_ESTIMATOR_VS_multTPC->Write();
    hist2D_ESTIMATOR_VS_multTPC->ProfileX()->Write();

    //    fileOutput->WriteObject(canv_ESTIMATOR_VS_multTPC);
    canv_ESTIMATOR_VS_multTPC->Write();

    // ########## MAIN PLOTTING FOR CORRS:

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

            grNN[cW][etaW]->SetName(   Form( "grNN_c%d_eta%d", cW, etaW ) );
            grPtPt[cW][etaW]->SetName( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
            grPtN[cW][etaW]->SetName(  Form( "grPtN_c%d_eta%d", cW, etaW ) );

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
                    float centr = -1;
                    if ( useCentrPercOrMult == 0 )
                        centr = w->cBinMin + (w->cBinMax - w->cBinMin)/2;
                    else
                        centr = multBinCenters[cW][cBin];

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

    //    for ( int cW = 0; cW < nCW; cW++ )
    //        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
    //        {
    //            grNN[cW][etaW]->Write();
    //            grPtPt[cW][etaW]->Write();
    //            grPtN[cW][etaW]->Write();
    //        }

    //    fileOutput->Close();





    // ########### draw graphs ###########
    int colors[] = { kBlack, kBlack, kBlue, kRed, kMagenta };
    int markers[] = { 20, 24, 5, 2, 21 };

    const int etaId = 0;


    // ########## NN
    TCanvas *canv_grNN = new TCanvas("canv_grNN","canv_grNN",20,50,700,600 );
    tuneCanvas(canv_grNN);
    grNN[0][etaId]->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels(grNN[0][etaId]);


    //centr
    for ( int cW = 0; cW < nCW; cW++ )
        drawGraph(grNN[cW][etaId], markers[cW], colors[cW], cW == 0 ? "AP" : "P");

    grNN[0][etaId]->SetMinimum( 0 );

    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grNN[cW][etaId], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();

    TLatex *tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);

    tex = new TLatex(0.4,0.89, "NN");
    drawTex(tex, 0.045);


    TString strPostfix;

    if ( flag_V0M_ZDC_CL1 == 0 )
        strPostfix = Form("V0M.eps");
    else if ( flag_V0M_ZDC_CL1 == 1 )
        strPostfix = Form("ZDCZEM.eps");
    else if ( flag_V0M_ZDC_CL1 == 2 )
        strPostfix = Form("CL1.eps");

    canv_grNN->SaveAs( Form("NN_%s", strPostfix.Data() ) );

    // ########## PtPt
    TCanvas *canv_grPtPt = new TCanvas("canv_grPtPt","canv_grPtPt",250,50,700,600 );
    tuneCanvas(canv_grPtPt);
    grPtPt[0][etaId]->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels(grPtPt[0][etaId]);


    //centr
    for ( int cW = 0; cW < nCW; cW++ )
        drawGraph(grPtPt[cW][etaId], markers[cW], colors[cW], cW == 0 ? "AP" : "P");

    if ( doBootstrap )
    {
        gr_BS_PtPt[0][0][0]->Draw("P");
        //        gr_PtPt_c10_BS->Write();
    }

    grPtPt[0][etaId]->SetMinimum( 0 );

    leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grPtPt[cW][etaId], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();

    tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);

    tex = new TLatex(0.4,0.89, "PtPt");
    drawTex(tex, 0.045);


    canv_grNN->SaveAs( Form("PtPt_%s", strPostfix.Data() ) );



    // ########## PtN
    TCanvas *canv_grPtN = new TCanvas("canv_grPtN","canv_grPtN",450,50,700,600 );
    tuneCanvas(canv_grPtN);
    grPtN[0][etaId]->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels(grPtN[0][etaId]);


    //centr 10
    for ( int cW = 0; cW < nCW; cW++ )
        drawGraph(grPtN[cW][etaId], markers[cW], colors[cW], cW == 0 ? "AP" : "P");

    grPtN[0][etaId]->SetMinimum( 0 );
    leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grPtN[cW][etaId], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();

    tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);

    tex = new TLatex(0.4,0.89, "PtN");
    drawTex(tex, 0.045);


    canv_grPtN->SaveAs( Form("PtN_%s", strPostfix.Data() ) );



    TCanvas *canv_grPtN_2D = new TCanvas("canv_grPtN_2D","canv_grPtN_2D",450,50,700,600 );
    tuneCanvas(canv_grPtN_2D);


    //    wins[0][0][0][0].hist2D_PtN->DrawCopy();
    wins[0][0][0][0].hist2D_PtN->ProfileX()->DrawCopy();



    // ########## CENTR ESTIMATOR EVENT RATIO:
    TCanvas *canv_grCentrRatio = new TCanvas("canv_grCentrRatio","canv_grCentrRatio",450,250,700,600 );
    tuneCanvas(canv_grCentrRatio);
    grFractEstByV0M[0]->SetTitle(";centrality percentile;ratio");
    tuneGraphAxisLabels(grFractEstByV0M[0]);

    //centr
    for ( int cW = 0; cW < nCW; cW++ )
        drawGraph(grFractEstByV0M[cW], markers[cW], colors[cW], cW == 0 ? "AP" : "P");

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
    gPad->SetLogy();
    for ( int cBin = 0; cBin < nCentrBins[0]; cBin++ )
    {
        TH1D *h = wins[0][cBin][0][0].hist1D_multDistrF;
        h->SetLineColor(kOrange-9+cBin);

        if ( cBin == 0 )
        {
            h->SetLineColor(kRed);
            h->GetYaxis()->SetRangeUser(1,100000);
        }

        if ( cBin == 0 )
            h->DrawCopy();
        else
            h->DrawCopy("same");
    }


    // Check entries in centrality bins:
    for ( int cW = 0; cW < nCW; cW++ )
    {
        cout << " ###### cW = " << cW << endl;
        for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
        {
            TH1D *h = wins[cW][cBin][0][0].hist1D_multDistrF;
            cout << "entries in cBin " << cBin << " = " << h->GetEntries() << endl;
            //            h->SetLineColor(kOrange-9+cBin);
            //            if ( cBin == 0 )
            //                h->DrawCopy();
            //            else
            //                h->DrawCopy("same");
        }
    }



    // ###### write graphs to output file
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
        {
            grNN[cW][etaW]->Write();
            grPtPt[cW][etaW]->Write();
            grPtN[cW][etaW]->Write();
        }

    fileOutput->Close();

    for (Int_t i=0; i<nFiles; i++)
        myFiles[i]->Close();

    return;
}



