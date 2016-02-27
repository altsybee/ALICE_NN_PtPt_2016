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

#include "SupplementaryClasses.h"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;

#include "utils.C"


void get_nn_ptpt_FROM_FB_TREE() // TString inputFileName = "MergedOutput.root")
{
    // FULL STATISTICS
    // LHC10h
    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFplus_all/MergedOutput.root" );




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

    if (!myFile)
    {
        cout << "No input file!" << endl;
        return;
    }
    myFile->ls();

    // !!! important for branch V0M below!
    bool isAnalysing502 = 0;

    int listId = 0;//1;
    // !!! listId can be different for 2.76 and 5.02!
    //    if ( isAnalysing502 )
    //        listId = 1;

    TList *listKeys = myFile->GetListOfKeys();
    TString listName = Form("%s",listKeys->At(listId)->GetName() );
    while ( !listName.Contains("PWGCFLRC") )
    {
        listId++;
        listName = Form("%s",listKeys->At(listId)->GetName() );
    }
    //    cout << ".... " << listKeys->At(listId)->GetName() << endl;
    //return;
    cout << "going into list: " << listKeys->At(listId)->GetName() << endl;

    //return;
    myFile->cd( listKeys->At(listId)->GetName() );
    TList *listKeys2 = gDirectory->GetListOfKeys();
    TList *myTask = (TList*) gDirectory->Get( listKeys2->At(0)->GetName() );
    TTree *t1 = (TTree*)  myTask->FindObject( "treeFB" );

//    int nEvents = t1->GetEntries();
        int nEvents = 20000;
    cout <<"nEvents = " << nEvents << endl;


    Float_t brV0M;
    if ( !isAnalysing502 )
        t1->SetBranchAddress( "centr_V0M", &brV0M );
    else
        t1->SetBranchAddress( "centrV0M_NEW_MULT_SEL", &brV0M );

    Float_t brZEMvsZDC;
    t1->SetBranchAddress( "centr_ZEMvsZDC", &brZEMvsZDC );

    Float_t brCL1;
    t1->SetBranchAddress( "centr_CL1", &brCL1 );


    Float_t br_vertexZ = 0;
    t1->SetBranchAddress( "vertexZ", &br_vertexZ );

    const float kVertexZcut = 8; // for LHC11h
    //    if ( isAnalysing502 )
    //        kVertexZcut = 7;
    bool USE_VERTEXZ_CUT = 0;
    int flag_V0M_ZDC_CL1 = 0;


    // ########## Number of eta-phi wins
    //    const int nEtaWins = 3;
    const int nEtaWins = 8;
    const int nPhiWins = 1;

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

    // ##### QA pre-loop (for mult binning)

    TH1D *hist1D_QA_percentilesEstimator = new TH1D( "hist1D_QA_percentilesEstimator", "hist1D_QA_percentilesEstimator;percentile;entries", 3001, -0.5, 300.5);
    TH1D *hist1D_QA_multALL = new TH1D( "hist1D_QA_multALL", "hist1D_QA_multALL;mult;entries", 3001, -0.5, 3000.5);

    TH2D *hist2D_ESTIMATOR_VS_multTPC = new TH2D( "hist2D_ESTIMATOR_VS_multTPC", "hist2D_ESTIMATOR_VS_multTPC;estimator;mult in TPC", 4080, -2, 100, 301, -0.5, 3000.5);


    // ##### pre-loop over events for mult bins
    int nAccepted_PRE_LOOP = 0;
    for (Long64_t i=0; i < nEvents; i++)
    {
        if ( i % 100000 == 0 )
            cout << "pre-loop: getting " << (int)i << endl;

        t1->GetEntry( i );

        float cEstimator = -1;
        if ( flag_V0M_ZDC_CL1 == 0 )
        {
            if ( brV0M > 90 ) //V0M cut
                continue;
            cEstimator = brV0M;
        }
        else if ( flag_V0M_ZDC_CL1 == 1 )
        {
            if ( brZEMvsZDC > 50 ) //ZDCvsZEM cut
                continue;
            cEstimator = brZEMvsZDC;
        }
        else if ( flag_V0M_ZDC_CL1 == 2 )
        {
            if ( brCL1 > 90 ) //CL1 cut
                continue;
            cEstimator = brCL1;
        }

        // vertex Z cut!!!
        if ( USE_VERTEXZ_CUT && fabs(br_vertexZ) > kVertexZcut )
            continue;

        hist1D_QA_percentilesEstimator->Fill(cEstimator);


        //calc mult in whole TPC using wins:
        int multTPC = 0;
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
            for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                multTPC += br[etaW][phiW].nF + br[etaW][phiW].nB;  //(br[0][phiW].nF + br[2][phiW].nF) + (br[0][phiW].nB + br[2][phiW].nB);

        //        cout << "multTPC=" << multTPC << endl;
        hist1D_QA_multALL->Fill(multTPC);

        hist2D_ESTIMATOR_VS_multTPC->Fill( cEstimator, multTPC );

        nAccepted_PRE_LOOP++;
    } // end of pre-loop





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
    TF1 *fBorderToCutOutliers = 0x0;

    if ( isAnalysing502 == 0 )
        fBorderToCutOutliers = new TF1("fBorderToCutOutliers","-100+1600*exp(-0.042*x)",0,90);
    else
        //        fBorderToCutOutliers = new TF1("fBorderToCutOutliers","-150+2050*exp(-0.042*x)",0,90);
        fBorderToCutOutliers = new TF1("fBorderToCutOutliers","-100+2150*exp(-0.042*x)",0,90);

    fBorderToCutOutliers->SetLineColor(kRed+1);
    fBorderToCutOutliers->DrawCopy("same");

    //return;

    // ########## Centrality bins:
    //    const int nCW = 2; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 10, 5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 9, 18 }; //n centrality bins

    //    const int nCW = 3; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 2.5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 10, 5, 2.5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 9, 17, 35 }; //n centrality bins


    const int nCW = 2; //nCentrWidths
    const double cWidths[nCW] = { 10, 5 }; //width of the centrality bins
    const double cStep[nCW] = { 5, 2.5 }; //centrality bins step
    const int nCentrBins[nCW] = { 17, 35 }; //n centrality bins

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

    //    const int nCW = 1; //nCentrWidths
    //    const double cWidths[nCW] = { 10 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17 }; //n centrality bins

    //    const int nCW = 2; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5.001 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35 }; //n centrality bins





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
//    const nEtaWins
    const int maxNCentrBins = 100; //TMath::MaxElement(nCW, &nCentrBins);
    WinPair wins[nCW][maxNCentrBins][nEtaWins][nPhiWins];
    CentralityOccupancy cOccupancy[nCW][maxNCentrBins];

    bool useCentrPercOrMult = 0; // 0 - perc bins, 1-mult bins
    bool doBootstrap = 0;

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
    for (Long64_t i=0; i < nEvents; i++)
    {
        if ( i % 100000 == 0 )
            cout << "getting " << (int)i << endl;
        //                cout <<"getting " << (int)i << "\r"; cout.flush();

        t1->GetEntry( i );
        //            t1->GetEntry( TMath::Nint( gRandom->Uniform(-0.5,nEvents-0.5) ) );

        float cEstimator = -1;
        if ( flag_V0M_ZDC_CL1 == 0 )
        {
            if ( brV0M > 90 ) //V0M cut
                continue;
            cEstimator = brV0M;
        }
        else if ( flag_V0M_ZDC_CL1 == 1 )
        {
            if ( brZEMvsZDC > 50 ) //ZDCvsZEM cut
                continue;
            cEstimator = brZEMvsZDC;
        }
        else if ( flag_V0M_ZDC_CL1 == 2 )
        {
            if ( brCL1 > 90 ) //CL1 cut
                continue;
            cEstimator = brCL1;
        }


        // vertex Z cut!!!
        if ( USE_VERTEXZ_CUT && fabs(br_vertexZ) > kVertexZcut )
            continue;



        //calc mult in whole TPC using wins:
        int multTPC = 0;
        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
            for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                multTPC += br[etaW][phiW].nF + br[etaW][phiW].nB;  //(br[0][phiW].nF + br[2][phiW].nF) + (br[0][phiW].nB + br[2][phiW].nB);


        // !!!! test cut by line on Perc_vs_mult plot:
        if ( multTPC < fBorderToCutOutliers->Eval(cEstimator) )
            continue;


        //assign "centrality" for this event: either cPerc or multTPC
        float centrValue = ( useCentrPercOrMult==0 ? cEstimator : multTPC );

        //will wins
        for ( int cW = 0; cW < nCW; cW++ )
            for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
            {
                cOccupancy[cW][cBin].fill(brV0M, brZEMvsZDC);
                for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                    for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                    {
                        //br[etaW][phiW].nF, br[etaW][phiW].nB, br[etaW][phiW].PtF, br[etaW][phiW].PtB
                        wins[cW][cBin][etaW][phiW].fill( centrValue, br[etaW][phiW].nF, br[etaW][phiW].nB, br[etaW][phiW].PtF, br[etaW][phiW].PtB );
                    }
            }
        nAccepted++;
    } // end of events
    cout << "nAccepted = " << nAccepted << endl;
    cout << "nAccepted/nAll = " << (float)nAccepted/nEvents << endl;


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
                gr_BS_PtPt[cW][etaW][phiW]->Write();
            } // end of phiW
        } // end of cW

        TCanvas *canv_GrCoeff = new TCanvas("canv_GrCoeff","canv_GrCoeff",20,150,900,700 );
        tuneCanvas(canv_GrCoeff);
        //        grC2->Draw("APL");

        //        grFromFit2D->SetLineColor(kRed);
        //        grFromFit2D->DrawClone("PL");

        tuneGraphAxisLabels(gr_BS_PtPt[0][0][0]);
        gr_BS_PtPt[0][0][0]->SetLineColor(kMagenta);
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


    return;
}



