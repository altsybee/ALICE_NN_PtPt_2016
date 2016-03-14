void draw_graphs_both_energies_ptpt()
{
    gROOT->ProcessLine(".L utils.C");

    TFile *f[8];

    f[0] = new TFile( "data_graphs/energy_502_fieldPP_graphPtPt_eta0_c1.root" );
    f[1] = new TFile( "data_graphs/energy_502_fieldPP_graphPtPt_eta0_c0.root" );

    f[2] = new TFile( "data_graphs/energy_276_graphPtPt_eta0_c1.root" );
    f[3] = new TFile( "data_graphs/energy_276_graphPtPt_eta0_c0.root" );

    f[4] = new TFile( "data_graphs/energy_502_fieldMM_graphPtPt_eta0_c1.root" );
    f[5] = new TFile( "data_graphs/energy_502_fieldMM_graphPtPt_eta0_c0.root" );

    f[6] = new TFile( "data_graphs/energy_502_MERGED_fields_graphPtPt_eta0_c1.root" );
    f[7] = new TFile( "data_graphs/energy_502_MERGED_fields_graphPtPt_eta0_c0.root" );



    TGraphErrors *gr[8];
    gr[0] = (TGraphErrors*)f[0]->Get("graphPtPt_eta0_c1");
    gr[1] = (TGraphErrors*)f[1]->Get("graphPtPt_eta0_c0");

    gr[2] = (TGraphErrors*)f[2]->Get("graphPtPt_eta0_c1");
    gr[3] = (TGraphErrors*)f[3]->Get("graphPtPt_eta0_c0");

    gr[4] = (TGraphErrors*)f[4]->Get("graphPtPt_eta0_c1");
    gr[5] = (TGraphErrors*)f[5]->Get("graphPtPt_eta0_c0");

    gr[6] = (TGraphErrors*)f[6]->Get("graphPtPt_eta0_c1");
    gr[7] = (TGraphErrors*)f[7]->Get("graphPtPt_eta0_c0");

//    drawGraph( gr[0], 24, kBlue, "PL" );
    drawGraph( gr[1], 20, kBlue, "APL" );

//    drawGraph( gr[2], 24, kBlack, "PL" );
    drawGraph( gr[3], 20, kBlack, "PL" );

//    drawGraph( gr[4], 24, kGreen, "PL" );
    drawGraph( gr[5], 20, kGreen, "PL" );

//    drawGraph( gr[6], 24, kRed, "PL" );
    drawGraph( gr[7], 20, kRed, "PL" );

//    gROOT->ProcessLine( ".q");
}
