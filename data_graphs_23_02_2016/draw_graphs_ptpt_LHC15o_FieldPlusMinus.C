void draw_graphs_ptpt_LHC15o_FieldPlusMinus()
{
    gROOT->ProcessLine(".L ../utils.C");

    TFile *f[8];

    f[0] = new TFile( "output_classesByV0M_LHC10h_c10_5_25_1_05.root" );

    f[1] = new TFile( "output_classesByV0M_LHC15o_fieldMM_c10_5_CUT_OUTLIERS.root" );
    f[2] = new TFile( "output_classesByV0M_LHC15o_fieldPP_c10_5_CUT_OUTLIERS.root" );

    f[3] = new TFile( "output_classesByV0M_LHC11h_FemtoMinus_c10_5_CUT_OUTLIERS.root" );
    f[4] = new TFile( "output_classesByV0M_LHC11h_FemtoPlus_c10_5_CUT_OUTLIERS.root" );

    f[5] = new TFile( "output_classesByV0M_LHC15o_fieldCOMBINED_c10_5_CUT_OUTLIERS.root" );
    f[6] = new TFile( "output_classesByV0M_LHC11h_FemtoCOMBINED_c10_5_CUT_OUTLIERS.root" );



    TGraphErrors *gr[10][10][10];
    for ( int cW = 0; cW < 2; cW++ )
            for ( int etaW = 0; etaW < 3; etaW++ )
            {
                gr[0][cW][etaW] = (TGraphErrors*)f[0]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );

                gr[1][cW][etaW] = (TGraphErrors*)f[1]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
                gr[2][cW][etaW] = (TGraphErrors*)f[2]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );

                gr[3][cW][etaW] = (TGraphErrors*)f[3]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
                gr[4][cW][etaW] = (TGraphErrors*)f[4]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );

                gr[5][cW][etaW] = (TGraphErrors*)f[5]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
                gr[6][cW][etaW] = (TGraphErrors*)f[6]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
            }


    int etaId = 0;

    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",150,50,700,600 );
    tuneCanvas(canv_graphs);

    tuneGraphAxisLabels(gr[0][0][etaId]);

    gr[0][0][etaId]->SetTitle( ";centrality percentile;b_{corr}" );

    //LHC10h
    drawGraph( gr[0][0][etaId], 20, kBlack, "APL" );
//    drawGraph( gr[0][1][etaId], 24, kBlack, "PL" );


    //LHC15o
    //MM
//    drawGraph( gr[1][0][etaId], 20, kRed, "PL" );
//    drawGraph( gr[1][1][etaId], 24, kRed, "PL" );

    //PP
//    drawGraph( gr[2][0][etaId], 21, kRed+2, "PL" );
//    drawGraph( gr[2][1][etaId], 24, kRed+2, "PL" );

    // COMBINED MM+PP:
    drawGraph( gr[5][0][etaId], 20, kRed, "PL" );
//    drawGraph( gr[5][1][etaId], 24, kRed, "PL" );

    //LHC11h
    //FemtoMinus
//    drawGraph( gr[3][0][etaId], 20, kBlue, "PL" );
//    drawGraph( gr[3][1][etaId], 24, kBlue, "PL" );

    //FemtoPlus
//    drawGraph( gr[4][0][etaId], 21, kBlue+2, "PL" );
//    drawGraph( gr[4][1][etaId], 24, kBlue+2, "PL" );

    // COMBINED FemtoPlusMinus:
    drawGraph( gr[6][0][etaId], 20, kBlue, "PL" );
//    drawGraph( gr[6][1][etaId], 24, kBlue, "PL" );

    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    leg->AddEntry( gr[0][0][etaId], "LHC10h", "p");
    leg->AddEntry( gr[6][0][etaId], "LHC11h MF++,--", "p");
    leg->AddEntry( gr[5][0][etaId], "LHC15o MF++,--", "p");
//    leg->AddEntry( gr[3][0][etaId], "LHC11h MF--", "p");
//    leg->AddEntry( gr[4][0][etaId], "LHC11h MF++", "p");
//    leg->AddEntry( gr[1][0][etaId], "LHC15o MF--", "p");
//    leg->AddEntry( gr[2][0][etaId], "LHC15o MF++", "p");

    leg->Draw();


//    gROOT->ProcessLine( ".q");
}
