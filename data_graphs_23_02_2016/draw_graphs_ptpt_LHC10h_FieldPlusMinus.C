void draw_graphs_ptpt_LHC10h_FieldPlusMinus()
{
    gROOT->ProcessLine(".L ../utils.C");

    TFile *f[8];

    f[0] = new TFile( "new_eW_binning_27_02_2016/output_histos_graphs_LHC10h_eW04_MFplus.root" );
    f[1] = new TFile( "new_eW_binning_27_02_2016/output_histos_graphs_LHC10h_eW04_MFminus.root" );



    TGraphErrors *gr[10][10][10];
    for ( int cW = 0; cW < 2; cW++ )
            for ( int etaW = 0; etaW < 2; etaW++ )
            {
                gr[0][cW][etaW] = (TGraphErrors*)f[0]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
                gr[1][cW][etaW] = (TGraphErrors*)f[1]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
            }


    int etaId = 0;

    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",150,50,700,600 );
    tuneCanvas(canv_graphs);

    tuneGraphAxisLabels(gr[0][0][etaId]);

    gr[0][0][etaId]->SetTitle( ";centrality percentile;b_{corr}" );

    //LHC10h MF plus
    drawGraph( gr[0][0][etaId], 20, kBlack, "APL" );
    drawGraph( gr[0][1][etaId], 24, kBlack, "PL" );

    //LHC10h MF minus
    drawGraph( gr[1][0][etaId], 20, kBlue, "PL" );
    drawGraph( gr[1][1][etaId], 24, kBlue, "PL" );

    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    leg->AddEntry( gr[0][0][etaId], "LHC10h MF++, class 10", "p");
    leg->AddEntry( gr[0][1][etaId], "LHC10h MF++, class 5", "p");
    leg->AddEntry( gr[1][0][etaId], "LHC10h MF--, class 10", "p");
    leg->AddEntry( gr[1][1][etaId], "LHC10h MF--, class 5", "p");

    leg->Draw();


//    gROOT->ProcessLine( ".q");
}
