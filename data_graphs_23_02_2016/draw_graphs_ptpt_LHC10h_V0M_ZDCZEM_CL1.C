void draw_graphs_ptpt_LHC10h_V0M_ZDCZEM_CL1()
{
    gROOT->ProcessLine(".L ../utils.C");

    TFile *f[8];

    f[0] = new TFile( "output_classesByV0M_LHC10h_c10_5_25_1_05_CUT_OUTLIERS.root" );

    f[1] = new TFile( "output_classesByZDCZEM_LHC10h_c10_5_25_1_05.root" );
    f[2] = new TFile( "output_classesByCL1_LHC10h_c10_5_25_1_05.root" );


    TGraphErrors *gr[10][10][10];
    for ( int cW = 0; cW < 2; cW++ )
        for ( int etaW = 0; etaW < 3; etaW++ )
        {
            gr[0][cW][etaW] = (TGraphErrors*)f[0]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );

            gr[1][cW][etaW] = (TGraphErrors*)f[1]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );
            gr[2][cW][etaW] = (TGraphErrors*)f[2]->Get( Form( "grPtPt_c%d_eta%d", cW, etaW ) );

        }


    int etaId = 0;

    TCanvas *canv_graphs = new TCanvas("canv_graphs","canv_graphs",150,250,700,600 );
    tuneCanvas(canv_graphs);

    tuneGraphAxisLabels(gr[0][0][etaId]);

    //LHC10h
    drawGraph( gr[0][0][etaId], 20, kBlack, "APL" );
//        drawGraph( gr[0][1][etaId], 24, kBlack, "PL" );


    drawGraph( gr[1][0][etaId], 20, kBlue, "PL" );
//        drawGraph( gr[1][1][etaId], 24, kBlue, "PL" );

    drawGraph( gr[2][0][etaId], 20, kOrange+9, "PL" );
//        drawGraph( gr[2][1][etaId], 24, kOrange+9, "PL" );


    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    leg->AddEntry(gr[0][0][etaId], "V0M", "p");
    leg->AddEntry(gr[1][0][etaId], "ZDCvsZEM", "p");
    leg->AddEntry(gr[2][0][etaId], "CL1", "p");

    leg->Draw();

    //    graphCorr[0][etaId]->GetYaxis()->SetRangeUser( 0, 0.92 );

    //    TLatex *tex = new TLatex(0.4,0.89, "#eta_{gap}=0.8, #delta#eta=0.4");
    TLatex *tex = new TLatex(0.5,0.3, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);



    //    gROOT->ProcessLine( ".q");
}
