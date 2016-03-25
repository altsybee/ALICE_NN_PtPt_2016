void tuneGraphAxisLabels( TGraphErrors *gr )
{
    gr->GetYaxis()->SetTitleOffset( 1.45 );
    gr->GetYaxis()->SetTitleSize( 0.048 );
    gr->GetYaxis()->SetLabelSize( 0.042 );

    gr->GetXaxis()->SetTitleOffset( 0.95 );
    gr->GetXaxis()->SetTitleSize( 0.048 );
    gr->GetXaxis()->SetLabelSize( 0.042 );

}

void tuneHist1D( TH1D *hist )
{
    hist->GetYaxis()->SetTitleOffset( 1.45 );
    hist->GetYaxis()->SetTitleSize( 0.048 );
    hist->GetYaxis()->SetLabelSize( 0.042 );

    hist->GetXaxis()->SetTitleOffset( 0.95 );
    hist->GetXaxis()->SetTitleSize( 0.048 );
    hist->GetXaxis()->SetLabelSize( 0.042 );
}

void tuneHist2D( TH2D *hist )
{
    hist->GetYaxis()->SetTitleOffset( 1.45 );
    hist->GetYaxis()->SetTitleSize( 0.048 );
    hist->GetYaxis()->SetLabelSize( 0.042 );

    hist->GetXaxis()->SetTitleOffset( 0.95 );
    hist->GetXaxis()->SetTitleSize( 0.048 );
    hist->GetXaxis()->SetLabelSize( 0.042 );
}

void tuneCanvas(TCanvas *canvas)
{
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.11);
//    canvas->SetGridy();
}

void drawGraph(TGraphErrors *gr, int marker, int color, const char* drawOpt =  "P same", double markerSize = -1 )
{
    gr->SetMarkerStyle( marker );
    gr->SetMarkerColor( color );
    gr->SetLineColor( color );
    if (markerSize>0)
        gr->SetMarkerSize( markerSize );
    gr->DrawClone( drawOpt );
}

void drawTex(TLatex *tex, double fontSize = 0.045)
{
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextSize(fontSize);
    //    tex->SetLineWidth(2);
    tex->DrawClone();

}

void calcPointsRatio(TGraphErrors *gr, TGraphErrors *grDenom)
{
    double x, y;
    double x1, y1;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        grDenom->GetPoint(i,x1,y1);
        gr->SetPoint(i,x,y/y1);
    }
}

void getQuantiles(TH1D *h, const int nq, double *yq)
{
    cout << "nq=" << nq << endl;
    //    const Int_t nq = 10;
    //    const Int_t nshots = 10;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    //    Double_t yq[nq];  // array to contain the quantiles
    for (Int_t i=0;i<nq;i++)
        xq[i] = Float_t(i+1)/nq;

    h->GetQuantiles( nq, yq, xq );

    for (Int_t i=0;i<nq;i++)
        cout << yq[i] << ", ";
    cout << endl;

}

void rearrangeBoundaries( const int nq, double *multBounds, double *boundsMin, double *boundsMax )
{
    for (Int_t i=0;i<nq;i++)
    {
        if ( i == -1 )
            boundsMin[i] = 0;
        else
            boundsMin[i] = multBounds[i-1];
        boundsMax[i] = multBounds[i+1];
        cout << "bounds for i=" << i << ": " << boundsMin[i] << " " << boundsMax[i] << endl;
    }
}


void drawCanvasWithClasses(TH1D *hist1D, TString label, const int nCentrBins, double *centrMultBounds, double *multBinCenters )
{
    TCanvas *canv_mult_withClasses = new TCanvas(Form("canv_mult_withClasses_%s",label.Data())
                                                 ,Form("canv_mult_withClasses_%s",label.Data()),300,100,800,600 );
    tuneCanvas(canv_mult_withClasses);

    hist1D->SetLineColor(kRed);
    hist1D->DrawCopy();

    gPad->SetLogy();

    //draw centrality classes on mult hist
    //    int nBinsForClassHists = fHistParticlesInCutConditionInEvent->GetNbinsX();
    TH1D **fHistCentrClass = new TH1D*[nCentrBins];
    for ( int iCentrClass = 0; iCentrClass < nCentrBins; iCentrClass++ )
    {
        fHistCentrClass[iCentrClass] = (TH1D*)hist1D->Clone( Form("fHistCentrClass%d", iCentrClass) );   //new TH1D( Form("fHistCentrClass%d", iCentrClass), Form("iCentrClass%d", iCentrClass)
        //                                                 , hist1D->GetNbinsX(), hist1D->GetBinLowEdge(1), hist1D->GetBin );
        fHistCentrClass[iCentrClass]->Reset();
    }

    int centrClassId = 0;
    for ( int iBin = 0; iBin < hist1D->GetNbinsX(); iBin++ )
    {
        //        cout << fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) << " " << centrMultBounds[centrClassId] << endl;
        if ( hist1D->GetBinCenter(iBin+1) > centrMultBounds[centrClassId] )
            centrClassId++;
        if ( centrClassId >= nCentrBins )
            break;
        double binContent = hist1D->GetBinContent( iBin+1 );
        fHistCentrClass[centrClassId]->SetBinContent( iBin+1, binContent );
    }
    for ( int iCentrClass = 0; iCentrClass < nCentrBins; iCentrClass++ )
    {
        fHistCentrClass[iCentrClass]->SetFillColor( kOrange - 5 + iCentrClass );
        fHistCentrClass[iCentrClass]->SetLineColor(kBlue);
        fHistCentrClass[iCentrClass]->DrawCopy( "same" );
        multBinCenters[iCentrClass] = fHistCentrClass[iCentrClass]->GetMean();

        //                double meanNchInCentrBin = fHistCentrClass[iCentrClass]->GetMean();
        //                fHistMultClassMeanNch->SetBinContent( iCentrClass+1, meanNchInCentrBin );
        //                cout << "meanNch in centrality bin " << iCentrClass << ": " << meanNchInCentrBin << endl;
    }

    canv_mult_withClasses->SaveAs( Form("output/canv_%s_classes_%d.eps", label.Data(), nCentrBins));

}




void shiftPointX(TGraphErrors *gr, double shift)
{
    double x, y;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        gr->SetPoint(i,x+shift,y);
    }
}
