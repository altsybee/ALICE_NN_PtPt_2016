void run_FB_analysis( int i = -1 )
{
    gROOT->ProcessLine( ".L SupplementaryClasses.cxx+");
    gROOT->ProcessLine( ".L analyse_FB_TREE.C+");


    //    int LIST_ID_FOR_INEFF = 7;
    //    for (Int_t ineffId = 0; ineffId < 9; ineffId++ )
    //    {


//    int nEventsByHand[] = { 14e6 };
    //2e4, 5e4};//,
                           // 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 14e6 };

    //for (Int_t i = 0; i < 91; i++ )
        //for (Int_t i = 0; i < 1; i++ )
    {
        gROOT->Reset();

        TStopwatch timer;
        timer.Start();

        //        analyse_FB_TREE();// ineffId );
        //analyse_FB_TREE( i );//nEventsByHand[i] );

        //int fileIdByHand = -1
        //int nEventsByHand = -1, // )
        //LIST_ID_FOR_INEFF=-1
        analyse_FB_TREE();// -1, -1, 1 ); //i );//nEventsByHand[i] );


        //time estimation
        timer.Stop();
        Double_t rtime = timer.RealTime();
        Double_t ctime = timer.CpuTime();

        printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
    }



//        gROOT->ProcessLine( ".q");
}
