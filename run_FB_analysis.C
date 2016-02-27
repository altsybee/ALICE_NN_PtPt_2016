void run_FB_analysis()
{
    gROOT->ProcessLine( ".L SupplementaryClasses.cxx+");
    gROOT->ProcessLine( ".L get_nn_ptpt_FROM_FB_TREE.C+");

    get_nn_ptpt_FROM_FB_TREE();

//    gROOT->ProcessLine( ".q");
}
