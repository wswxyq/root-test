#include <iostream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TChain.h"

using namespace std;

int mejn() {


///////////////////////////////////

    //Load a file

//////////////////////////////////

    TChain *chain= new TChain("DecayTreeTuple/DecayTree");
    chain->Add("Tuples/DVNtuple_mumu_tmp_MU.root");
    chain->Add("Tuples/DVNtuple_mumu_tmp_MD.root");
    //you can add as many files as needed

    cout<<"Number of entries: "<< chain->GetEntries()<<endl;
    //prints in the command line the number of entries

//////////////////////////////////
    //Step 1: Draw a histogram of the D0 mass and save it to a pdf
//////////////////////////////////

    TH1F *h1 = new TH1F("h1", "D0 mass", 50, 1700, 2000);
    chain->Project("h1", "D0_M") ;
    //this fills the histogram h1 with variable D0_M

    TCanvas *c1 = new TCanvas( "c1", "c1", 200, 10, 700, 500 );
    h1->Draw();
    c1->Print("test_MD.pdf");

    return 0;
}