#include <iostream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TChain.h"
#include "TCut.h"

using namespace std;

int subtr_common(const char* var_name, const char* var_name_display,
                 const char* TSubfolder_name, const char* TChain_name, const char* TChain_file, const char* trigger_content,
                 const char* SignalWindow_content, const char* SideBand_content, Double_t low, Double_t up, Int_t o_index
                ) {

///////////////////////////////////

    //Load a file
    const char* bkg_file="h_bkg_";
    const char* sgnl_file="h_sgnl_";
    const char* bkg_sgnl_file="h_bkg&sgnl_";

    string bkg_fullPath = "../plotdir/";
    string sgnl_fullPath = "../plotdir/";
    string bkg_sgnl_fullPath = "../plotdir/";

    bkg_fullPath+=TSubfolder_name;
    bkg_fullPath+="/background/";
    bkg_fullPath+=bkg_file;
    bkg_fullPath+=var_name_display;
    bkg_fullPath+=".pdf";

    sgnl_fullPath+=TSubfolder_name;
    sgnl_fullPath+="/signal/";
    sgnl_fullPath+=sgnl_file;
    sgnl_fullPath+=var_name_display;
    sgnl_fullPath+=".pdf";

    bkg_sgnl_fullPath+=TSubfolder_name;
    bkg_sgnl_fullPath+="/sig_vs_bkg/";
    bkg_sgnl_fullPath+=bkg_sgnl_file;
    bkg_sgnl_fullPath+=var_name_display;
    bkg_sgnl_fullPath+=".pdf";

    const char* xbkg_file= bkg_fullPath.c_str();
    const char* xsgnl_file= sgnl_fullPath.c_str();
    const char* xbkg_sgnl_file = bkg_sgnl_fullPath.c_str();

    //////////////////////////////////
    TChain *chain= new TChain(TChain_name);
    chain->Add(TChain_file);
    //you can add as many files as needed
    cout<<"Number of entries: "<< chain->GetEntries()<<endl;
    //prints in the command line the number of entries

    TCut trigger = trigger_content;
    TCut SignalWindow = SignalWindow_content;
    TCut SideBand = SideBand_content;
    //TH1F *h1 = new TH1F("h1", var_name, 100, low, up);
    TH1F *h2 = new TH1F("h2", var_name, 100, low, up);
    TH1F *h3 = new TH1F("h3", var_name, 100, low, up);

    chain->Project("h1", var_name, trigger&&SignalWindow) ;
    chain->Project("h2", var_name, trigger&&SideBand) ;
    chain->Project("h3", var_name, trigger&&SignalWindow) ;
    //this fills the histogram h*_D_PT with variable D_PT

    TCanvas *c1 = new TCanvas( "c1", var_name, 200, 10, 700, 500 );
    h2->SetFillColor(8);
    h2->SetFillStyle(3005);
    h2->Draw();
    c1->Print(xbkg_file);

    TCanvas *c2 = new TCanvas( "c2", var_name, 200, 10, 700, 500 );
    h3->Add(h2, -1);
    h3->SetFillColor(4);
    h3->SetFillStyle(3016);
    h3->Draw();
    c2->Print(xsgnl_file);

    TCanvas *c3 = new TCanvas( "c3", var_name, 200, 10, 700, 500 );

    if (o_index==0) {
        h2->DrawNormalized();
        h3->DrawNormalized("same");
    } else
    {
        h3->DrawNormalized();
        h2->DrawNormalized("same");
    }
    c3->Print(xbkg_sgnl_file);

    delete chain;
    delete h2;
    delete h3;
    delete c1;
    delete c2;
    delete c3;

    return 0;
}