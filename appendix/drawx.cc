#include <iostream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TChain.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"

using namespace std;

void drawx() {
    //gStyle->SetOptStat(mode);
    const char* Tree_before_cut="ReducedTree";
    const char* Tree_after_cut="MVATree";
    const char* file_DD_before="../exec-all/DD/D2KSKpPiPi_DD.root";
    const char* file_LL_before="../exec-all/LL/D2KSKpPiPi_LL.root";
    const char* file_DD_after="../exec-all/DD/TMVA_output_DD.root";
    const char* file_LL_after="../exec-all/LL/TMVA_output_LL.root";
    TChain *chain_DD_before= new TChain(Tree_before_cut);
    TChain *chain_DD_after= new TChain(Tree_after_cut);
    TChain *chain_LL_before= new TChain(Tree_before_cut);
    TChain *chain_LL_after= new TChain(Tree_after_cut);
    chain_DD_before->Add(file_DD_before);
    chain_DD_after->Add(file_LL_after);
    chain_LL_before->Add(file_LL_before);
    chain_LL_after->Add(file_LL_after);

    Double_t low=1799.65;
    Double_t up=1939.65;

    TH1F *h1 = new TH1F("h1", "D_M", 50, low, up);
    TH1F *h2 = new TH1F("h2", "D_M", 50, low, up);
    TH1F *h3 = new TH1F("h3", "D_M", 50, low, up);
    TH1F *h4 = new TH1F("h4", "D_M", 50, low, up);


    //this fills the histogram h*_D_PT with variable D_PT

    TCanvas *c1 = new TCanvas( "c1", "D_M", 200, 10, 700, 500 );
    chain_DD_before->SetFillColor(8);
    chain_DD_after->SetFillColor(4);
    chain_DD_before->SetFillStyle(3005);
    chain_DD_after->SetFillStyle(3016);
    chain_DD_before->Draw("D_M","");
    chain_DD_after->Draw("D_M","","same");
    c1->Print("DD.png");


    TCanvas *c2 = new TCanvas( "c2", "D_M", 200, 10, 700, 500 );
    chain_LL_before->SetFillColor(8);
    chain_LL_after->SetFillColor(4);
    chain_LL_before->SetFillStyle(3005);
    chain_LL_after->SetFillStyle(3016);
    chain_LL_before->Draw("D_M","");
    chain_LL_after->Draw("D_M","","same");
    c2->Print("LL.png");


    delete h1;
    delete h2;
    delete h3;
    delete h4;




}
