#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooHist.h"
#include "TMath.h"
void upd() { 
    const char* TChain_name="D2KSPiPiPi_DD/DecayTree";
	const char* TChain_file="/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root";
    TFile *f = new TFile("hs.root","update"); 
    TTree *T = (TTree*)f->Get("ntuple"); 
    float px,py; 
    float pt; 
    TBranch *bpt = T->Branch("pt",&pt,"pt/F"); 
    T->SetBranchAddress("px",&px); 
    T->SetBranchAddress("py",&py); 
    Long64_t nentries = T->GetEntries(); 
    for (Long64_t i=0;i<nentries;i++) 
    { 
        T->GetEntry(i); 
        pt = TMath::Sqrt(px*px+py*py); 
        bpt->Fill();
     } 
    T->Print(); 
    T->Write(); 
    delete f;
    }