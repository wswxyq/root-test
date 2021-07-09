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

using namespace RooFit ;
using namespace std;
void ll() {
    TFile f("tree.root");
    TTree *datatree = (TTree*)(f.Get("MyTree"));
    TBranch *branch = datatree->GetBranch("B_M");
    int numEntry = branch->GetEntries();
    Float_t mass[100];
    branch->SetAddress(mass);
    RooRealVar m("B_M","mass",5300., 5450.);
    RooArgSet m_arg(m,"m_args");
    RooDataSet * dataset = new RooDataSet("dataset","dataset",m_arg);
    for(int i=0; i<datatree->GetEntries(); i++) {
        datatree->GetEvent(i);
        m = mass[0];
        dataset->add(m_arg,1.0);
    }
}