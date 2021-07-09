#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TH1D.h"
#include "RooRealVar.h"
#include "RooFit.h"
#include "RooBlindTools.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooUnblindPrecision.h"
#include "RooFitResult.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooHist.h"
#include "TMath.h"
#include "RooExtendPdf.h"
using namespace std;
Double_t gaussian_blind_fit(Bool_t isBlind) {

    // Set a flag to blind or unblind, be cautious of this in practice you don't want
    // to unblind yourself by mistake...
    //Bool_t isBlind(kFALSE);

    // Expected yield, this parameter sets the max number of events to be generated in
    // the Gaussian, hence when we fit we "expect" to see a number of this magnitude
    // returned in our final result.
    Double_t Yield_value(10);
    RooRealVar Yield("Yield", "Yield", Yield_value);

    // First declare the regular variables
    RooRealVar x("x","x", -10.0, 10.0) ;
    RooRealVar sigma("sigma","sigma", 1.0, 0.1, 2.0) ;
    RooRealVar mean("mean","mean", 0.0) ;
    mean.setConstant();

    // RooCategory used for blind/unblind switching.
    TString blind("blind"), unblind("unblind");
    RooCategory blindCat("blindCat","blind state Category");
    blindCat.defineType(unblind, 0);
    blindCat.defineType(blind, 1);

    if(isBlind)
        blindCat.setLabel(blind);
    else
        blindCat.setLabel(unblind);

    // Generate a dataset with a true m of 0.70
    RooGaussian gaussian("gaussian","gaussian", x, mean, sigma) ;
    RooDataSet* data = gaussian.generate(x, Yield_value) ;

    RooUnblindPrecision m_unblind("m_unblind","m (unblind)", "BlindString", 0., 0.1, x, blindCat, 0);
    RooGaussian testgaussian("testgaussian", "testgaussian", m_unblind, mean, sigma) ;

    // Build a PDF feeding it the unblinded value of m instead of m itself
    RooAbsPdf* extended = new RooExtendPdf("extened","extended", testgaussian, Yield);

    // Fit data with gfit (using a blind deltam)
    RooFitResult* fitResult = extended->fitTo(*data, RooFit::Extended(kTRUE), RooFit::NumCPU(2), RooFit::Save(kTRUE), RooFit::Minos(kTRUE)) ;

    // m_unblind will not reveal its contents easily
    m_unblind.Print() ;
    fitResult->Print("V");
    Int_t fit_index= fitResult->status();

    TCanvas* c1 = new TCanvas("c1", "c1", 900, 500);
    RooPlot* frame = x.frame( RooFit::Title("Best Plot Ever"));
    data->plotOn(frame, RooFit::Binning(40), RooFit::DataError( RooAbsData::SumW2 ) );
    extended->plotOn(frame, RooFit::LineColor(4), RooFit::Normalization( 1.0, RooAbsReal::Relative));
    frame->Draw();
    frame->SetMinimum(1.e-5);
    c1->Draw();
    cout<<"======================="<<endl;
    cout<<fit_index<<endl;

    delete extended;
    extended=0;
    return m_unblind.getVal();

}

void testblind() {
    double compare1, compare2;
    compare1=gaussian_blind_fit(kTRUE);
    compare2=gaussian_blind_fit(kFALSE);
    cout<<compare1<<endl;
    cout<<compare2<<endl;

}