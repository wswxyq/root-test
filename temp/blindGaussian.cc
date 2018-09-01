#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
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
// Blinding the signal yield of a Gaussian fit.
using namespace std;
void blindGaussian()
{
  // Set a flag to blind or unblind, be cautious of this in practice you don't want 
  // to unblind yourself by mistake...
  Bool_t isBlind(kTRUE);
  //Bool_t isBlind(kFALSE);
   
  // Expected yield, this parameter sets the max number of events to be generated in
  // the Gaussian, hence when we fit we "expect" to see a number of this magnitude 
  // returned in our final result.
  double expectedYield(10000);
 
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
  RooDataSet* data = gaussian.generate(x, expectedYield) ;
   
  // Define the parameter containing the signal yields from the fit. Always allow floatation
  // below 0, in practice yields can be negative. I allow a range of 2 * the maximum expected, 
  // in general should stop the errors calculations reaching limits.
  RooRealVar nsig("nsig", "nsig", expectedYield/2.0, -expectedYield * 2.0, expectedYield * 2.0);
 
  // Create an unblinding transformation (a function object deriving from RooAbsReal).
  // This object will return the unblind value corresponding to the blind value in nsig.
  // In this example we use the 'precision blinding' technique as found in
  // RooBlindTools, which is a direct translation of the BlindTools package.
  // Please consult the BlindTools documentation for details on the blinding technique itself
  RooUnblindPrecision m_unblind_nsig("m_unblind_nsig","nsig (unblind)", "BlindString", nsig.getVal(), 100, nsig, blindCat, 0);
   
  // Build a PDF feeding it the unblinded value of m instead of m itself
  RooAbsPdf* extended = new RooExtendPdf("extened","extended", gaussian, m_unblind_nsig); 
   
  // Fit data with gfit (using a blind deltam)
  RooFitResult* fitResult = extended->fitTo(*data, RooFit::Extended(kTRUE), RooFit::NumCPU(2), RooFit::Save(kTRUE), RooFit::Minos(kTRUE)) ;
 
  // m_unblind will not reveal its contents easily
  m_unblind_nsig.Print() ;
  fitResult->Print("V");
  TCanvas* c1 = new TCanvas("c1", "c1", 900, 500);
  RooPlot* frame = x.frame( RooFit::Title("Best Plot Ever"));
  data->plotOn(frame, RooFit::Binning(40), RooFit::DataError( RooAbsData::SumW2 ) );
  extended->plotOn(frame, RooFit::LineColor(4), RooFit::Normalization( 1.0, RooAbsReal::Relative));
  frame->Draw();
  frame->SetMinimum(1.e-5); 
  c1->Draw();
  cout<<fitResult->status()<<endl;
 
  delete extended; extended=0;
}