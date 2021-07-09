#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
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
// Faked data for a blinded fit to two crystal ball functions
#include <math.h>
#include <string>
#include <sstream>

using namespace RooFit;

//-------------------------------------------------------------------------------
//-------------------------- Configure parameters -------------------------------
// Master blinding tool, if false everything is unblinded
Bool_t isBlind(kTRUE);
// Define blinding Bs or Bd regions
Bool_t blindBs(kTRUE);
Bool_t blindBd(kTRUE);
Bool_t isResidual(kFALSE);


// Set fileName and paths
TString massVariable("B_M_JpsiConstr");
TString decayMode("B #rightarrow J/#psi K_{S}^{0}");
Int_t bins(40);
// Mass Ranges
Double_t B_mass_lowest(5167.0), B_mass_highest(5478.0),
         Bd_mass_min(5251.0), Bd_mass_max(5309.0),
         Bs_mass_min(5330.0), Bs_mass_max(5393.0);
// Estimates of expected yields for signal and background
Double_t E_Bd_nsig(700.0), E_Bs_nsig(400.0), E_background(5000.0);
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
//
// FUNCTION DEFINITION:
//
//-------------------------------------------------------------------------------

// Evaluate the residual:
double residual( double datum, double pdf )
{
    double chi2 = 0.;

    if ( pdf > 0 )
        chi2 += 2. * ( pdf - datum );

    if ( datum > 0 && pdf > 0 )
        chi2 += 2. * datum * log( datum / pdf );

    return ( ( datum >= pdf ) ? sqrt( chi2 ) : -sqrt( chi2 ) );
}

// Make the RooHist of the residual.
TH1D* residualHist( const RooHist* rhist, const RooCurve* curve )
{
    double r = 0.2;
    double sr = 1. / r;

    // Grab info from the histogram.
    int     n = rhist->GetN();
    double* x = rhist->GetX();
    double* y = rhist->GetY();

    // Create residual histogram.
    double xMin = x[ 0     ];
    double xMax = x[ n - 1 ];
    TH1D* residuals_temp = new TH1D( "r", "", n, xMin, xMax );

    double datum = 0.;
    double pdf   = 0.;

    // Fill the histogram.
    if ( curve )
        for ( int bin = 0; bin < n; bin++ )
        {
            datum = y[ bin ];
            pdf   = curve->Eval( x[ bin ] );
            residuals_temp->SetBinContent( bin + 1, residual( datum, pdf ) );
            residuals_temp->SetBinError  ( bin + 1, 1. );
        }

    residuals_temp->SetMinimum    ( -5.   );
    residuals_temp->SetMaximum    (  5.   );
    residuals_temp->SetStats      ( false );
    residuals_temp->SetMarkerStyle( 8     );
    residuals_temp->SetMarkerSize ( .8    );

    TAxis* xAxis = residuals_temp->GetXaxis();
    xAxis->SetTickLength ( sr * xAxis->GetTickLength()  );
    xAxis->SetLabelSize  ( sr * xAxis->GetLabelSize()   );
    xAxis->SetTitleSize  ( sr * xAxis->GetTitleSize()   );
    xAxis->SetLabelOffset( sr * xAxis->GetLabelOffset() );

    TAxis* yAxis = residuals_temp->GetYaxis();
    yAxis->SetNdivisions ( 504                          );
    yAxis->SetLabelSize  ( sr * yAxis->GetLabelSize()   );

    return residuals_temp;
}

// To draw lines in the residual histograms:
TLine* midLine;
TLine* uppLine;
TLine* lowLine;

void Lines( double xMin, double xMax )
{
    midLine = new TLine( xMin,  0., xMax,  0. );
    uppLine = new TLine( xMin,  2., xMax,  2. );
    lowLine = new TLine( xMin, -2., xMax, -2. );

    uppLine->SetLineColor( kRed );
    lowLine->SetLineColor( kRed );

    midLine->Draw( "same" );
    uppLine->Draw( "same" );
    lowLine->Draw( "same" );
}

void Lines( const TH1D* resid )
{
    double xMin = resid->GetXaxis()->GetXmin();
    double xMax = resid->GetXaxis()->GetXmax();

    midLine = new TLine( xMin,  0., xMax,  0. );
    uppLine = new TLine( xMin,  2., xMax,  2. );
    lowLine = new TLine( xMin, -2., xMax, -2. );

    uppLine->SetLineColor( kRed );
    lowLine->SetLineColor( kRed );

    midLine->Draw( "same" );
    uppLine->Draw( "same" );
    lowLine->Draw( "same" );
}

TString convertToString(Double_t number)
{
    std::ostringstream ss;
    ss <<  number;
    //return std::to_string(number);
    return TString(ss.str());
}

Double_t convertToDouble( const std::string& s )
{
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        return 0;
    return x;
}

//-------------------------------------------------------------------------------
//
// Main fitting sequence
//
//-------------------------------------------------------------------------------
void fit2CrystalBall()
{
    gROOT->SetStyle("Plain");

    if(isBlind && (!blindBs && !blindBd))
    {
        std::cout << "WARNING:"<< std::endl;
        std::cout << "    isBlind is true, but you have a region you wish to unblind."<< std::endl;
        std::cout << "    for clarity blindBd=" << blindBd << " and blindBs=" << blindBs << std::endl;
        return;
    }

    RooRealVar  *MASS   = new RooRealVar(massVariable, "M_{"+decayMode+"}", B_mass_lowest, B_mass_highest, "MeV/c^{2}");

    Double_t dmass(5.), Bs_mean(5366.77), Bd_mean(5279.58);
    RooRealVar *cbmean_core = new RooRealVar("cbmean_core", "cbmean core", Bd_mean, Bd_mean - dmass, Bd_mean + dmass) ;
    cbmean_core->setConstant(kTRUE);
    RooRealVar *cbsigma_core = new RooRealVar("cbsigma_core", "cbsigma core", 6, 1, 10.0) ;
    RooRealVar *n_core = new RooRealVar("n_core","", 1.15304);
    RooRealVar *alpha_core = new RooRealVar("alpha_core","", 2.32123);
    RooAbsPdf *cball_core = new RooCBShape("cball_core", "crystal ball core", *MASS, *cbmean_core, *cbsigma_core, *alpha_core, *n_core);

    RooRealVar *cbmean_tail = new RooRealVar("cbmean_tail", "cbmean tail", Bd_mean, Bd_mean - dmass, Bd_mean + dmass);
    cbmean_tail->setConstant(kTRUE);
    RooRealVar *cbsigma_tail = new RooRealVar("cbsigma_tail", "cbsigma tail", 5, 1, 10.0) ;
    RooRealVar *n_tail = new RooRealVar("n_tail","", 0.396915);
    RooRealVar *alpha_tail = new RooRealVar("alpha_tail","", -3.11431);
    RooAbsPdf *cball_tail = new RooCBShape("cball_tail", "crystal ball tail", *MASS, *cbmean_core, *cbsigma_tail, *alpha_tail, *n_tail);

    RooRealVar *fmass_core = new RooRealVar("fmass_core", "Bd mass core fraction", 0.9,0.,1.0);
    RooAddPdf  *Bmass_signal = new RooAddPdf("Bmass_signal","Bd mass core+tail gauss", RooArgList(*cball_core,*cball_tail), *fmass_core);

    RooRealVar *Bscbmean_core  = new RooRealVar("Bscbmean_core", "Bscbmean core", Bs_mean, Bs_mean - dmass, Bs_mean + dmass);
    Bscbmean_core->setConstant(kTRUE);
    RooRealVar *Bscbsigma_core = new RooRealVar("Bscbsigma_core", "Bscbsigma core", 5, 1, 10.0) ;
    RooRealVar *Bsn_core       = new RooRealVar("Bsn_core","", 1.15304);
    RooRealVar *Bsalpha_core   = new RooRealVar("Bsalpha_core","", 2.32123);
    RooAbsPdf  *Bscball_core   = new RooCBShape("Bscball_core", "Bs crystal ball core", *MASS, *Bscbmean_core, *Bscbsigma_core, *Bsalpha_core, *Bsn_core);

    RooRealVar *Bscbmean_tail  = new RooRealVar("Bscbmean_tail", "Bs cbmean tail", Bs_mean, Bs_mean - dmass, Bs_mean + dmass);
    Bscbmean_tail->setConstant(kTRUE);
    RooRealVar *Bscbsigma_tail = new RooRealVar("Bscbsigma_tail", "Bs cbsigma tail", 5, 1, 10.0) ;
    RooRealVar *Bsn_tail       = new RooRealVar("Bsn_tail","", 0.396915);
    RooRealVar *Bsalpha_tail   = new RooRealVar("Bsalpha_tail","", -3.11431);
    RooAbsPdf  *Bscball_tail   = new RooCBShape("Bscball_tail", "Bs crystal ball tail", *MASS, *Bscbmean_core, *Bscbsigma_tail, *Bsalpha_tail, *Bsn_tail);

    RooRealVar *Bsfmass_core  = new RooRealVar("Bsfmass_core", "Bd mass core fraction", 0.9,0.,1.0);
    RooAddPdf  *Bsmass_signal = new RooAddPdf("Bsmass_signal","Bd mass core+tail gauss", RooArgList(*Bscball_core,*Bscball_tail),*Bsfmass_core);

    RooRealVar c0("c0","c0 parameter", 0.4, -1., 1.) ;
    RooRealVar c1("c1","c1 parameter", -0.07, -1., 1.) ;
    RooAbsPdf *pdf_back = new RooChebychev("pdf_back","Background", *MASS, RooArgList(c0,c1));

    // Generate a data set from the pdfs
    RooDataSet* nsigdata = Bmass_signal->generate( *MASS, E_Bd_nsig) ;
    RooDataSet* Bsnsigdata = Bsmass_signal->generate( *MASS, E_Bs_nsig) ;
    RooDataSet* bkgdata = pdf_back->generate( *MASS, E_background);
    RooDataSet *SignalData = new RooDataSet("data","data", RooArgSet(*MASS));
    SignalData->append(*nsigdata);
    SignalData->append(*Bsnsigdata);
    SignalData->append(*bkgdata);

    // To stop fit hitting bounds make sure to float yields with large range AND negative
    Double_t Nentries(static_cast<Double_t>(SignalData->sumEntries()));
    RooRealVar nsig("nsig","B_{d}^{0} Signal Yield", E_Bd_nsig, -2. * Nentries, 2. * Nentries);
    RooRealVar Bsnsig("Bsnsig","B_{s}^{0} Signal Yield", E_Bs_nsig, -2. * Nentries, 2. * Nentries);
    RooRealVar nbkg("nbkg","Background Yield", E_background, -2. * Nentries, 2. * Nentries);

    /////////////////////////////////////
    ///// Blind Procedure
    /////////////////////////////////////

    // RooCategory used for blind/unblind switching.
    TString blind("blind"), unblind("unblind");
    RooCategory blindCatBd("blindCatBd","Bd blind state Category");
    RooCategory blindCatBs("blindCatBs","Bs blind state Category");
    blindCatBd.defineType(unblind, 0);
    blindCatBs.defineType(unblind, 0);
    blindCatBd.defineType(blind, 1);
    blindCatBs.defineType(blind, 1);
    if(isBlind && blindBd)
        blindCatBd.setLabel(blind);
    else
        blindCatBd.setLabel(unblind);
    if(isBlind && blindBs)
        blindCatBs.setLabel(blind);
    else
        blindCatBs.setLabel(unblind);

    // Unblinding transformation
    RooUnblindPrecision B_nsig("B_nsig","nsig B0 blind","nsigB0blinded", nsig.getVal(), 1000.,  nsig, blindCatBd);
    RooUnblindPrecision B_Bsnsig("B_Bsnsig","nsig Bs blind","nsigBsblinded", Bsnsig.getVal(), 1000., Bsnsig, blindCatBs);

    RooAbsPdf *ext_pdf_sig   = new RooExtendPdf("ext_pdf_sig","Ext. PDF Signal",*Bmass_signal, B_nsig);
    RooAbsPdf *ext_pdf_Bssig = new RooExtendPdf("ext_pdf_Bssig","Ext. PDF Signal",*Bsmass_signal, B_Bsnsig);
    RooAbsPdf *ext_pdf_bkg     = new RooExtendPdf("ext_pdf_bkg","Ext. PDF back",*pdf_back, nbkg);

    // The total extended pdf
    RooAbsPdf *pdf_total = new RooAddPdf("pdf_total","Total PDF", RooArgList(*ext_pdf_sig, *ext_pdf_Bssig, *ext_pdf_bkg));


    MASS->setRange("lowersideband", B_mass_lowest, Bd_mass_min);
    MASS->setRange("uppersideband", Bs_mass_max, B_mass_highest);
    MASS->setRange("unblind_signal_Bd", B_mass_lowest, Bs_mass_min);
    MASS->setRange("middleregion", Bd_mass_max, Bs_mass_min);
    MASS->setRange("unblind_signal_Bs", Bd_mass_max, B_mass_highest);
    TCut lowersideband = TCut(massVariable + "<" + convertToString(Bd_mass_min));
    TCut uppersideband = TCut(TString( massVariable + ">" + convertToString(Bs_mass_max) ) );
    TCut middleregion = TCut(TString( massVariable + "<" + convertToString(Bs_mass_min) + "&&" + massVariable + ">" + convertToString(Bd_mass_max) ) );
    TCut unblind_signal_Bd = TCut(TString( massVariable + "<" + convertToString(Bs_mass_min) ) );
    TCut unblind_signal_Bs = TCut(TString( massVariable + ">" + convertToString(Bd_mass_max) ) );

    RooFitResult* fitResult = pdf_total->fitTo(*SignalData, NumCPU(4), Extended(kTRUE), Save(kTRUE), Minos(kFALSE));
    RooMsgService::instance().setSilentMode(kTRUE);

    // Plot PDF and toy data overlaid
    TCanvas* canvas = new TCanvas("canvas","B_{d}^{0} Mass pdf",1000,700) ;

    if(isResidual) canvas->SetWindowSize(1000, 800);
    Double_t chi2(0.);

    Double_t nd1(0.), nd2(0.), nd3(0.);
    RooPlot *p3 = MASS->frame(Title("Invariant Mass Plot M_{"+decayMode+"}"));
    if(!isBlind) {
        SignalData->plotOn(p3, Binning(bins), Name("data_hist"), DataError( RooAbsData::SumW2 ));
        pdf_total->plotOn(p3, Range("fitRange"), Normalization(1.0,RooAbsReal::RelativeExpected), Name("main_curve"));
        pdf_total->plotOn(p3, Range("fitRange"), Components(*ext_pdf_bkg),LineStyle(kDashed),LineColor(kMagenta),Normalization(1.0,RooAbsReal::RelativeExpected));
        pdf_total->plotOn(p3, Range("fitRange"), Components(*ext_pdf_sig),LineStyle(kDashed),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
        pdf_total->plotOn(p3, Range("fitRange"), Components(*ext_pdf_Bssig),LineStyle(kDashed),LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
    }
    else
    {
        if(blindBd && blindBs)
        {
            std::cout << "Keeping both Bs and Bd regions blind!" << std::endl;
            nd1=(static_cast<Double_t>((SignalData->reduce(lowersideband)->numEntries())))/(static_cast<Double_t>(Nentries));
            nd2=(static_cast<Double_t>((SignalData->reduce(middleregion)->numEntries())))/(static_cast<Double_t>(Nentries));
            nd3=(static_cast<Double_t>((SignalData->reduce(uppersideband)->numEntries())))/(static_cast<Double_t>(Nentries));

            SignalData->plotOn(p3, CutRange("lowersideband"), Name("data_hist1"), Binning(bins) );
            pdf_total->plotOn( p3, Range("lowersideband"), Normalization(nd1,RooAbsReal::Relative), LineWidth(3), LineColor( kBlue), Name("main_curve1"));
            SignalData->plotOn(p3, CutRange("middleregion"), Binning(bins), DataError( RooAbsData::SumW2 ) );
            pdf_total->plotOn( p3, Range("middleregion"), Normalization(nd2,RooAbsReal::Relative), LineWidth(3), LineColor( kBlue));
            SignalData->plotOn( p3, CutRange("uppersideband"), Binning(bins), DataError( RooAbsData::SumW2 ) );
            pdf_total->plotOn( p3, Range("uppersideband"), Normalization(nd3,RooAbsReal::Relative), LineWidth(3), LineColor( kBlue));
        }
        else if(!blindBd && blindBs)
        {
            std::cout << "Unblinding the Bd mass region, keep Bs blinded!" << std::endl;
            nd1=(static_cast<Double_t>((SignalData->reduce(unblind_signal_Bd)->numEntries())))/(static_cast<Double_t>(Nentries));
            nd2=(static_cast<Double_t>((SignalData->reduce(uppersideband)->numEntries())))/(static_cast<Double_t>(Nentries));

            SignalData->plotOn(p3, CutRange("unblind_signal_Bd"), Name("data_hist1"), Binning(bins), DataError( RooAbsData::SumW2 ) );
            pdf_total->plotOn( p3, Range("unblind_signal_Bd"), Normalization(nd1, RooAbsReal::Relative), LineWidth(3), LineColor( kBlue), Name("main_curve1"));
            pdf_total->plotOn(p3, Range("unblind_signal_Bd"), Components(*ext_pdf_sig),LineStyle(kDashed),LineColor(kRed),Normalization(nd1,RooAbsReal::RelativeExpected));
            SignalData->plotOn( p3, CutRange("uppersideband"), Binning(bins), DataError( RooAbsData::SumW2 ) );
            pdf_total->plotOn( p3, Range("uppersideband"), Normalization(nd2, RooAbsReal::Relative), LineWidth(3), LineColor( kBlue));

        }
        else
        {
            std::cout << "Unblinding the Bs mass region, keep Bd blinded!" << std::endl;
            nd1=(static_cast<Double_t>((SignalData->reduce(lowersideband)->numEntries())))/(static_cast<Double_t>(Nentries));
            nd2=(static_cast<Double_t>((SignalData->reduce(unblind_signal_Bs)->numEntries())))/(static_cast<Double_t>(Nentries));

            SignalData->plotOn(p3, CutRange("lowersideband"), Name("data_hist1"), Binning(bins), DataError( RooAbsData::SumW2 ) );
            pdf_total->plotOn( p3, Range("lowersideband"), Normalization(nd1,RooAbsReal::Relative), LineWidth(3), LineColor( kBlue), Name("main_curve1"));
            SignalData->plotOn(p3, CutRange("unblind_signal_Bs"), Binning(bins), DataError( RooAbsData::SumW2 ) );
            pdf_total->plotOn( p3, Range("unblind_signal_Bs"), Normalization(nd2, RooAbsReal::Relative), LineWidth(3), LineColor( kBlue));
            pdf_total->plotOn(p3, Range("unblind_signal_Bs"), Components(*ext_pdf_Bssig),LineStyle(kDashed),LineColor(kGreen),Normalization(nd1,RooAbsReal::RelativeExpected));
        }
    }

    Float_t xmin(0.12), xmax(.37), ymin(.68), ymax(.88);

    if (isResidual && !isBlind)
    {
        p3->GetXaxis()->SetLabelSize(0);
        p3->GetXaxis()->SetTitleOffset(1.15);
        p3->GetYaxis()->SetTitleOffset(1.26);
        RooHist* histogram = p3->getHist("data_hist");
        RooCurve* curve = p3->getCurve("main_curve");
        TH1D* hresidual  = residualHist(histogram,curve);
        canvas->Divide( 1, 2, .1, .1 );
        double r  = .25;
        double sl = 1. / ( 1. - r );
        p3->GetYaxis()->SetLabelSize(sl * p3->GetYaxis()->GetLabelSize());
        //double labelSize = 0.2;
        hresidual->GetXaxis()->SetLabelSize((1./(1.+r)) * hresidual->GetXaxis()->GetLabelSize());
        hresidual->GetYaxis()->SetLabelSize((1./(1.+r)) * hresidual->GetYaxis()->GetLabelSize());
        TPad* padHisto = (TPad*) canvas->cd(1);
        TPad* padResid = (TPad*) canvas->cd(2);
        double small = 0.1;
        padHisto->SetPad( 0., r, 1., 1. );
        padHisto->SetBottomMargin( small );
        padResid->SetPad( 0., 0., 1., r  );
        padResid->SetBottomMargin( 0.3  );
        padResid->SetTopMargin   ( small );
        padHisto->cd();
        TPaveText *box= new TPaveText( xmin, ymin, xmax, ymax, "NDC");
        box->SetFillColor(10);
        box->SetBorderSize(0);
        box->SetTextAlign(12);
        box->SetTextSize(0.04F);
        box->SetFillStyle(1001);
        box->SetFillColor(10);
        Char_t buf[30], buf2[30], buf3[30], bufflhcb[30];
        sprintf( buf,  "B^{0} nsig = %3.2f", nsig.getVal() );
        sprintf( buf2,  "B_{s}^{0} nsig= %3.2f", Bsnsig.getVal() );
        sprintf( buf3,  "nbkg= %3.2f", nbkg.getVal() );
        TText* textbd = box->AddText( buf );
        textbd->SetTextSize(0.035);
        textbd->SetTextAlign(11);
        //textbd->SetTextColor(kRed);
        TText* textbs = box->AddText( buf2 );
        textbs->SetTextSize(0.035);
        textbs->SetTextAlign(11);
        //textbs->SetTextColor(kGreen);
        TText* textbkg = box->AddText( buf3 );
        textbkg->SetTextSize(0.035);
        textbkg->SetTextAlign(11);
        //textbkg->SetTextColor(kMagenta);
        p3->addObject(box);
        p3->Draw() ;
        //myPave->Draw("SAME");
        padResid->cd();
        hresidual->Draw();
        Lines( hresidual );
        hresidual->Draw( "SAME" );
    }
    else
    {
        p3->GetXaxis()->SetTitleOffset(1.15);
        p3->GetYaxis()->SetTitleOffset(1.26);
        gPad->SetLeftMargin(0.1);
        TPaveText *box= new TPaveText( xmin, ymin, xmax, ymax, "NDC");
        // Or set the exact coordinatres -> TPaveText *box= new TPaveText( 5180., 160., 5260., 200., "br");
        box->SetFillColor(10);
        box->SetBorderSize(0);
        box->SetTextAlign(12);
        box->SetTextSize(0.04F);
        box->SetFillStyle(1001);
        Char_t buf[30], buf2[30], buf3[30];
        sprintf( buf,  "B_{d}^{0} yield = %3.2f", nsig.getVal() );
        sprintf( buf2,  "B_{s}^{0} yield = %3.2f", Bsnsig.getVal() );
        sprintf( buf3,  "bkg yield = %3.2f", nbkg.getVal() );
        TText* textbd = box->AddText( buf );
        textbd->SetTextSize(0.035);
        textbd->SetTextAlign(11);
        //textbd->SetTextColor(kRed);
        TText* textbs = box->AddText( buf2 );
        textbs->SetTextSize(0.035);
        textbs->SetTextAlign(11);
        //textbs->SetTextColor(kGreen);
        TText* textbkg = box->AddText( buf3 );
        textbkg->SetTextSize(0.035);
        textbkg->SetTextAlign(11);
        //textbkg->SetTextColor(kMagenta);
        p3->addObject(box);
        p3->Draw() ;
    }

    p3->SetMinimum(1.e-5);
    //canvas->Modify();
    canvas->Draw();
    canvas->SaveAs("MassPlot.eps");
    canvas->SaveAs("MassPlot.png");
    fitResult->Print("V");

    // to get a chi2 value we need to plot the complete pdf, for this as long as we never draw it we are ok.
    RooPlot *p = MASS->frame();
    SignalData->plotOn(p, Binning(bins), Name("data_chi2"), DataError( RooAbsData::SumW2 ));
    pdf_total->plotOn(p, Range("fitRange"), Normalization(1.0,RooAbsReal::RelativeExpected), Name("curve_chi2"));
    // There is a method within RooPlot called chiSquare(). This returns the reduced chi2 which when multiplied out
    // by the number of bins gives the actual chi2 fit between the .
    chi2 = p->chiSquare( "curve_chi2", "data_chi2") * bins;
    // The total number of degrees of freedom is defined as the number of bins subtracted by the free floating parameters
    // in the fit. We already know the number of bins as we specified it, so the free params (we also know) but can be
    // obtained using:-
    Int_t floated = static_cast<Int_t>(fitResult->floatParsFinal().getSize());

    std::cout << "Chi2 = " << chi2 << "  (for "<< bins << " bins and " << floated << " floating params)" << std::endl;
    std::cout<<" Prob(chi2,ndof) of above = " << TMath::Prob(chi2, bins-floated) << std::endl;
    return;
}