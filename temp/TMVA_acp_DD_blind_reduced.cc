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
#include "RooUnblindPrecision.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooHist.h"
#include "TMath.h"
using namespace RooFit ;
using namespace std;

void TMVA_acp_DD_blind_reduced(){
	
	Bool_t isBlind(kTRUE);
	//Bool_t isBlind(kFALSE);
	RooCategory blindCat("blindCat","blind state Category");
	TString blind("blind"), unblind("unblind");
	blindCat.defineType(unblind, 0);
	blindCat.defineType(blind, 1);
   
	if(isBlind)
		blindCat.setLabel(blind);
	else
		blindCat.setLabel(unblind);
    //define TChain, TFile 
    const char* TChain_name="MVATree";
    const char* TChain_file="../saved/TMVA_root_files/TMVApp_DDx.root";
    TChain *chain= new TChain(TChain_name);
    chain->Add(TChain_file);



	RooRealVar D_M("D_M", "D+ mass", 1779.65, 1940.65, "MeV") ;
    RooRealVar BDT_response("BDT_response", "BDT_response value", 0.243397, RooNumber::infinity());
    RooRealVar Ct("Ct", "Ct value", -RooNumber::infinity(), RooNumber::infinity());
    RooRealVar P0_ID("P0_ID", "P0_ID value", -1000., 1000.);
	RooArgSet variables(D_M, BDT_response, Ct, P0_ID); 

	Double_t D_M_double ;
    Double_t BDT_response_double ;
    chain->SetBranchAddress("D_M",&D_M_double); 
    chain->SetBranchAddress("BDT_response",&BDT_response_double); 

	//creat dataset and write
    RooDataSet *ds_plus_plus=new RooDataSet("ds_plus_plus", "ds_plus_plus", variables, Import(*chain), 
		Cut(" BDT_response > 0.2433 && P0_ID>0 && Ct>0 "));

    
    RooRealVar mean("mean","mean of gaussian", 1870.15) ;
	RooRealVar sigma("sigma","width of gaussian", 7.50013) ;
	RooRealVar alpha_var("alpha_var","alpha_var", 1.36544) ;
	RooRealVar n_var("n_var","n_var", 42.2671) ;
	RooCBShape gauss_all("CBfunction", "CBfunction", D_M, mean, sigma, alpha_var, n_var);

	RooRealVar x1("x1", "para1", 0.448623);
	RooRealVar x2("x2", "para2", -0.102754);
	RooRealVar x3("x3", "para3", -0.0502192);
	RooChebychev poly("poly", "3-para-poly", D_M, RooArgList(x1, x2, x3));
	//num total : 94333


	//signal------------------------------------------------------------------
	RooRealVar Ntot_plus("Ntot_plus", "Ntot_plus", 0., 2000.); //number of D+

	RooRealVar AT_plus("AT_plus", "AT_plus", 0, -0.2, 0.2);
	RooUnblindPrecision AT_plus_unblind("AT_plus", "AT_plus", "Atblind", 0., 0.02, AT_plus, blindCat, 1);
	RooRealVar AT_minus("AT_minus", "AT_minus", 0, -0.2, 0.2);
	RooUnblindPrecision AT_minus_unblind("AT_minus", "AT_minus", "Atbarblind", 0., 0.02, AT_minus, blindCat, 1);

	RooFormulaVar Acp_unblind("Acp", "0.5*(@0-@1)", RooArgList(AT_plus, AT_minus));

    RooFormulaVar Nsig_plus_plus("Nsig_plus_plus", "0.5*@0*(1+@1)", RooArgList(Ntot_plus, AT_plus_unblind));

	//background----------------------------------------------------------
	RooRealVar Ntot_bkg_plus("Ntot_bkg_plus", "Ntot_bkg_plus", 0., 2000.);

	RooRealVar AT_bkg_plus("AT_bkg_plus", "AT_bkg_plus", -0.2, 0.2);

	RooFormulaVar Nbkg_plus_plus("Nbkg_plus_plus", "0.5*@0*(1+@1)", RooArgList(Ntot_bkg_plus, AT_bkg_plus));

	RooAddPdf  model_plus_plus("model_plus_plus","model_plus_plus",
		RooArgList(gauss_all, poly), RooArgList(Nsig_plus_plus, Nbkg_plus_plus)) ;

	//---define simultaneous PDF----


	model_plus_plus.fitTo(*ds_plus_plus);


	Ntot_plus.Print();  Ntot_bkg_plus.Print(); 
	AT_plus.Print(); AT_minus.Print(); AT_bkg_plus.Print(); 
	cout<<"============================="<<endl;
	cout<<AT_plus.getVal()<<"	"<<AT_minus.getVal()<<endl;
	cout<<AT_plus_unblind.getVal()<<"	"<<AT_minus_unblind.getVal()<<endl;
	AT_plus_unblind.Print(); AT_minus_unblind.Print();


	Double_t at_plus_double=AT_plus.getVal();
	Double_t at_minus_double=AT_minus.getVal();

	RooFormulaVar Acp("Acp", "0.5*(@0-@1)", RooArgList(AT_plus, AT_minus));


	cout<<Acp_unblind.getVal()<<endl;

}
