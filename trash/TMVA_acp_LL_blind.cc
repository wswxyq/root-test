#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <fstream>
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
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
using namespace std ;

void TMVA_acp_LL_blind(){
	
	int i ;

    cout << "Using Blind Analysis? ('1' for yes, others for no)"<<"\n" ;
    cin >> i ;
	Bool_t isBlind ;

	if (i==1) {
		isBlind=kTRUE;
    }
	else
	{
		isBlind=kFALSE;
	}
	
	//Bool_t isBlind(kTRUE);
	//Bool_t isBlind(kFALSE);
	RooCategory blindCat("blindCat","blind state Category");
	TString blind("blind"), unblind("unblind");
	blindCat.defineType(unblind, 0);
	blindCat.defineType(blind, 1);
   
	if(isBlind)
		blindCat.setLabel(blind);
	else
		blindCat.setLabel(unblind);
	//blindCat.setIndex(int(isBlind));

    //define TChain, TFile 
    const char* TChain_name="MVATree";
    const char* TChain_file="./LL/TMVA_output_LL.root";
    TChain *chain= new TChain(TChain_name);
    chain->Add(TChain_file);


	RooRealVar D_M("D_M", "D+ mass", 1799.65, 1939.65, "MeV") ;
    RooRealVar BDT_response("BDT_response", "BDT_response value", 0.243397, RooNumber::infinity());
    RooRealVar Ct("Ct", "Ct value", -RooNumber::infinity(), RooNumber::infinity());
    RooRealVar P0_ID("P0_ID", "P0_ID value", -1000., 1000.);
	RooArgSet variables(D_M, BDT_response, Ct, P0_ID); 

	//Double_t D_M_double ;
    //Double_t BDT_response_double ;
    //chain->SetBranchAddress("D_M", &D_M_double); 
    //chain->SetBranchAddress("BDT_response", &BDT_response_double); 

	//creat dataset and write
    RooDataSet *ds_plus_plus=new RooDataSet("ds_plus_plus", "ds_plus_plus", variables, Import(*chain), 
		Cut(" P0_ID>0 && Ct>0 "));
	RooDataSet *ds_plus_minus=new RooDataSet("ds_plus_minus", "ds_plus_minus", variables, Import(*chain), 
		Cut(" P0_ID>0 && Ct<0 "));
	RooDataSet *ds_minus_plus=new RooDataSet("ds_minus_plus", "ds_minus_plus", variables, Import(*chain), 
		Cut(" P0_ID<0 && Ct>0 "));
	RooDataSet *ds_minus_minus=new RooDataSet("ds_minus_minus", "ds_minus_minus", variables, Import(*chain), 
		Cut(" P0_ID<0 && Ct<0 "));

    
	RooRealVar mean("mean","mean of gaussian",1870.21) ;
	RooRealVar mean_1("mean_1","mean of gaussian(1)",1864.91) ;
	RooRealVar sigma("sigma","width of gaussian",6.72256) ;
	RooRealVar sigma_1("sigma_1","width of gaussian(1)",26.1125) ;
	RooGaussian gauss("gauss","gaussian PDF", D_M, mean, sigma) ;  
	RooGaussian gauss_1("gauss_1","gaussian PDF(1)", D_M, mean_1, sigma_1) ;  
	RooRealVar gauss_frac("gauss_frac","fraction of component 1 in signal",0.782636) ;
	RooAddPdf gauss_all("gauss_all", "combination of gaussian function", RooArgList(gauss, gauss_1), RooArgList(gauss_frac));


	RooRealVar x1("x1", "para1", 0.00420437);
	//RooRealVar x2("x2", "para2", -1., 1.);
	//RooRealVar x3("x3", "para3", -1., 1.);
	//RooChebychev poly("poly", "3-para-poly", D_M, RooArgList(x1, x2, x3));
	//RooPolynomial poly("poly", "3-para-poly", D_M, RooArgList(x1));
	RooExponential poly("poly", "3-para-poly", D_M, x1);


	//signal------------------------------------------------------------------
	RooRealVar Ntot_plus("Ntot_plus", "Ntot_plus", 0., 200000.); //number of D+
	RooRealVar Ntot_minus("Ntot_minus", "Ntot_minus", 0., 200000.); //number od D-

	RooRealVar AT_plus("AT_plus", "AT_plus", 0, -0.1, 0.1);
	RooUnblindPrecision AT_plus_unblind("AT_plus", "AT_plus", "Atblind", 0., 0.1, AT_plus, blindCat, 1);
	RooRealVar AT_minus("AT_minus", "AT_minus", 0, -0.1, 0.1);
	RooUnblindPrecision AT_minus_unblind("AT_minus", "AT_minus", "Atbarblind", 0., 0.1, AT_minus, blindCat, 1);


    RooFormulaVar Nsig_plus_plus("Nsig_plus_plus", "0.5*@0*(1+@1)", RooArgList(Ntot_plus, AT_plus_unblind));
    RooFormulaVar Nsig_plus_minus("Nsig_plus_minus", "0.5*@0*(1-@1)", RooArgList(Ntot_plus, AT_plus_unblind));
    RooFormulaVar Nsig_minus_minus("Nsig_minus_minus", "0.5*@0*(1+@1)", RooArgList(Ntot_minus, AT_minus_unblind));
	RooFormulaVar Nsig_minus_plus("Nsig_minus_plus", "0.5*@0*(1-@1)", RooArgList(Ntot_minus, AT_minus_unblind));


	//background----------------------------------------------------------
	RooRealVar Ntot_bkg_plus("Ntot_bkg_plus", "Ntot_bkg_plus", 0., 200000.);
	RooRealVar Ntot_bkg_minus("Ntot_bkg_minus", "Ntot_bkg_minus", 0., 200000.);

	RooRealVar AT_bkg_plus("AT_bkg_plus", "AT_bkg_plus", -0.1, 0.1);
	RooRealVar AT_bkg_minus("AT_bkg_minus", "AT_bkg_minus", -0.1, 0.1);

	RooFormulaVar Nbkg_plus_plus("Nbkg_plus_plus", "0.5*@0*(1+@1)", RooArgList(Ntot_bkg_plus, AT_bkg_plus));
    RooFormulaVar Nbkg_plus_minus("Nbkg_plus_minus", "0.5*@0*(1-@1)", RooArgList(Ntot_bkg_plus, AT_bkg_plus));
    RooFormulaVar Nbkg_minus_minus("Nbkg_minus_minus", "0.5*@0*(1+@1)", RooArgList(Ntot_bkg_minus, AT_bkg_minus));
	RooFormulaVar Nbkg_minus_plus("Nbkg_minus_plus", "0.5*@0*(1-@1)", RooArgList(Ntot_bkg_minus, AT_bkg_minus));


	RooAddPdf  model_plus_plus("model_plus_plus","model_plus_plus",
		RooArgList(gauss_all, poly), RooArgList(Nsig_plus_plus, Nbkg_plus_plus)) ;
	RooAddPdf  model_minus_plus("model_minus_plus","model_minus_plus",
		RooArgList(gauss_all, poly), RooArgList(Nsig_minus_plus, Nbkg_minus_plus)) ;
	RooAddPdf  model_plus_minus("model_plus_minus","model_plus_minus",
		RooArgList(gauss_all, poly), RooArgList(Nsig_plus_minus, Nbkg_plus_minus)) ;
	RooAddPdf  model_minus_minus("model_minus_minus","model_minus_minus",
		RooArgList(gauss_all, poly), RooArgList(Nsig_minus_minus, Nbkg_minus_minus)) ;



	//---create combined datasample for simultaneous fit----
	RooCategory alldata("alldata","alldata") ;
	alldata.defineType("plus_plus") ;
	alldata.defineType("plus_minus") ;
	alldata.defineType("minus_plus") ;
	alldata.defineType("minus_minus") ;

  	RooDataSet combData("combData","combined data", variables, Index(alldata),
		Import("plus_plus",*ds_plus_plus), Import("plus_minus",*ds_plus_minus), 
		Import("minus_plus",*ds_minus_plus), Import("minus_minus",*ds_minus_minus)) ;

	//here, dh1 and dh1m are the two datasets and mass is the variable in which you fit

	//---define simultaneous PDF----
	RooSimultaneous simPdf("simPdf","simultaneous pdf", alldata) ;
	simPdf.addPdf(model_plus_plus,"plus_plus") ;
	simPdf.addPdf(model_plus_minus,"plus_minus") ;
	simPdf.addPdf(model_minus_plus,"minus_plus") ;
	simPdf.addPdf(model_minus_minus,"minus_minus") ;
	cout<<"begin fit"<<endl;
	//---do the sim. fit---
	auto result= simPdf.fitTo(combData, RooFit::NumCPU(4), RooFit::Save(kTRUE), RooFit::Minos(kTRUE)) ;
	Int_t fit_index = result->status();
	//model_plus_plus.fitTo(*ds_plus_plus);

	cout<<"end fit"<<endl;

	RooFormulaVar Acp_blind("Acp", "0.5*(@0-@1)", RooArgList(AT_plus_unblind, AT_minus_unblind));

	//---draw a component----
	TCanvas* c = new TCanvas("total_plot","total_plot", 1800, 1200) ;
	c->Divide(2,2) ;


	RooPlot* frame1 = D_M.frame(Title("D+(Ct>0)"));
	combData.plotOn(frame1,Cut("alldata==alldata::plus_plus")) ;
	simPdf.plotOn(frame1,Slice(alldata,"plus_plus"),ProjWData(alldata, combData)) ;
	//simPdf.paramOn(frame1);
	combData.statOn(frame1);
	c->cd(1) ; 	
	frame1->Draw();

	RooPlot* frame2 = D_M.frame(Title("D+(Ct<0)"));
	combData.plotOn(frame2,Cut("alldata==alldata::plus_minus")) ;
	simPdf.plotOn(frame2,Slice(alldata,"plus_minus"),ProjWData(alldata, combData)) ;
	//simPdf.paramOn(frame2);
	combData.statOn(frame2);
	c->cd(2) ; 	
	frame2->Draw();

	RooPlot* frame3 = D_M.frame(Title("D-(Ct>0)"));
	combData.plotOn(frame3,Cut("alldata==alldata::minus_plus")) ;
	simPdf.plotOn(frame3,Slice(alldata,"minus_plus"),ProjWData(alldata, combData)) ;
	//simPdf.paramOn(frame3);
	combData.statOn(frame3);
	c->cd(3) ; 	
	frame3->Draw();

	RooPlot* frame4 = D_M.frame(Title("D-(Ct<0)"));
	combData.plotOn(frame4,Cut("alldata==alldata::minus_minus")) ;
	simPdf.plotOn(frame4,Slice(alldata,"minus_minus"),ProjWData(alldata, combData)) ;
	//simPdf.paramOn(frame4);
	combData.statOn(frame4);
	c->cd(4) ; 	
	frame4->Draw();
	c->Print("./LL/At.eps");
	c->Print("./LL/At.pdf");

	Ntot_plus.Print(); Ntot_minus.Print(); Ntot_bkg_plus.Print(); Ntot_bkg_minus.Print();
	AT_plus.Print(); AT_minus.Print(); AT_bkg_plus.Print(); AT_bkg_minus.Print();
	AT_plus_unblind.Print(); AT_minus_unblind.Print();
	//result->Print("V");


	Double_t at_plus_double=AT_plus.getVal();
	Double_t at_minus_double=AT_minus.getVal();

	RooFormulaVar Acp("Acp", "0.5*(@0-@1)", RooArgList(AT_plus, AT_minus));
	cout<< "Acp = " << 0.5*(at_plus_double-at_minus_double) <<" +/- "
		<<sqrt(TMath::Power(0.5*AT_plus.getError(), 2)+pow(0.5*AT_minus.getError(), 2))<<endl;

	//cout<< "Acp_unblind = " << Acp_unblind.getVal() <<" +/- "
	//	<< Acp_unblind.getPropagatedError(*result)<<endl;
	//Acp.getPropagatedError(*result);


	std::ofstream myfile;
    myfile.open ("./LL/fitresult.txt");
	myfile<<"Fitted number of events: "<< 
		Ntot_plus.getVal()+Ntot_minus.getVal()+Ntot_bkg_plus.getVal()+Ntot_bkg_minus.getVal() <<endl;
	myfile<<"Number of total entries in dataset: "<<ds_plus_plus->numEntries()+
	ds_plus_minus->numEntries()+ds_minus_plus->numEntries()+ds_minus_minus->numEntries()<<endl;
	myfile<<"signal yield:"<< Ntot_plus.getVal()+Ntot_minus.getVal()<<endl;
	myfile<<"background yield:"<< Ntot_bkg_plus.getVal()+Ntot_bkg_minus.getVal()<<endl;
	myfile<<"============================="<<endl;
	myfile<<"At(unblind):"<<AT_plus.getVal()<<endl;
	myfile<<"At_bar(unblind):"<<AT_minus.getVal()<<endl;
	myfile<<"At(blind):"<<AT_plus_unblind.getVal()<<endl;
	myfile<<"At_bar(blind):"<<AT_minus_unblind.getVal()<<endl;
	myfile<<"Acp_blind: "<< Acp_blind.getVal() << " +/- "
		<<sqrt(TMath::Power(0.5*AT_plus.getError(), 2)+pow(0.5*AT_minus.getError(), 2))<<endl;
	myfile<<"fit status: "<<fit_index<<endl;
    myfile.close();

}
