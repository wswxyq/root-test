#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <fstream>
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "TMath.h"
using namespace RooFit ;
using namespace std;

int TMVA_fit_LL_Dplus(){

    //define TChain, TFile 
    const char* TChain_name="MVATree";
    const char* TChain_file="./LL/TMVA_output_LL.root";
    TChain *chain= new TChain(TChain_name);
    chain->Add(TChain_file);
    
	RooRealVar D_M("D_M", "D+ mass", 1799.65, 1939.65, "MeV") ;
    RooRealVar BDT_response("BDT_response", "BDT_response value", 0.901139, RooNumber::infinity());
    RooRealVar P0_ID("P0_ID", "P0_ID value", -1000., 1000.);
	RooArgSet variables(D_M, BDT_response, P0_ID); 

	Double_t D_M_double ;
    Double_t BDT_response_double ;
    chain->SetBranchAddress("D_M",&D_M_double); 
    chain->SetBranchAddress("BDT_response",&BDT_response_double); 

	//creat dataset and write

    RooDataSet *ds=new RooDataSet("ds_plus", "ds_plus", variables, Import(*chain), 
		Cut(" P0_ID>0 "));


	RooRealVar mean("mean","mean of gaussian",1869.,1859.,1879.) ;
	RooRealVar mean_1("mean_1","mean of gaussian(1)",1869.,1859.,1879.) ;
	RooRealVar sigma("sigma","width of gaussian",10.,0.1,20.) ;
	RooRealVar sigma_1("sigma_1","width of gaussian(1)",20.,0.1,30.) ;
	RooGaussian gauss("gauss","gaussian PDF", D_M, mean, sigma) ;  
	RooGaussian gauss_1("gauss_1","gaussian PDF(1)", D_M, mean_1, sigma_1) ;  
	RooRealVar gauss_frac("gauss_frac","fraction of component 1 in signal",0.8,0.,1.) ;
	RooAddPdf gauss_all("gauss_all", "combination of gaussian function", RooArgList(gauss, gauss_1), RooArgList(gauss_frac));


	RooRealVar x1("x1", "para1", -1., 1.);
	RooRealVar x2("x2", "para2", -1., 1.);
	RooRealVar x3("x3", "para3", -1., 1.);
	//RooChebychev poly("poly", "3-para-poly", D_M, RooArgList(x1, x2, x3));
	//RooPolynomial poly("poly", "3-para-poly", D_M, RooArgList(x1));
	RooExponential poly("poly", "3-para-poly", D_M, x1);



	//num total:94333
	RooRealVar gauss_poly_frac("gauss_poly_frac","fraction of component 1 in signal",1000., 600000) ;
	RooRealVar gauss_poly_frac_1("gauss_poly_frac_1","fraction of component 2 in signal",1000., 600000) ;
	RooAddPdf event("even","event",RooArgList(gauss_all, poly), RooArgList(gauss_poly_frac, gauss_poly_frac_1)) ;

    auto result= event.fitTo(*ds, RooFit::NumCPU(10), RooFit::Save(kTRUE), RooFit::Minos(kTRUE));
	Int_t fit_index = result->status();
    
    // Construct plot frame in 'D_M'
	RooPlot* xframe = D_M.frame(Title("Event p.d.f.")) ;
	RooPlot* xframe_2 = D_M.frame(Title(" ")) ;
	xframe_2->SetYTitle("pull distribution");

    // Draw all frames on a canvas
	ds->plotOn(xframe);
	event.plotOn(xframe) ;
	RooHist* hpull = xframe->pullHist() ;
	event.plotOn(xframe,Components(gauss_all),LineColor(kRed),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(poly),LineColor(kBlue),LineStyle(kDashed)) ;
	xframe_2->addPlotable(hpull, "P") ;
	
	TCanvas* c = new TCanvas("total_plot","total_plot", 800, 1200) ;
	c->Divide(1,2) ;
	c->cd(1) ; 
	gPad->SetLeftMargin(0.15) ; 
	xframe->GetYaxis()->SetTitleOffset(1.6) ; 
	xframe->Draw() ;

	c->cd(2) ; 
	gPad->SetLeftMargin(0.15) ; 
	xframe_2->GetYaxis()->SetTitleOffset(1.6) ; 
	xframe_2->GetYaxis()->SetRangeUser(-5., 5.);
	xframe_2->Draw() ;

	c->Print("./LL/fit_Dplus.pdf");

	x1.Print(); x2.Print(); x3.Print(); mean.Print(); sigma.Print(); /*alpha_var.Print(); n_var.Print();*/
	gauss_poly_frac.Print(); gauss_poly_frac_1.Print();
	cout<<"Number of total entries in dataset: "<<ds->numEntries()<<endl;



	std::ofstream myfile;
    myfile.open ("./LL/fit_Dplus.txt");
	myfile<<"number of entries in dataset: "<<ds->numEntries()<<endl;
	myfile<<"x1: "<<x1.getVal()<<endl;
	myfile<<"mass: "<<mean.getVal()<<endl;
	myfile<<"mass_1: "<<mean_1.getVal()<<endl;
	myfile<<"sigma: "<<sigma.getVal()<<endl;
	myfile<<"sigma_1: "<<sigma_1.getVal()<<endl;
	myfile<<"gauss_frac: "<<gauss_frac.getVal()<<endl;
	myfile<<"gauss_poly_frac: "<<gauss_poly_frac.getVal()<<endl;
	myfile<<"gauss_poly_frac_1: "<<gauss_poly_frac_1.getVal()<<endl;
	myfile<<"fit status: "<<fit_index<<endl;
    myfile.close();


    return 0;
}
