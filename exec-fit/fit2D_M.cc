//fit mass spectrum using splot and get s-weight
// generate root file(reduced_DD_1/reduced_DD_2/reduced_LL_1/reduced_LL_2)
// 	with s-weight after mass cut and trigger cut
// also plot the fit result for mass etc.
//follwed by TMVA_general

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
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooHist.h"
#include "TMath.h"

using namespace RooFit ;
using namespace std;
using namespace TMath;

//for DD
int fit_M_D2KSPiPiPi_DD()
{

	//define TChain, TFile 
	const char* TChain_name="D2KSPiPiPi_DD/DecayTree";
	const char* TChain_file="/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root";
	TChain *chain= new TChain(TChain_name);
	chain->Add(TChain_file);
	TChain *chainx= new TChain(TChain_name);
	chainx->Add(TChain_file);
    TFile *newfile = new TFile("../saved/TMVA_root_files/reduced_DD_1.root","recreate");


	//variables
	RooRealVar D_M("D_M", "D+ mass", 1779.65, 1959.65, "MeV") ;
	RooRealVar D_L0HadronDecision_TOS("D_L0HadronDecision_TOS","D_L0HadronDecision_TOS",0,1);
	RooRealVar D_L0MuonDecision_TIS("D_L0MuonDecision_TIS","D_L0MuonDecision_TIS",0,1);
	RooRealVar D_L0DiMuonDecision_TIS("D_L0DiMuonDecision_TIS","D_L0DiMuonDecision_TIS",0,1);
	RooRealVar D_L0ElectronDecision_TIS("D_L0ElectronDecision_TIS","D_L0ElectronDecision_TIS",0,1);
	RooRealVar D_L0PhotonDecision_TIS("D_L0PhotonDecision_TIS","D_L0PhotonDecision_TIS",0,1);
	RooRealVar D_Hlt1TrackMVADecision_TOS("D_Hlt1TrackMVADecision_TOS","D_Hlt1TrackMVADecision_TOS",0,1);
	RooRealVar D_Hlt1TwoTrackMVADecision_TOS("D_Hlt1TwoTrackMVADecision_TOS","D_Hlt1TwoTrackMVADecision_TOS",0,1);
	RooArgSet variables(D_M,D_L0HadronDecision_TOS, D_L0MuonDecision_TIS,D_L0DiMuonDecision_TIS,D_L0ElectronDecision_TIS,
	D_L0PhotonDecision_TIS,D_Hlt1TrackMVADecision_TOS,D_Hlt1TwoTrackMVADecision_TOS); 
	RooRealVar D_PT("D_PT", "D+ transverse momentum", 0., 1000000.) ;
	RooRealVar D_ReFit_P("D_ReFit_P_double", "D_ReFit_P", 0., 100000000.) ;
	RooRealVar D_ReFit_chi2("D_ReFit_chi2_double", "D_ReFit_chi2", 0., 100000.) ;
	RooRealVar KS_PT("KS_PT_double", "KS_PT",  -100000., 100000.) ;
	RooRealVar D_ReFit_KS0_P("D_ReFit_KS0_P_double", "D_ReFit_KS0_P", 0., 100000000.) ;
	RooRealVar P0_ProbNNpi("P0_ProbNNpi", "P0_ProbNNpi", 0., 100.) ;
	RooRealVar P1_ProbNNpi("P1_ProbNNpi", "P1_ProbNNpi", 0., 100.) ;
	RooRealVar P2_ProbNNpi("P2_ProbNNpi", "P2_ProbNNpi", 0., 100.) ;
	RooRealVar D_ReFit_decayLength_vs_Err("D_ReFit_decayLength_vs_Err", "D_ReFit_decayLength_vs_Err", 0, RooNumber::infinity()) ;
	RooRealVar LOG_D_IPCHI2_OWNPV("LOG_D_IPCHI2_OWNPV", "LOG_D_IPCHI2_OWNPV", 0., RooNumber::infinity());
	variables.add(RooArgSet(D_PT, D_ReFit_P, D_ReFit_chi2, KS_PT, D_ReFit_KS0_P));
	variables.add(RooArgSet(P0_ProbNNpi, P1_ProbNNpi, P2_ProbNNpi, D_ReFit_decayLength_vs_Err, LOG_D_IPCHI2_OWNPV));


	//values to write
	Double_t D_M_double ;
	Double_t D_PT_double ;
	Double_t P0_ProbNNpi_double ;
	Double_t P1_ProbNNpi_double ;
	Double_t P2_ProbNNpi_double ;
	Double_t D_IPCHI2_OWNPV_double ;
	Bool_t D_L0HadronDecision_TOS_bool;
	Bool_t D_L0MuonDecision_TIS_bool;
	Bool_t D_L0DiMuonDecision_TIS_bool;
	Bool_t D_L0ElectronDecision_TIS_bool;
	Bool_t D_L0PhotonDecision_TIS_bool;
	Bool_t D_Hlt1TrackMVADecision_TOS_bool;
	Bool_t D_Hlt1TwoTrackMVADecision_TOS_bool;
	Float_t D_ReFit_P_list[1];
	Float_t D_ReFit_chi2_list[1];
	Float_t KS_PT_list[1];
	Float_t D_ReFit_KS0_P_list[1];
	Float_t D_ReFit_decayLength_list[1];
	Float_t D_ReFit_decayLengthErr_list[1];


	//set branch address
    chain->SetBranchAddress("D_M",&D_M_double); 
    chain->SetBranchAddress("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS_bool); 
    chain->SetBranchAddress("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0DiMuonDecision_TIS",&D_L0DiMuonDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0ElectronDecision_TIS",&D_L0ElectronDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0PhotonDecision_TIS",&D_L0PhotonDecision_TIS_bool); 
    chain->SetBranchAddress("D_Hlt1TrackMVADecision_TOS",&D_Hlt1TrackMVADecision_TOS_bool); 
    chain->SetBranchAddress("D_Hlt1TwoTrackMVADecision_TOS",&D_Hlt1TwoTrackMVADecision_TOS_bool); 
    chain->SetBranchAddress("D_PT",&D_PT_double); 
    chain->SetBranchAddress("D_ReFit_P",&D_ReFit_P_list); 
    chain->SetBranchAddress("D_ReFit_chi2",&D_ReFit_chi2_list); 
    chain->SetBranchAddress("KS_PT",&KS_PT_list); 
    chain->SetBranchAddress("D_ReFit_KS0_P",&D_ReFit_KS0_P_list); 
    chain->SetBranchAddress("D_ReFit_decayLength",&D_ReFit_decayLength_list); 
    chain->SetBranchAddress("D_ReFit_decayLengthErr",&D_ReFit_decayLengthErr_list); 
    chain->SetBranchAddress("P0_ProbNNpi",&P0_ProbNNpi_double); 
    chain->SetBranchAddress("P1_ProbNNpi",&P1_ProbNNpi_double); 
    chain->SetBranchAddress("P2_ProbNNpi",&P2_ProbNNpi_double); 
    chain->SetBranchAddress("D_IPCHI2_OWNPV",&D_IPCHI2_OWNPV_double); 


	//creat dataset and write
	RooDataSet *ds=new RooDataSet("ds", "ds", variables);
	cout << "******Start importing....\n";
    Long64_t nentries = chain->GetEntries(); 
    Long64_t nevents=0; //count number of dataset events
	cout << "<<<<Number of Entries in chain:"<<nentries <<endl;
    for (Long64_t i=0;i<nentries;i++) 
    { 
        chain->GetEntry(i);
		
		if (i<100) {/*
			cout << D_M_double<<"\t"<<D_PT_double<<"\t"
			<<P0_ProbNNpi_double<<"\t"<<P1_ProbNNpi_double<<"\t"<<P2_ProbNNpi_double<<"\t"
			<<D_L0HadronDecision_TOS_bool<<"\t"<<D_L0MuonDecision_TIS_bool<<"\t"
			<<D_L0DiMuonDecision_TIS_bool<<"\t"<<D_L0ElectronDecision_TIS_bool<<"\t"
			<<D_L0PhotonDecision_TIS_bool<<"\t"<<D_Hlt1TrackMVADecision_TOS_bool<<"\t"
			<<D_Hlt1TwoTrackMVADecision_TOS_bool<<"\t"<<D_ReFit_P_list[0]<<"\t"
			<<D_ReFit_chi2_list[0]<<"\t"<<KS_PT_list[0]<<"\t"<<D_ReFit_KS0_P_list[0]<<"\t"
			<<D_ReFit_decayLength_list[0]<<"\t"<<D_ReFit_decayLengthErr_list[0]<<"\t"
			<<"\n";		*/
		}
		
		D_M = D_M_double;
		D_L0HadronDecision_TOS = (Double_t)D_L0HadronDecision_TOS_bool;
		D_L0MuonDecision_TIS = (Double_t)D_L0MuonDecision_TIS_bool;
		D_L0DiMuonDecision_TIS = (Double_t)D_L0DiMuonDecision_TIS_bool;
		D_L0ElectronDecision_TIS = (Double_t)D_L0ElectronDecision_TIS_bool;
		D_L0PhotonDecision_TIS = (Double_t)D_L0PhotonDecision_TIS_bool;
		D_Hlt1TrackMVADecision_TOS = (Double_t)D_Hlt1TrackMVADecision_TOS_bool;
		D_Hlt1TwoTrackMVADecision_TOS = (Double_t)D_Hlt1TwoTrackMVADecision_TOS_bool;
        D_PT = D_PT_double; 
        D_ReFit_P = (Double_t)D_ReFit_P_list[0]; 
        D_ReFit_chi2 = (Double_t)D_ReFit_chi2_list[0]; 
        KS_PT = (Double_t)KS_PT_list[0]; 
        D_ReFit_KS0_P = (Double_t)D_ReFit_KS0_P_list[0]; 
		P0_ProbNNpi = P0_ProbNNpi_double;
		P1_ProbNNpi = P1_ProbNNpi_double;
		P2_ProbNNpi = P2_ProbNNpi_double;
        D_ReFit_decayLength_vs_Err = (Double_t)D_ReFit_decayLength_list[0]/(Double_t)D_ReFit_decayLengthErr_list[0]; 
		LOG_D_IPCHI2_OWNPV = TMath::Log(D_IPCHI2_OWNPV_double);

		if (((D_L0HadronDecision_TOS_bool==1) || (D_L0MuonDecision_TIS_bool==1) || (D_L0DiMuonDecision_TIS_bool==1) || 
				(D_L0ElectronDecision_TIS_bool==1) || (D_L0PhotonDecision_TIS_bool==1) ) && ( (D_Hlt1TrackMVADecision_TOS_bool==1) || 
				(D_Hlt1TwoTrackMVADecision_TOS_bool==1)) && (D_M_double>1779.65 && D_M_double<1959.65)&&(D_IPCHI2_OWNPV_double>0) ) {
			ds->add(variables,1.0);
			nevents++;
		}
     } 
		//chain->Print();
		cout << "<<<<number of events on dataset:" << nevents <<endl;
		cout << "******Importing complete. RooDataSet ds created!\n";


	//creat newtree	
	chainx->SetBranchStatus("*",0);
    chainx->SetBranchStatus("D_M",1);
    chainx->SetBranchStatus("D_L0HadronDecision_TOS",1);
    chainx->SetBranchStatus("D_L0MuonDecision_TIS",1);
    chainx->SetBranchStatus("D_L0DiMuonDecision_TIS",1);
    chainx->SetBranchStatus("D_L0ElectronDecision_TIS",1);
    chainx->SetBranchStatus("D_L0PhotonDecision_TIS",1);
    chainx->SetBranchStatus("D_Hlt1TrackMVADecision_TOS",1);
    chainx->SetBranchStatus("D_Hlt1TwoTrackMVADecision_TOS",1);
    chainx->SetBranchStatus("D_PT",1);
    chainx->SetBranchStatus("D_ReFit_P",1);
    chainx->SetBranchStatus("D_ReFit_chi2",1);
    chainx->SetBranchStatus("KS_PT",1);
    chainx->SetBranchStatus("D_ReFit_KS0_P",1);
    chainx->SetBranchStatus("D_ReFit_decayLength",1);
    chainx->SetBranchStatus("D_ReFit_decayLengthErr",1);
    chainx->SetBranchStatus("P0_ProbNNpi",1);
    chainx->SetBranchStatus("P1_ProbNNpi",1);
    chainx->SetBranchStatus("P2_ProbNNpi",1);
    chainx->SetBranchStatus("D_IPCHI2_OWNPV",1);
	chainx->SetBranchStatus("P0_ID",1);
	chainx->SetBranchStatus("P0_P",1);
	chainx->SetBranchStatus("P0_PX",1);
	chainx->SetBranchStatus("P0_PY",1);
	chainx->SetBranchStatus("P0_PZ",1);
	chainx->SetBranchStatus("P1_P",1);
	chainx->SetBranchStatus("P1_PX",1);
	chainx->SetBranchStatus("P1_PY",1);
	chainx->SetBranchStatus("P1_PZ",1);
	chainx->SetBranchStatus("P2_P",1);
	chainx->SetBranchStatus("P2_PX",1);
	chainx->SetBranchStatus("P2_PY",1);
	chainx->SetBranchStatus("P2_PZ",1);

	TTree *newtree = (TTree*)chainx->CopyTree("(   ( (D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
        (D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1) ) && ( (D_Hlt1TrackMVADecision_TOS==1) || \
        (D_Hlt1TwoTrackMVADecision_TOS==1) )   )&& (D_M>1779.65&&D_M<1959.65) &&(D_IPCHI2_OWNPV>0)");
    //newtree->Print();


	//check if nentries_newtree & nevents matches
	Long64_t nentries_newtree=newtree->GetEntries();
	if(nentries_newtree!=nevents) { 
		cout << "nentries_newtree!=nevents" <<endl;
		return 0; 
		}


	RooRealVar mean("mean","mean of gaussian",1869.,1859.,1879.) ;
	RooRealVar mean_1("mean_1","mean of gaussian(1)",1869.,1859.,1879.) ;
	RooRealVar sigma("sigma","width of gaussian",10.,0.1,20.) ;
	RooRealVar sigma_1("sigma_1","width of gaussian(1)",20.,0.1,30.) ;
	RooGaussian gauss("gauss","gaussian PDF", D_M, mean, sigma) ;  
	RooGaussian gauss_1("gauss_1","gaussian PDF(1)", D_M, mean_1, sigma_1) ;  
	RooRealVar gauss_frac("gauss_frac","fraction of component 1 in signal",0.8,0.,1.) ;
	RooAddPdf gauss_all("gauss_all", "combination of gaussian function", RooArgList(gauss, gauss_1), RooArgList(gauss_frac));

	RooRealVar x1("x1", "para1", 0., 100.);
	RooRealVar x2("x2", "para2", 0., 100.);
	RooRealVar x3("x3", "para3", 0., 100.);
	RooPolynomial poly("poly", "3-para-poly", D_M, RooArgList(x1, x2, x3),0);

	RooRealVar gauss_poly_frac("gauss_poly_frac","fraction of signal",20000,300000) ;
	RooRealVar gauss_poly_frac_1("gauss_poly_frac_1","fraction of background",20000,300000) ;
	RooAddPdf event("even","event",RooArgList(gauss_all, poly), RooArgList(gauss_poly_frac, gauss_poly_frac_1)) ;


	// Construct plot frame in 'D_M'
	RooPlot* xframe = D_M.frame(Title("Event p.d.f.")) ;
	RooPlot* xframe_2 = D_M.frame(Title(" ")) ;
	xframe_2->SetYTitle("pull distribution");

	event.fitTo(*ds) ;
	cout<<"fitting complete!\n";
	//print variable values
	cout<<"para list:\n";
	mean.Print() ;  mean_1.Print() ;  sigma.Print() ;  sigma_1.Print() ;  x1.Print();  x2.Print();  x3.Print();
	gauss_frac.Print(); 
	gauss_poly_frac.Print();  gauss_poly_frac_1.Print(); 

	// Draw all frames on a canvas
	// Draw part 1
	ds->plotOn(xframe);
	event.plotOn(xframe) ;
	RooHist* hpull = xframe->pullHist() ;
	event.plotOn(xframe,Components(poly),LineColor(kRed),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(gauss),LineColor(kBlue),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(gauss_1),LineColor(kGreen),LineStyle(kDashed)) ;
	xframe_2->addPlotable(hpull, "P") ;
	
	TCanvas* c = new TCanvas("total_plot","total_plot",600,1200) ;
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

	c->Print("../plotdir/D2KSPiPiPi_DD/fit_mass_D2KSPiPiPi_DD.pdf");



	// Draw part 2
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *ds, &event, RooArgList(gauss_poly_frac, gauss_poly_frac_1) );

	std::cout << "Check SWeights:" << std::endl;
	std::cout << std::endl <<  "Yield of gauss_poly_frac is "
	<< gauss_poly_frac.getVal() << ".  From sWeights it is "
	<< sData->GetYieldFromSWeight("gauss_poly_frac") << std::endl;
	std::cout << std::endl <<  "Yield of gauss_poly_frac_1 is "
	<< gauss_poly_frac_1.getVal() << ".  From sWeights it is "
	<< sData->GetYieldFromSWeight("gauss_poly_frac_1") << std::endl;


	
    Double_t swei, bwei;
	swei= 1. ;
    bwei= 1. ;
    TBranch *swei_branch = newtree->Branch("swei", &swei, "swei/D");
    TBranch *bwei_branch = newtree->Branch("bwei", &bwei, "bwei/D");
	
	for(Int_t sdata_index = 0; sdata_index < nentries_newtree; sdata_index++)
	{
		newtree->GetEntry(sdata_index);
		swei=sData->GetSWeight(sdata_index, "gauss_poly_frac");
		bwei=sData->GetSWeight(sdata_index, "gauss_poly_frac_1");
		//cout<<"this is the "<<sdata_index<<"th item of weights"<<endl;
		//cout<<swei<<endl;
		//cout<<bwei<<endl;
		swei_branch->Fill();
		bwei_branch->Fill();

	}
	//newtree->Print();
	newfile->cd();
    newtree->Write("ReducedTree");

	// create weightfed data set
	RooDataSet * signalpart = new RooDataSet(ds->GetName(),ds->GetTitle(),ds,*ds->get(),0,"gauss_poly_frac_sw") ;
	RooDataSet * backgroundpart = new RooDataSet(ds->GetName(),ds->GetTitle(),ds,*ds->get(),0,"gauss_poly_frac_1_sw") ;

	//D_PT
	TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	RooPlot *frame = D_PT.frame(Range(0.,20000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
	RooPlot* frame1 = D_PT.frame(Range(0.,20000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_D_PT_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_P
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_P.frame(Range(0.,300000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_P.frame(Range(0.,300000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_D_ReFit_P_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_chi2
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_chi2.frame(Range(0.,100.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_chi2.frame(Range(0.,100.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_D_ReFit_chi2_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//KS_PT
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = KS_PT.frame(Range(0.,10000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = KS_PT.frame(Range(0.,10000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_KS_PT_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_KS0_P
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_KS0_P.frame(Range(0.,100000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_KS0_P.frame(Range(0.,100000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_D_ReFit_KS0_P_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//P1_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P1_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P1_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_P1_ProbNNpi_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//P2_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P2_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P2_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_P2_ProbNNpi_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//P0_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P0_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P0_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_P0_ProbNNpi_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_decayLength_vs_Err
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_decayLength_vs_Err.frame(Range(0.,300.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_decayLength_vs_Err.frame(Range(0.,300.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_D_ReFit_decayLength_vs_Err_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//LOG_D_IPCHI2_OWNPV
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = LOG_D_IPCHI2_OWNPV.frame(Range(-2.,8.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = LOG_D_IPCHI2_OWNPV.frame(Range(-2.,8.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_DD/SPlot/splot_fit_LOG_D_IPCHI2_OWNPV_D2KSPiPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//release pointers
	delete xframe;
	delete xframe_2;
	delete chain;
	delete c;
	return 0;
}



int fit_M_D2KSKpPiPi_DD()
{

	//define TChain, TFile 
	const char* TChain_name="D2KSKpPiPi_DD/DecayTree";
	const char* TChain_file="/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root";
	TChain *chain= new TChain(TChain_name);
	chain->Add(TChain_file);
	TChain *chainx= new TChain(TChain_name);
	chainx->Add(TChain_file);
    TFile *newfile = new TFile("../saved/TMVA_root_files/reduced_DD_2.root","recreate");


	//variables
	RooRealVar D_M("D_M", "D+ mass", 1779.65, 1959.65, "MeV") ;
	RooRealVar D_L0HadronDecision_TOS("D_L0HadronDecision_TOS","D_L0HadronDecision_TOS",0,1);
	RooRealVar D_L0MuonDecision_TIS("D_L0MuonDecision_TIS","D_L0MuonDecision_TIS",0,1);
	RooRealVar D_L0DiMuonDecision_TIS("D_L0DiMuonDecision_TIS","D_L0DiMuonDecision_TIS",0,1);
	RooRealVar D_L0ElectronDecision_TIS("D_L0ElectronDecision_TIS","D_L0ElectronDecision_TIS",0,1);
	RooRealVar D_L0PhotonDecision_TIS("D_L0PhotonDecision_TIS","D_L0PhotonDecision_TIS",0,1);
	RooRealVar D_Hlt1TrackMVADecision_TOS("D_Hlt1TrackMVADecision_TOS","D_Hlt1TrackMVADecision_TOS",0,1);
	RooRealVar D_Hlt1TwoTrackMVADecision_TOS("D_Hlt1TwoTrackMVADecision_TOS","D_Hlt1TwoTrackMVADecision_TOS",0,1);
	RooArgSet variables(D_M,D_L0HadronDecision_TOS, D_L0MuonDecision_TIS,D_L0DiMuonDecision_TIS,D_L0ElectronDecision_TIS,
	D_L0PhotonDecision_TIS,D_Hlt1TrackMVADecision_TOS,D_Hlt1TwoTrackMVADecision_TOS); 
	RooRealVar D_PT("D_PT", "D+ transverse momentum", 0., 1000000.) ;
	RooRealVar D_ReFit_P("D_ReFit_P_double", "D_ReFit_P", 0., 100000000.) ;
	RooRealVar D_ReFit_chi2("D_ReFit_chi2_double", "D_ReFit_chi2", 0., 100000.) ;
	RooRealVar KS_PT("KS_PT_double", "KS_PT",  -100000., 100000.) ;
	RooRealVar D_ReFit_KS0_P("D_ReFit_KS0_P_double", "D_ReFit_KS0_P", 0., 100000000.) ;
	RooRealVar P0_ProbNNpi("P0_ProbNNpi", "P0_ProbNNpi", 0., 100.) ;
	RooRealVar P1_ProbNNpi("P1_ProbNNpi", "P1_ProbNNpi", 0., 100.) ;
	RooRealVar P2_ProbNNpi("P2_ProbNNpi", "P2_ProbNNpi", 0., 100.) ;
	RooRealVar D_ReFit_decayLength_vs_Err("D_ReFit_decayLength_vs_Err", "D_ReFit_decayLength_vs_Err", 0, RooNumber::infinity()) ;
	RooRealVar LOG_D_IPCHI2_OWNPV("LOG_D_IPCHI2_OWNPV", "LOG_D_IPCHI2_OWNPV", 0., RooNumber::infinity());
	variables.add(RooArgSet(D_PT, D_ReFit_P, D_ReFit_chi2, KS_PT, D_ReFit_KS0_P));
	variables.add(RooArgSet(P0_ProbNNpi, P1_ProbNNpi, P2_ProbNNpi, D_ReFit_decayLength_vs_Err, LOG_D_IPCHI2_OWNPV));


	//values to write
	Double_t D_M_double ;
	Double_t D_PT_double ;
	Double_t P0_ProbNNpi_double ;
	Double_t P1_ProbNNpi_double ;
	Double_t P2_ProbNNpi_double ;
	Double_t D_IPCHI2_OWNPV_double ;
	Bool_t D_L0HadronDecision_TOS_bool;
	Bool_t D_L0MuonDecision_TIS_bool;
	Bool_t D_L0DiMuonDecision_TIS_bool;
	Bool_t D_L0ElectronDecision_TIS_bool;
	Bool_t D_L0PhotonDecision_TIS_bool;
	Bool_t D_Hlt1TrackMVADecision_TOS_bool;
	Bool_t D_Hlt1TwoTrackMVADecision_TOS_bool;
	Float_t D_ReFit_P_list[1];
	Float_t D_ReFit_chi2_list[1];
	Float_t KS_PT_list[1];
	Float_t D_ReFit_KS0_P_list[1];
	Float_t D_ReFit_decayLength_list[1];
	Float_t D_ReFit_decayLengthErr_list[1];


	//set branch address
    chain->SetBranchAddress("D_M",&D_M_double); 
    chain->SetBranchAddress("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS_bool); 
    chain->SetBranchAddress("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0DiMuonDecision_TIS",&D_L0DiMuonDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0ElectronDecision_TIS",&D_L0ElectronDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0PhotonDecision_TIS",&D_L0PhotonDecision_TIS_bool); 
    chain->SetBranchAddress("D_Hlt1TrackMVADecision_TOS",&D_Hlt1TrackMVADecision_TOS_bool); 
    chain->SetBranchAddress("D_Hlt1TwoTrackMVADecision_TOS",&D_Hlt1TwoTrackMVADecision_TOS_bool); 
    chain->SetBranchAddress("D_PT",&D_PT_double); 
    chain->SetBranchAddress("D_ReFit_P",&D_ReFit_P_list); 
    chain->SetBranchAddress("D_ReFit_chi2",&D_ReFit_chi2_list); 
    chain->SetBranchAddress("KS_PT",&KS_PT_list); 
    chain->SetBranchAddress("D_ReFit_KS0_P",&D_ReFit_KS0_P_list); 
    chain->SetBranchAddress("D_ReFit_decayLength",&D_ReFit_decayLength_list); 
    chain->SetBranchAddress("D_ReFit_decayLengthErr",&D_ReFit_decayLengthErr_list); 
    chain->SetBranchAddress("P0_ProbNNpi",&P0_ProbNNpi_double); 
    chain->SetBranchAddress("P1_ProbNNpi",&P1_ProbNNpi_double); 
    chain->SetBranchAddress("P2_ProbNNpi",&P2_ProbNNpi_double); 
    chain->SetBranchAddress("D_IPCHI2_OWNPV",&D_IPCHI2_OWNPV_double); 
	


	//creat dataset and write
	RooDataSet *ds=new RooDataSet("ds", "ds", variables);
	cout << "******Start importing....\n";
    Long64_t nentries = chain->GetEntries(); 
    Long64_t nevents=0; //count number of dataset events
	cout << "<<<<Number of Entries in chain:"<<nentries <<endl;
    for (Long64_t i=0;i<nentries;i++) 
    { 
        chain->GetEntry(i);
		
		if (i<100) {/*
			cout << D_M_double<<"\t"<<D_PT_double<<"\t"
			<<P0_ProbNNpi_double<<"\t"<<P1_ProbNNpi_double<<"\t"<<P2_ProbNNpi_double<<"\t"
			<<D_L0HadronDecision_TOS_bool<<"\t"<<D_L0MuonDecision_TIS_bool<<"\t"
			<<D_L0DiMuonDecision_TIS_bool<<"\t"<<D_L0ElectronDecision_TIS_bool<<"\t"
			<<D_L0PhotonDecision_TIS_bool<<"\t"<<D_Hlt1TrackMVADecision_TOS_bool<<"\t"
			<<D_Hlt1TwoTrackMVADecision_TOS_bool<<"\t"<<D_ReFit_P_list[0]<<"\t"
			<<D_ReFit_chi2_list[0]<<"\t"<<KS_PT_list[0]<<"\t"<<D_ReFit_KS0_P_list[0]<<"\t"
			<<D_ReFit_decayLength_list[0]<<"\t"<<D_ReFit_decayLengthErr_list[0]<<"\t"
			<<"\n";		*/
		}
		
		D_M = D_M_double;
		D_L0HadronDecision_TOS = (Double_t)D_L0HadronDecision_TOS_bool;
		D_L0MuonDecision_TIS = (Double_t)D_L0MuonDecision_TIS_bool;
		D_L0DiMuonDecision_TIS = (Double_t)D_L0DiMuonDecision_TIS_bool;
		D_L0ElectronDecision_TIS = (Double_t)D_L0ElectronDecision_TIS_bool;
		D_L0PhotonDecision_TIS = (Double_t)D_L0PhotonDecision_TIS_bool;
		D_Hlt1TrackMVADecision_TOS = (Double_t)D_Hlt1TrackMVADecision_TOS_bool;
		D_Hlt1TwoTrackMVADecision_TOS = (Double_t)D_Hlt1TwoTrackMVADecision_TOS_bool;
        D_PT = D_PT_double; 
        D_ReFit_P = (Double_t)D_ReFit_P_list[0]; 
        D_ReFit_chi2 = (Double_t)D_ReFit_chi2_list[0]; 
        KS_PT = (Double_t)KS_PT_list[0]; 
        D_ReFit_KS0_P = (Double_t)D_ReFit_KS0_P_list[0]; 
		P0_ProbNNpi = P0_ProbNNpi_double;
		P1_ProbNNpi = P1_ProbNNpi_double;
		P2_ProbNNpi = P2_ProbNNpi_double;
        D_ReFit_decayLength_vs_Err = (Double_t)D_ReFit_decayLength_list[0]/(Double_t)D_ReFit_decayLengthErr_list[0]; 
		LOG_D_IPCHI2_OWNPV = TMath::Log(D_IPCHI2_OWNPV_double);

		if (((D_L0HadronDecision_TOS_bool==1) || (D_L0MuonDecision_TIS_bool==1) || (D_L0DiMuonDecision_TIS_bool==1) || 
				(D_L0ElectronDecision_TIS_bool==1) || (D_L0PhotonDecision_TIS_bool==1) ) && ( (D_Hlt1TrackMVADecision_TOS_bool==1) || 
				(D_Hlt1TwoTrackMVADecision_TOS_bool==1)) && (D_M_double>1779.65 && D_M_double<1959.65)&&(D_IPCHI2_OWNPV_double>0) ) {
			ds->add(variables,1.0);
			nevents++;
		}
     } 
		//chain->Print();
		cout << "<<<<number of events on dataset:" << nevents <<endl;
		cout << "******Importing complete. RooDataSet ds created!\n";


	//creat newtree	
	chainx->SetBranchStatus("*",0);
    chainx->SetBranchStatus("D_M",1);
    chainx->SetBranchStatus("D_L0HadronDecision_TOS",1);
    chainx->SetBranchStatus("D_L0MuonDecision_TIS",1);
    chainx->SetBranchStatus("D_L0DiMuonDecision_TIS",1);
    chainx->SetBranchStatus("D_L0ElectronDecision_TIS",1);
    chainx->SetBranchStatus("D_L0PhotonDecision_TIS",1);
    chainx->SetBranchStatus("D_Hlt1TrackMVADecision_TOS",1);
    chainx->SetBranchStatus("D_Hlt1TwoTrackMVADecision_TOS",1);
    chainx->SetBranchStatus("D_PT",1);
    chainx->SetBranchStatus("D_ReFit_P",1);
    chainx->SetBranchStatus("D_ReFit_chi2",1);
    chainx->SetBranchStatus("KS_PT",1);
    chainx->SetBranchStatus("D_ReFit_KS0_P",1);
    chainx->SetBranchStatus("D_ReFit_decayLength",1);
    chainx->SetBranchStatus("D_ReFit_decayLengthErr",1);
    chainx->SetBranchStatus("P0_ProbNNpi",1);
    chainx->SetBranchStatus("P1_ProbNNpi",1);
    chainx->SetBranchStatus("P2_ProbNNpi",1);
    chainx->SetBranchStatus("D_IPCHI2_OWNPV",1);
	chainx->SetBranchStatus("P0_ID",1);
	chainx->SetBranchStatus("P0_P",1);
	chainx->SetBranchStatus("P0_PX",1);
	chainx->SetBranchStatus("P0_PY",1);
	chainx->SetBranchStatus("P0_PZ",1);
	chainx->SetBranchStatus("P1_P",1);
	chainx->SetBranchStatus("P1_PX",1);
	chainx->SetBranchStatus("P1_PY",1);
	chainx->SetBranchStatus("P1_PZ",1);
	chainx->SetBranchStatus("P2_P",1);
	chainx->SetBranchStatus("P2_PX",1);
	chainx->SetBranchStatus("P2_PY",1);
	chainx->SetBranchStatus("P2_PZ",1);
	TTree *newtree = (TTree*)chainx->CopyTree("(   ( (D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
        (D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1) ) && ( (D_Hlt1TrackMVADecision_TOS==1) || \
        (D_Hlt1TwoTrackMVADecision_TOS==1) )   )&& (D_M>1779.65&&D_M<1959.65) &&(D_IPCHI2_OWNPV>0)");
    //newtree->Print();


	//check if nentries_newtree & nevents matches
	Long64_t nentries_newtree=newtree->GetEntries();
	if(nentries_newtree!=nevents) { 
		cout << "nentries_newtree!=nevents" <<endl;
		return 0; 
		}


	RooRealVar mean("mean","mean of gaussian",1869.,1859.,1879.) ;
	RooRealVar sigma("sigma","width of gaussian",10.,0.1,100.) ;
	RooRealVar alpha_var("alpha_var","alpha_var",1,0.1,10.) ;
	RooRealVar n_var("n_var","n_var",30.7631, 0. , 100.) ;
	RooCBShape gauss_all("CBfunction", "CBfunction", D_M, mean, sigma, alpha_var, n_var);

	RooRealVar x1("x1", "para1", -100., 100.);
	RooRealVar x2("x2", "para2", -100., 100.);
	RooRealVar x3("x3", "para3", -100., 100.);
	RooChebychev poly("poly", "3-para-poly", D_M, RooArgList(x1, x2, x3));
	//num total:94333
	RooRealVar gauss_poly_frac("gauss_poly_frac","fraction of component 1 in signal",1000., 10000) ;
	RooRealVar gauss_poly_frac_1("gauss_poly_frac_1","fraction of component 2 in signal",10000., 100000) ;
	RooAddPdf event("even","event",RooArgList(gauss_all, poly), RooArgList(gauss_poly_frac, gauss_poly_frac_1)) ;

	// Construct plot frame in 'D_M'
	RooPlot* xframe = D_M.frame(Title("Event p.d.f.")) ;
	RooPlot* xframe_2 = D_M.frame(Title(" ")) ;
	xframe_2->SetYTitle("pull distribution");

	event.fitTo(*ds) ;
	cout<<"fitting complete!\n";
	//print variable values
	cout<<"para list:\n";
	mean.Print() ;  sigma.Print() ;  alpha_var.Print() ;  n_var.Print() ;  x1.Print();  x2.Print();  x3.Print();
	gauss_poly_frac.Print(); 
	gauss_poly_frac_1.Print(); 

	// Draw all frames on a canvas
	// Draw part 1
	ds->plotOn(xframe);
	event.plotOn(xframe) ;
	RooHist* hpull = xframe->pullHist() ;
	event.plotOn(xframe,Components(poly),LineColor(kRed),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(gauss_all),LineColor(kBlue),LineStyle(kDashed)) ;
	xframe_2->addPlotable(hpull, "P") ;
	
	TCanvas* c = new TCanvas("total_plot","total_plot",600,1200) ;
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

	c->Print("../plotdir/D2KSKpPiPi_DD/fit_mass_D2KSKpPiPi_DD.pdf");



	// Draw part 2
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *ds, &event, RooArgList(gauss_poly_frac, gauss_poly_frac_1) );

	std::cout << "Check SWeights:" << std::endl;
	std::cout << std::endl <<  "Yield of gauss_poly_frac is "
	<< gauss_poly_frac.getVal() << ".  From sWeights it is "
	<< sData->GetYieldFromSWeight("gauss_poly_frac") << std::endl;
	std::cout << std::endl <<  "Yield of gauss_poly_frac_1 is "
	<< gauss_poly_frac_1.getVal() << ".  From sWeights it is "
	<< sData->GetYieldFromSWeight("gauss_poly_frac_1") << std::endl;


	
    Double_t swei, bwei;
	swei= 1. ;
    bwei= 1. ;
    TBranch *swei_branch = newtree->Branch("swei", &swei, "swei/D");
    TBranch *bwei_branch = newtree->Branch("bwei", &bwei, "bwei/D");
	
	for(Int_t sdata_index = 0; sdata_index < nentries_newtree; sdata_index++)
	{
		newtree->GetEntry(sdata_index);
		swei=sData->GetSWeight(sdata_index, "gauss_poly_frac");
		bwei=sData->GetSWeight(sdata_index, "gauss_poly_frac_1");
		//cout<<"this is the "<<sdata_index<<"th item of weights"<<endl;
		//cout<<swei<<endl;
		//cout<<bwei<<endl;
		swei_branch->Fill();
		bwei_branch->Fill();

	}
	//newtree->Print();
	newfile->cd();
    newtree->Write("ReducedTree");

	// create weightfed data set
	RooDataSet * signalpart = new RooDataSet(ds->GetName(),ds->GetTitle(),ds,*ds->get(),0,"gauss_poly_frac_sw") ;
	RooDataSet * backgroundpart = new RooDataSet(ds->GetName(),ds->GetTitle(),ds,*ds->get(),0,"gauss_poly_frac_1_sw") ;

	//D_PT
	TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	RooPlot *frame = D_PT.frame(Range(0.,20000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
	RooPlot* frame1 = D_PT.frame(Range(0.,20000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_D_PT_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_P
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_P.frame(Range(0.,300000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_P.frame(Range(0.,300000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_D_ReFit_P_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_chi2
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_chi2.frame(Range(0.,100.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_chi2.frame(Range(0.,100.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_D_ReFit_chi2_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//KS_PT
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = KS_PT.frame(Range(0.,10000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = KS_PT.frame(Range(0.,10000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_KS_PT_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_KS0_P
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_KS0_P.frame(Range(0.,100000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_KS0_P.frame(Range(0.,100000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_D_ReFit_KS0_P_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//P1_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P1_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P1_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_P1_ProbNNpi_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//P2_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P2_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P2_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_P2_ProbNNpi_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//P0_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P0_ProbNNpi.frame(Range(0.,0.2), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P0_ProbNNpi.frame(Range(0.,0.2), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_P0_ProbNNpi_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_decayLength_vs_Err
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_decayLength_vs_Err.frame(Range(0.,300.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_decayLength_vs_Err.frame(Range(0.,300.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_D_ReFit_decayLength_vs_Err_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//LOG_D_IPCHI2_OWNPV
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = LOG_D_IPCHI2_OWNPV.frame(Range(-2.,8.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = LOG_D_IPCHI2_OWNPV.frame(Range(-2.,8.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_DD/SPlot/splot_fit_LOG_D_IPCHI2_OWNPV_D2KSKpPiPi_DD.pdf");
	delete cdata; delete frame; delete frame1;

	//release pointers
	delete xframe;
	delete xframe_2;
	delete chain;
	delete c;
	return 0;
}







//for LL
int fit_M_D2KSPiPiPi_LL()
{
	
	cout << "<<<< define TChain, TFile ..."<<endl;

	const char* TChain_name="D2KSPiPiPi_LL/DecayTree";
	const char* TChain_file="/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root";
	TChain *chain= new TChain(TChain_name);
	chain->Add(TChain_file);
	TChain *chainx= new TChain(TChain_name);
	chainx->Add(TChain_file);
    TFile *newfile = new TFile("../saved/TMVA_root_files/reduced_LL_1.root","recreate");


	cout << "<<<< define variables..."<<endl;

	RooRealVar D_M("D_M", "D+ mass", 1779.65, 1959.65, "MeV") ;
	RooRealVar D_L0HadronDecision_TOS("D_L0HadronDecision_TOS","D_L0HadronDecision_TOS",0,1);
	RooRealVar D_L0MuonDecision_TIS("D_L0MuonDecision_TIS","D_L0MuonDecision_TIS",0,1);
	RooRealVar D_L0DiMuonDecision_TIS("D_L0DiMuonDecision_TIS","D_L0DiMuonDecision_TIS",0,1);
	RooRealVar D_L0ElectronDecision_TIS("D_L0ElectronDecision_TIS","D_L0ElectronDecision_TIS",0,1);
	RooRealVar D_L0PhotonDecision_TIS("D_L0PhotonDecision_TIS","D_L0PhotonDecision_TIS",0,1);
	RooRealVar D_Hlt1TrackMVADecision_TOS("D_Hlt1TrackMVADecision_TOS","D_Hlt1TrackMVADecision_TOS",0,1);
	RooRealVar D_Hlt1TwoTrackMVADecision_TOS("D_Hlt1TwoTrackMVADecision_TOS","D_Hlt1TwoTrackMVADecision_TOS",0,1);
	RooArgSet variables(D_M,D_L0HadronDecision_TOS, D_L0MuonDecision_TIS,D_L0DiMuonDecision_TIS,D_L0ElectronDecision_TIS,
	D_L0PhotonDecision_TIS,D_Hlt1TrackMVADecision_TOS,D_Hlt1TwoTrackMVADecision_TOS); 
	RooRealVar D_PT("D_PT", "D+ transverse momentum", 0., 1000000.) ;
	RooRealVar D_ReFit_P("D_ReFit_P_double", "D_ReFit_P", 0., 100000000.) ;
	RooRealVar D_ReFit_chi2("D_ReFit_chi2_double", "D_ReFit_chi2", 0., 100000.) ;
	RooRealVar KS_PT("KS_PT_double", "KS_PT",  -100000., 100000.) ;
	RooRealVar D_ReFit_KS0_P("D_ReFit_KS0_P_double", "D_ReFit_KS0_P", 0., 100000000.) ;
	RooRealVar P0_ProbNNpi("P0_ProbNNpi", "P0_ProbNNpi", 0., 100.) ;
	RooRealVar P1_ProbNNpi("P1_ProbNNpi", "P1_ProbNNpi", 0., 100.) ;
	RooRealVar P2_ProbNNpi("P2_ProbNNpi", "P2_ProbNNpi", 0., 100.) ;
	RooRealVar D_ReFit_decayLength_vs_Err("D_ReFit_decayLength_vs_Err", "D_ReFit_decayLength_vs_Err", 0, RooNumber::infinity()) ;
	RooRealVar LOG_D_IPCHI2_OWNPV("LOG_D_IPCHI2_OWNPV", "LOG_D_IPCHI2_OWNPV", 0., RooNumber::infinity());
	variables.add(RooArgSet(D_PT, D_ReFit_P, D_ReFit_chi2, KS_PT, D_ReFit_KS0_P));
	variables.add(RooArgSet(P0_ProbNNpi, P1_ProbNNpi, P2_ProbNNpi, D_ReFit_decayLength_vs_Err, LOG_D_IPCHI2_OWNPV));


	cout << "<<<< define values to write..."<<endl;

	Double_t D_M_double ;
	Double_t D_PT_double ;
	Double_t P0_ProbNNpi_double ;
	Double_t P1_ProbNNpi_double ;
	Double_t P2_ProbNNpi_double ;
	Double_t D_IPCHI2_OWNPV_double ;
	Bool_t D_L0HadronDecision_TOS_bool;
	Bool_t D_L0MuonDecision_TIS_bool;
	Bool_t D_L0DiMuonDecision_TIS_bool;
	Bool_t D_L0ElectronDecision_TIS_bool;
	Bool_t D_L0PhotonDecision_TIS_bool;
	Bool_t D_Hlt1TrackMVADecision_TOS_bool;
	Bool_t D_Hlt1TwoTrackMVADecision_TOS_bool;
	Float_t D_ReFit_P_list[1];
	Float_t D_ReFit_chi2_list[1];
	Float_t KS_PT_list[1];
	Float_t D_ReFit_KS0_P_list[1];
	Float_t D_ReFit_decayLength_list[1];
	Float_t D_ReFit_decayLengthErr_list[1];


	cout << "<<<< set branch address..."<<endl;
	
    chain->SetBranchAddress("D_M",&D_M_double); 
    chain->SetBranchAddress("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS_bool); 
    chain->SetBranchAddress("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0DiMuonDecision_TIS",&D_L0DiMuonDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0ElectronDecision_TIS",&D_L0ElectronDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0PhotonDecision_TIS",&D_L0PhotonDecision_TIS_bool); 
    chain->SetBranchAddress("D_Hlt1TrackMVADecision_TOS",&D_Hlt1TrackMVADecision_TOS_bool); 
    chain->SetBranchAddress("D_Hlt1TwoTrackMVADecision_TOS",&D_Hlt1TwoTrackMVADecision_TOS_bool); 
    chain->SetBranchAddress("D_PT",&D_PT_double); 
    chain->SetBranchAddress("D_ReFit_P",&D_ReFit_P_list); 
    chain->SetBranchAddress("D_ReFit_chi2",&D_ReFit_chi2_list); 
    chain->SetBranchAddress("KS_PT",&KS_PT_list); 
    chain->SetBranchAddress("D_ReFit_KS0_P",&D_ReFit_KS0_P_list); 
    chain->SetBranchAddress("D_ReFit_decayLength",&D_ReFit_decayLength_list); 
    chain->SetBranchAddress("D_ReFit_decayLengthErr",&D_ReFit_decayLengthErr_list); 
    chain->SetBranchAddress("P0_ProbNNpi",&P0_ProbNNpi_double); 
    chain->SetBranchAddress("P1_ProbNNpi",&P1_ProbNNpi_double); 
    chain->SetBranchAddress("P2_ProbNNpi",&P2_ProbNNpi_double); 
    chain->SetBranchAddress("D_IPCHI2_OWNPV",&D_IPCHI2_OWNPV_double); 


	cout << "<<<< creat dataset and write..."<<endl;
	
	RooDataSet *ds=new RooDataSet("ds", "ds", variables);
    Long64_t nentries = chain->GetEntries(); 
    Long64_t nevents=0; //count number of dataset events
	cout << "<<<< Number of Entries in chain(original):"<<nentries <<endl;
	cout << "<<<< Start importing data to *ds....\n";

    for (Long64_t i=0;i<nentries;i++) 
    { 
        chain->GetEntry(i);
		/*
		if (i<100) {
			cout << D_M_double<<"\t"<<D_PT_double<<"\t"
			<<P0_ProbNNpi_double<<"\t"<<P1_ProbNNpi_double<<"\t"<<P2_ProbNNpi_double<<"\t"
			<<D_L0HadronDecision_TOS_bool<<"\t"<<D_L0MuonDecision_TIS_bool<<"\t"
			<<D_L0DiMuonDecision_TIS_bool<<"\t"<<D_L0ElectronDecision_TIS_bool<<"\t"
			<<D_L0PhotonDecision_TIS_bool<<"\t"<<D_Hlt1TrackMVADecision_TOS_bool<<"\t"
			<<D_Hlt1TwoTrackMVADecision_TOS_bool<<"\t"<<D_ReFit_P_list[0]<<"\t"
			<<D_ReFit_chi2_list[0]<<"\t"<<KS_PT_list[0]<<"\t"<<D_ReFit_KS0_P_list[0]<<"\t"
			<<D_ReFit_decayLength_list[0]<<"\t"<<D_ReFit_decayLengthErr_list[0]<<"\t"
			<<"\n";		
		}*/
		
		D_M = D_M_double;
		D_L0HadronDecision_TOS = (Double_t)D_L0HadronDecision_TOS_bool;
		D_L0MuonDecision_TIS = (Double_t)D_L0MuonDecision_TIS_bool;
		D_L0DiMuonDecision_TIS = (Double_t)D_L0DiMuonDecision_TIS_bool;
		D_L0ElectronDecision_TIS = (Double_t)D_L0ElectronDecision_TIS_bool;
		D_L0PhotonDecision_TIS = (Double_t)D_L0PhotonDecision_TIS_bool;
		D_Hlt1TrackMVADecision_TOS = (Double_t)D_Hlt1TrackMVADecision_TOS_bool;
		D_Hlt1TwoTrackMVADecision_TOS = (Double_t)D_Hlt1TwoTrackMVADecision_TOS_bool;
        D_PT = D_PT_double; 
        D_ReFit_P = (Double_t)D_ReFit_P_list[0]; 
        D_ReFit_chi2 = (Double_t)D_ReFit_chi2_list[0]; 
        KS_PT = (Double_t)KS_PT_list[0]; 
        D_ReFit_KS0_P = (Double_t)D_ReFit_KS0_P_list[0]; 
		P0_ProbNNpi = P0_ProbNNpi_double;
		P1_ProbNNpi = P1_ProbNNpi_double;
		P2_ProbNNpi = P2_ProbNNpi_double;
        D_ReFit_decayLength_vs_Err = (Double_t)D_ReFit_decayLength_list[0]/(Double_t)D_ReFit_decayLengthErr_list[0]; 
		LOG_D_IPCHI2_OWNPV = TMath::Log(D_IPCHI2_OWNPV_double);

		if (((D_L0HadronDecision_TOS_bool==1) || (D_L0MuonDecision_TIS_bool==1) || (D_L0DiMuonDecision_TIS_bool==1) || 
				(D_L0ElectronDecision_TIS_bool==1) || (D_L0PhotonDecision_TIS_bool==1) ) && ( (D_Hlt1TrackMVADecision_TOS_bool==1) || 
				(D_Hlt1TwoTrackMVADecision_TOS_bool==1)) && (D_M_double>1779.65 && D_M_double<1959.65)&&(D_IPCHI2_OWNPV_double>0)
				&&(!(IsNaN(KS_PT_list[0]) ) )
				&&(!(IsNaN(D_PT_double)))
				&&(!(IsNaN(D_ReFit_P_list[0])))
				&&(!(IsNaN(D_ReFit_KS0_P_list[0]))) 
				&&(!(IsNaN(P0_ProbNNpi_double))) 
				&&(!(IsNaN(P1_ProbNNpi_double))) 
				&&(!(IsNaN(P2_ProbNNpi_double)))
				&&(!(IsNaN(D_ReFit_decayLength_list[0])))
				&&(D_ReFit_decayLengthErr_list[0]!=0)
				 ) {
			ds->add(variables,1.0);
			nevents++;
		}
     } 
	//chain->Print();
	cout << "<<<<number of events on dataset:" << nevents <<endl;
	cout << "<<<<Importing complete. RooDataSet ds created!\n";




	//creat newtree	
	chainx->SetBranchStatus("*",0);
    chainx->SetBranchStatus("D_M",1);
    chainx->SetBranchStatus("D_L0HadronDecision_TOS",1);
    chainx->SetBranchStatus("D_L0MuonDecision_TIS",1);
    chainx->SetBranchStatus("D_L0DiMuonDecision_TIS",1);
    chainx->SetBranchStatus("D_L0ElectronDecision_TIS",1);
    chainx->SetBranchStatus("D_L0PhotonDecision_TIS",1);
    chainx->SetBranchStatus("D_Hlt1TrackMVADecision_TOS",1);
    chainx->SetBranchStatus("D_Hlt1TwoTrackMVADecision_TOS",1);
    chainx->SetBranchStatus("D_PT",1);
    chainx->SetBranchStatus("D_ReFit_P",1);
    chainx->SetBranchStatus("D_ReFit_chi2",1);
    chainx->SetBranchStatus("KS_PT",1);
    chainx->SetBranchStatus("D_ReFit_KS0_P",1);
    chainx->SetBranchStatus("D_ReFit_decayLength",1);
    chainx->SetBranchStatus("D_ReFit_decayLengthErr",1);
    chainx->SetBranchStatus("P0_ProbNNpi",1);
    chainx->SetBranchStatus("P1_ProbNNpi",1);
    chainx->SetBranchStatus("P2_ProbNNpi",1);
    chainx->SetBranchStatus("D_IPCHI2_OWNPV",1);
	chainx->SetBranchStatus("P0_ID",1);
	chainx->SetBranchStatus("P0_P",1);
	chainx->SetBranchStatus("P0_PX",1);
	chainx->SetBranchStatus("P0_PY",1);
	chainx->SetBranchStatus("P0_PZ",1);
	chainx->SetBranchStatus("P1_P",1);
	chainx->SetBranchStatus("P1_PX",1);
	chainx->SetBranchStatus("P1_PY",1);
	chainx->SetBranchStatus("P1_PZ",1);
	chainx->SetBranchStatus("P2_P",1);
	chainx->SetBranchStatus("P2_PX",1);
	chainx->SetBranchStatus("P2_PY",1);
	chainx->SetBranchStatus("P2_PZ",1);

	TTree *newtree = (TTree*)chainx->CopyTree("(   ( (D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
        (D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1) ) && ( (D_Hlt1TrackMVADecision_TOS==1) || \
        (D_Hlt1TwoTrackMVADecision_TOS==1) )   )&& (D_M>1779.65&&D_M<1959.65) &&(D_IPCHI2_OWNPV>0)&&( (KS_PT>0.)||(KS_PT<=0.) )\
        &&( (D_PT>0.)||(D_PT<=0.) )\
        &&( (D_ReFit_P>0.)||(D_ReFit_P<=0.) )\
        &&( (D_ReFit_KS0_P>0.)||(D_ReFit_KS0_P<=0.) )\
        &&( (P0_ProbNNpi>0.)||(P0_ProbNNpi<=0.) )\
        &&( (P1_ProbNNpi>0.)||(P1_ProbNNpi<=0.) )\
        &&( (P2_ProbNNpi>0.)||(P2_ProbNNpi<=0.) )\
        &&( (D_ReFit_decayLength>0.)||(D_ReFit_decayLength<=0.) )\
        &&( (D_ReFit_decayLengthErr>0.)||(D_ReFit_decayLengthErr<0.) )\
        ");
    //newtree->Print();


	//check if nentries_newtree & nevents matches
	Long64_t nentries_newtree=newtree->GetEntries();
	if(nentries_newtree!=nevents) { 
		cout << "nentries_newtree!=nevents" <<endl;
		return 0; 
		}



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
	//RooPolynomial poly("poly", "3-para-poly", D_M, RooArgList(x1), 0);
	RooExponential poly("poly", "3-para-poly", D_M, x1);

	RooRealVar gauss_poly_frac("gauss_poly_frac","fraction of signal",20000,500000) ;
	RooRealVar gauss_poly_frac_1("gauss_poly_frac_1","fraction of background",20000,500000) ;
	RooAddPdf event("even","event",RooArgList(gauss, poly), RooArgList(gauss_poly_frac, gauss_poly_frac_1)) ;


	// Construct plot frame in 'D_M'
	RooPlot* xframe = D_M.frame(Title("Event p.d.f.")) ;
	RooPlot* xframe_2 = D_M.frame(Title(" ")) ;
	xframe_2->SetYTitle("pull distribution");

	event.fitTo(*ds, NumCPU(10)) ;
	cout<<"fitting complete!\n";
	//print variable values
	cout<<"para list:\n";

	mean.Print() ; 
	//mean_1.Print() ;  
	sigma.Print() ;  
	//sigma_1.Print() ;  
	x1.Print();  
	x2.Print();  
	x3.Print();
	//gauss_frac.Print(); 
	gauss_poly_frac.Print();  
	gauss_poly_frac_1.Print(); 

	// Draw all frames on a canvas
	// Draw part 1
	ds->plotOn(xframe);
	event.plotOn(xframe) ;
	RooHist* hpull = xframe->pullHist() ;
	event.plotOn(xframe,Components(poly),LineColor(kRed),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(gauss_all),LineColor(kBlue),LineStyle(kDashed)) ;
	//event.plotOn(xframe,Components(gauss),LineColor(kBlue),LineStyle(kDashed)) ;
	//event.plotOn(xframe,Components(gauss_1),LineColor(kGreen),LineStyle(kDashed)) ;
	xframe_2->addPlotable(hpull, "P") ;
	
	TCanvas* c = new TCanvas("total_plot","total_plot",600,1200) ;
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

	c->Print("../plotdir/D2KSPiPiPi_LL/fit_mass_D2KSPiPiPi_LL.pdf");



	// Draw part 2
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *ds, &event, RooArgList(gauss_poly_frac, gauss_poly_frac_1) );

	std::cout << "Check SWeights:" << std::endl;
	std::cout << std::endl <<  "Yield of gauss_poly_frac is "
	<< gauss_poly_frac.getVal() << ".  From sWeights it is "
	<< sData->GetYieldFromSWeight("gauss_poly_frac") << std::endl;
	std::cout << std::endl <<  "Yield of gauss_poly_frac_1 is "
	<< gauss_poly_frac_1.getVal() << ".  From sWeights it is "
	<< sData->GetYieldFromSWeight("gauss_poly_frac_1") << std::endl;


	
    Double_t swei, bwei;
	swei= 1. ;
    bwei= 1. ;
    TBranch *swei_branch = newtree->Branch("swei", &swei, "swei/D");
    TBranch *bwei_branch = newtree->Branch("bwei", &bwei, "bwei/D");
	
	for(Int_t sdata_index = 0; sdata_index < nentries_newtree; sdata_index++)
	{
		newtree->GetEntry(sdata_index);
		swei=sData->GetSWeight(sdata_index, "gauss_poly_frac");
		bwei=sData->GetSWeight(sdata_index, "gauss_poly_frac_1");
		//cout<<"this is the "<<sdata_index<<"th item of weights"<<endl;
		//cout<<swei<<endl;
		//cout<<bwei<<endl;
		swei_branch->Fill();
		bwei_branch->Fill();

	}
	//newtree->Print();
	newfile->cd();
    newtree->Write("ReducedTree");

	// create weightfed data set
	RooDataSet * signalpart = new RooDataSet(ds->GetName(),ds->GetTitle(),ds,*ds->get(),0,"gauss_poly_frac_sw") ;
	RooDataSet * backgroundpart = new RooDataSet(ds->GetName(),ds->GetTitle(),ds,*ds->get(),0,"gauss_poly_frac_1_sw") ;

	//D_PT
	TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	RooPlot *frame = D_PT.frame(Range(0.,20000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
	RooPlot* frame1 = D_PT.frame(Range(0.,20000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_D_PT_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_P
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_P.frame(Range(0.,300000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_P.frame(Range(0.,300000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_D_ReFit_P_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_chi2
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_chi2.frame(Range(0.,100.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_chi2.frame(Range(0.,100.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_D_ReFit_chi2_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//KS_PT
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = KS_PT.frame(Range(0.,10000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = KS_PT.frame(Range(0.,10000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_KS_PT_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_KS0_P
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_KS0_P.frame(Range(0.,100000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_KS0_P.frame(Range(0.,100000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_D_ReFit_KS0_P_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//P1_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P1_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P1_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_P1_ProbNNpi_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//P2_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P2_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P2_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_P2_ProbNNpi_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//P0_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P0_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P0_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_P0_ProbNNpi_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_decayLength_vs_Err
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_decayLength_vs_Err.frame(Range(0.,300.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_decayLength_vs_Err.frame(Range(0.,300.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_D_ReFit_decayLength_vs_Err_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//LOG_D_IPCHI2_OWNPV
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = LOG_D_IPCHI2_OWNPV.frame(Range(-2.,8.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = LOG_D_IPCHI2_OWNPV.frame(Range(-2.,8.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSPiPiPi_LL/SPlot/splot_fit_LOG_D_IPCHI2_OWNPV_D2KSPiPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//release pointers
	delete xframe;
	delete xframe_2;
	delete chain;
	delete c;
	return 0;
}



int fit_M_D2KSKpPiPi_LL()
{

	cout << "<<<< define TChain, TFile ..."<<endl;

	const char* TChain_name="D2KSKpPiPi_LL/DecayTree";
	const char* TChain_file="/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root";
	TChain *chain= new TChain(TChain_name);
	chain->Add(TChain_file);
	TChain *chainx= new TChain(TChain_name);
	chainx->Add(TChain_file);
    TFile *newfile = new TFile("../saved/TMVA_root_files/reduced_LL_2.root","recreate");


	cout << "<<<< define variables..."<<endl;

	RooRealVar D_M("D_M", "D+ mass", 1779.65, 1959.65, "MeV") ;
	RooRealVar D_L0HadronDecision_TOS("D_L0HadronDecision_TOS","D_L0HadronDecision_TOS",0,1);
	RooRealVar D_L0MuonDecision_TIS("D_L0MuonDecision_TIS","D_L0MuonDecision_TIS",0,1);
	RooRealVar D_L0DiMuonDecision_TIS("D_L0DiMuonDecision_TIS","D_L0DiMuonDecision_TIS",0,1);
	RooRealVar D_L0ElectronDecision_TIS("D_L0ElectronDecision_TIS","D_L0ElectronDecision_TIS",0,1);
	RooRealVar D_L0PhotonDecision_TIS("D_L0PhotonDecision_TIS","D_L0PhotonDecision_TIS",0,1);
	RooRealVar D_Hlt1TrackMVADecision_TOS("D_Hlt1TrackMVADecision_TOS","D_Hlt1TrackMVADecision_TOS",0,1);
	RooRealVar D_Hlt1TwoTrackMVADecision_TOS("D_Hlt1TwoTrackMVADecision_TOS","D_Hlt1TwoTrackMVADecision_TOS",0,1);
	RooArgSet variables(D_M,D_L0HadronDecision_TOS, D_L0MuonDecision_TIS,D_L0DiMuonDecision_TIS,D_L0ElectronDecision_TIS,
	D_L0PhotonDecision_TIS,D_Hlt1TrackMVADecision_TOS,D_Hlt1TwoTrackMVADecision_TOS); 
	RooRealVar D_PT("D_PT", "D+ transverse momentum", 0., 1000000.) ;
	RooRealVar D_ReFit_P("D_ReFit_P_double", "D_ReFit_P", 0., 100000000.) ;
	RooRealVar D_ReFit_chi2("D_ReFit_chi2_double", "D_ReFit_chi2", 0., 100000.) ;
	RooRealVar KS_PT("KS_PT_double", "KS_PT",  -100000., 100000.) ;
	RooRealVar D_ReFit_KS0_P("D_ReFit_KS0_P_double", "D_ReFit_KS0_P", 0., 100000000.) ;
	RooRealVar P0_ProbNNpi("P0_ProbNNpi", "P0_ProbNNpi", 0., 100.) ;
	RooRealVar P1_ProbNNpi("P1_ProbNNpi", "P1_ProbNNpi", 0., 100.) ;
	RooRealVar P2_ProbNNpi("P2_ProbNNpi", "P2_ProbNNpi", 0., 100.) ;
	RooRealVar D_ReFit_decayLength_vs_Err("D_ReFit_decayLength_vs_Err", "D_ReFit_decayLength_vs_Err", 0, RooNumber::infinity()) ;
	RooRealVar LOG_D_IPCHI2_OWNPV("LOG_D_IPCHI2_OWNPV", "LOG_D_IPCHI2_OWNPV", 0., RooNumber::infinity());
	variables.add(RooArgSet(D_PT, D_ReFit_P, D_ReFit_chi2, KS_PT, D_ReFit_KS0_P));
	variables.add(RooArgSet(P0_ProbNNpi, P1_ProbNNpi, P2_ProbNNpi, D_ReFit_decayLength_vs_Err, LOG_D_IPCHI2_OWNPV));


	cout << "<<<< define values to write..."<<endl;

	Double_t D_M_double ;
	Double_t D_PT_double ;
	Double_t P0_ProbNNpi_double ;
	Double_t P1_ProbNNpi_double ;
	Double_t P2_ProbNNpi_double ;
	Double_t D_IPCHI2_OWNPV_double ;
	Bool_t D_L0HadronDecision_TOS_bool;
	Bool_t D_L0MuonDecision_TIS_bool;
	Bool_t D_L0DiMuonDecision_TIS_bool;
	Bool_t D_L0ElectronDecision_TIS_bool;
	Bool_t D_L0PhotonDecision_TIS_bool;
	Bool_t D_Hlt1TrackMVADecision_TOS_bool;
	Bool_t D_Hlt1TwoTrackMVADecision_TOS_bool;
	Float_t D_ReFit_P_list[1];
	Float_t D_ReFit_chi2_list[1];
	Float_t KS_PT_list[1];
	Float_t D_ReFit_KS0_P_list[1];
	Float_t D_ReFit_decayLength_list[1];
	Float_t D_ReFit_decayLengthErr_list[1];


	cout << "<<<< set branch address..."<<endl;

    chain->SetBranchAddress("D_M",&D_M_double); 
    chain->SetBranchAddress("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS_bool); 
    chain->SetBranchAddress("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0DiMuonDecision_TIS",&D_L0DiMuonDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0ElectronDecision_TIS",&D_L0ElectronDecision_TIS_bool); 
    chain->SetBranchAddress("D_L0PhotonDecision_TIS",&D_L0PhotonDecision_TIS_bool); 
    chain->SetBranchAddress("D_Hlt1TrackMVADecision_TOS",&D_Hlt1TrackMVADecision_TOS_bool); 
    chain->SetBranchAddress("D_Hlt1TwoTrackMVADecision_TOS",&D_Hlt1TwoTrackMVADecision_TOS_bool); 
    chain->SetBranchAddress("D_PT",&D_PT_double); 
    chain->SetBranchAddress("D_ReFit_P",&D_ReFit_P_list); 
    chain->SetBranchAddress("D_ReFit_chi2",&D_ReFit_chi2_list); 
    chain->SetBranchAddress("KS_PT",&KS_PT_list); 
    chain->SetBranchAddress("D_ReFit_KS0_P",&D_ReFit_KS0_P_list); 
    chain->SetBranchAddress("D_ReFit_decayLength",&D_ReFit_decayLength_list); 
    chain->SetBranchAddress("D_ReFit_decayLengthErr",&D_ReFit_decayLengthErr_list); 
    chain->SetBranchAddress("P0_ProbNNpi",&P0_ProbNNpi_double); 
    chain->SetBranchAddress("P1_ProbNNpi",&P1_ProbNNpi_double); 
    chain->SetBranchAddress("P2_ProbNNpi",&P2_ProbNNpi_double); 
    chain->SetBranchAddress("D_IPCHI2_OWNPV",&D_IPCHI2_OWNPV_double); 
	

	cout << "<<<< creat dataset and write..."<<endl;

	RooDataSet *ds=new RooDataSet("ds", "ds", variables);
    Long64_t nentries = chain->GetEntries(); 
    Long64_t nevents=0; //count number of dataset events
	cout << "<<<<Number of Entries in chain:"<<nentries <<endl;
	cout << "<<<< Start importing data to *ds....\n";

    for (Long64_t i=0;i<nentries;i++) 
    { 
        chain->GetEntry(i);
		
		if (i<100) {/*
			cout << D_M_double<<"\t"<<D_PT_double<<"\t"
			<<P0_ProbNNpi_double<<"\t"<<P1_ProbNNpi_double<<"\t"<<P2_ProbNNpi_double<<"\t"
			<<D_L0HadronDecision_TOS_bool<<"\t"<<D_L0MuonDecision_TIS_bool<<"\t"
			<<D_L0DiMuonDecision_TIS_bool<<"\t"<<D_L0ElectronDecision_TIS_bool<<"\t"
			<<D_L0PhotonDecision_TIS_bool<<"\t"<<D_Hlt1TrackMVADecision_TOS_bool<<"\t"
			<<D_Hlt1TwoTrackMVADecision_TOS_bool<<"\t"<<D_ReFit_P_list[0]<<"\t"
			<<D_ReFit_chi2_list[0]<<"\t"<<KS_PT_list[0]<<"\t"<<D_ReFit_KS0_P_list[0]<<"\t"
			<<D_ReFit_decayLength_list[0]<<"\t"<<D_ReFit_decayLengthErr_list[0]<<"\t"
			<<"\n";		*/
		}
		
		D_M = D_M_double;
		D_L0HadronDecision_TOS = (Double_t)D_L0HadronDecision_TOS_bool;
		D_L0MuonDecision_TIS = (Double_t)D_L0MuonDecision_TIS_bool;
		D_L0DiMuonDecision_TIS = (Double_t)D_L0DiMuonDecision_TIS_bool;
		D_L0ElectronDecision_TIS = (Double_t)D_L0ElectronDecision_TIS_bool;
		D_L0PhotonDecision_TIS = (Double_t)D_L0PhotonDecision_TIS_bool;
		D_Hlt1TrackMVADecision_TOS = (Double_t)D_Hlt1TrackMVADecision_TOS_bool;
		D_Hlt1TwoTrackMVADecision_TOS = (Double_t)D_Hlt1TwoTrackMVADecision_TOS_bool;
        D_PT = D_PT_double; 
        D_ReFit_P = (Double_t)D_ReFit_P_list[0]; 
        D_ReFit_chi2 = (Double_t)D_ReFit_chi2_list[0]; 
        KS_PT = (Double_t)KS_PT_list[0]; 
        D_ReFit_KS0_P = (Double_t)D_ReFit_KS0_P_list[0]; 
		P0_ProbNNpi = P0_ProbNNpi_double;
		P1_ProbNNpi = P1_ProbNNpi_double;
		P2_ProbNNpi = P2_ProbNNpi_double;
        D_ReFit_decayLength_vs_Err = (Double_t)D_ReFit_decayLength_list[0]/(Double_t)D_ReFit_decayLengthErr_list[0]; 
		LOG_D_IPCHI2_OWNPV = TMath::Log(D_IPCHI2_OWNPV_double);

		if (((D_L0HadronDecision_TOS_bool==1) || (D_L0MuonDecision_TIS_bool==1) || (D_L0DiMuonDecision_TIS_bool==1) || 
				(D_L0ElectronDecision_TIS_bool==1) || (D_L0PhotonDecision_TIS_bool==1) ) && ( (D_Hlt1TrackMVADecision_TOS_bool==1) || 
				(D_Hlt1TwoTrackMVADecision_TOS_bool==1)) && (D_M_double>1779.65 && D_M_double<1959.65)&&(D_IPCHI2_OWNPV_double>0) ) {
			ds->add(variables,1.0);
			nevents++;
		}
     } 
	//chain->Print();
	cout << "<<<<number of events on dataset:" << nevents <<endl;
	cout << "******Importing complete. RooDataSet ds created!\n";


	//creat newtree	
	chainx->SetBranchStatus("*",0);
    chainx->SetBranchStatus("D_M",1);
    chainx->SetBranchStatus("D_L0HadronDecision_TOS",1);
    chainx->SetBranchStatus("D_L0MuonDecision_TIS",1);
    chainx->SetBranchStatus("D_L0DiMuonDecision_TIS",1);
    chainx->SetBranchStatus("D_L0ElectronDecision_TIS",1);
    chainx->SetBranchStatus("D_L0PhotonDecision_TIS",1);
    chainx->SetBranchStatus("D_Hlt1TrackMVADecision_TOS",1);
    chainx->SetBranchStatus("D_Hlt1TwoTrackMVADecision_TOS",1);
    chainx->SetBranchStatus("D_PT",1);
    chainx->SetBranchStatus("D_ReFit_P",1);
    chainx->SetBranchStatus("D_ReFit_chi2",1);
    chainx->SetBranchStatus("KS_PT",1);
    chainx->SetBranchStatus("D_ReFit_KS0_P",1);
    chainx->SetBranchStatus("D_ReFit_decayLength",1);
    chainx->SetBranchStatus("D_ReFit_decayLengthErr",1);
    chainx->SetBranchStatus("P0_ProbNNpi",1);
    chainx->SetBranchStatus("P1_ProbNNpi",1);
    chainx->SetBranchStatus("P2_ProbNNpi",1);
    chainx->SetBranchStatus("D_IPCHI2_OWNPV",1);
	chainx->SetBranchStatus("P0_ID",1);
	chainx->SetBranchStatus("P0_P",1);
	chainx->SetBranchStatus("P0_PX",1);
	chainx->SetBranchStatus("P0_PY",1);
	chainx->SetBranchStatus("P0_PZ",1);
	chainx->SetBranchStatus("P1_P",1);
	chainx->SetBranchStatus("P1_PX",1);
	chainx->SetBranchStatus("P1_PY",1);
	chainx->SetBranchStatus("P1_PZ",1);
	chainx->SetBranchStatus("P2_P",1);
	chainx->SetBranchStatus("P2_PX",1);
	chainx->SetBranchStatus("P2_PY",1);
	chainx->SetBranchStatus("P2_PZ",1);
	TTree *newtree = (TTree*)chainx->CopyTree("(   ( (D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
        (D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1) ) && ( (D_Hlt1TrackMVADecision_TOS==1) || \
        (D_Hlt1TwoTrackMVADecision_TOS==1) )   )&& (D_M>1779.65&&D_M<1959.65) &&(D_IPCHI2_OWNPV>0)");
    //newtree->Print();


	//check if nentries_newtree & nevents matches
	Long64_t nentries_newtree=newtree->GetEntries();
	if(nentries_newtree!=nevents) { 
		cout << "nentries_newtree!=nevents" <<endl;
		return 0; 
		}


	RooRealVar mean("mean","mean of gaussian",1869.,1859.,1879.) ;
	RooRealVar sigma("sigma","width of gaussian",10.,0.1,20.) ;
	RooRealVar alpha_var("alpha_var","alpha_var",1.4,1.0,10.) ;
	RooRealVar n_var("n_var","n_var",30.7631, 0. , 100.) ;
	//RooCBShape gauss_all("CBfunction", "CBfunction", D_M, mean, sigma, alpha_var, n_var);
	RooGaussian gauss_all("Gfunction", "Gfunction", D_M, mean, sigma);

	RooRealVar x1("x1", "para1", -100., 100.);
	RooRealVar x2("x2", "para2", -100., 100.);
	RooRealVar x3("x3", "para3", -100., 100.);
	RooChebychev poly("poly", "3-para-poly", D_M, RooArgList(x1, x2, x3));
	//num total:94333
	RooRealVar gauss_poly_frac("gauss_poly_frac","fraction of component 1 in signal",1000., 200000) ;
	RooRealVar gauss_poly_frac_1("gauss_poly_frac_1","fraction of component 2 in signal",10000., 2000000) ;
	RooAddPdf event("even","event",RooArgList(gauss_all, poly), RooArgList(gauss_poly_frac, gauss_poly_frac_1)) ;


	// Construct plot frame in 'D_M'
	RooPlot* xframe = D_M.frame(Title("Event p.d.f.")) ;
	RooPlot* xframe_2 = D_M.frame(Title(" ")) ;
	xframe_2->SetYTitle("pull distribution");

	event.fitTo(*ds) ;
	cout<<"fitting complete!\n";
	//print variable values
	cout<<"para list:\n";
	mean.Print() ;  sigma.Print() ;  alpha_var.Print() ;  n_var.Print() ;  x1.Print();  x2.Print();  x3.Print();
	gauss_poly_frac.Print(); 
	gauss_poly_frac_1.Print(); 

	// Draw all frames on a canvas
	// Draw part 1
	ds->plotOn(xframe);
	event.plotOn(xframe) ;
	RooHist* hpull = xframe->pullHist() ;
	event.plotOn(xframe,Components(poly),LineColor(kRed),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(gauss_all),LineColor(kBlue),LineStyle(kDashed)) ;
	xframe_2->addPlotable(hpull, "P") ;
	
	TCanvas* c = new TCanvas("total_plot","total_plot",600,1200) ;
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

	c->Print("../plotdir/D2KSKpPiPi_LL/fit_mass_D2KSKpPiPi_LL.pdf");



	// Draw part 2
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *ds, &event, RooArgList(gauss_poly_frac, gauss_poly_frac_1) );

	std::cout << "Check SWeights:" << std::endl;
	std::cout << std::endl <<  "Yield of gauss_poly_frac is "
	<< gauss_poly_frac.getVal() << ".  From sWeights it is "
	<< sData->GetYieldFromSWeight("gauss_poly_frac") << std::endl;
	std::cout << std::endl <<  "Yield of gauss_poly_frac_1 is "
	<< gauss_poly_frac_1.getVal() << ".  From sWeights it is "
	<< sData->GetYieldFromSWeight("gauss_poly_frac_1") << std::endl;


	
    Double_t swei, bwei;
	swei= 1. ;
    bwei= 1. ;
    TBranch *swei_branch = newtree->Branch("swei", &swei, "swei/D");
    TBranch *bwei_branch = newtree->Branch("bwei", &bwei, "bwei/D");
	
	for(Int_t sdata_index = 0; sdata_index < nentries_newtree; sdata_index++)
	{
		newtree->GetEntry(sdata_index);
		swei=sData->GetSWeight(sdata_index, "gauss_poly_frac");
		bwei=sData->GetSWeight(sdata_index, "gauss_poly_frac_1");
		//cout<<"this is the "<<sdata_index<<"th item of weights"<<endl;
		//cout<<swei<<endl;
		//cout<<bwei<<endl;
		swei_branch->Fill();
		bwei_branch->Fill();

	}
	//newtree->Print();
	newfile->cd();
    newtree->Write("ReducedTree");

	// create weightfed data set
	RooDataSet * signalpart = new RooDataSet(ds->GetName(),ds->GetTitle(),ds,*ds->get(),0,"gauss_poly_frac_sw") ;
	RooDataSet * backgroundpart = new RooDataSet(ds->GetName(),ds->GetTitle(),ds,*ds->get(),0,"gauss_poly_frac_1_sw") ;

	//D_PT
	TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	RooPlot *frame = D_PT.frame(Range(0.,20000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
	RooPlot* frame1 = D_PT.frame(Range(0.,20000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_D_PT_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_P
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_P.frame(Range(0.,300000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_P.frame(Range(0.,300000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_D_ReFit_P_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_chi2
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_chi2.frame(Range(0.,100.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_chi2.frame(Range(0.,100.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_D_ReFit_chi2_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//KS_PT
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = KS_PT.frame(Range(0.,10000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = KS_PT.frame(Range(0.,10000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_KS_PT_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_KS0_P
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_KS0_P.frame(Range(0.,100000.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_KS0_P.frame(Range(0.,100000.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_D_ReFit_KS0_P_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//P1_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P1_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P1_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_P1_ProbNNpi_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//P2_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P2_ProbNNpi.frame(Range(0.95,1.01), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P2_ProbNNpi.frame(Range(0.95,1.01), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_P2_ProbNNpi_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//P0_ProbNNpi
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = P0_ProbNNpi.frame(Range(0.,0.2), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = P0_ProbNNpi.frame(Range(0.,0.2), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_P0_ProbNNpi_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//D_ReFit_decayLength_vs_Err
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = D_ReFit_decayLength_vs_Err.frame(Range(0.,300.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = D_ReFit_decayLength_vs_Err.frame(Range(0.,300.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_D_ReFit_decayLength_vs_Err_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//LOG_D_IPCHI2_OWNPV
	cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
	cdata->Divide(1,2);
	frame = LOG_D_IPCHI2_OWNPV.frame(Range(-2.,8.), Title("signal"));
	signalpart->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(1);  frame->Draw();
    frame1 = LOG_D_IPCHI2_OWNPV.frame(Range(-2.,8.), Title("background")) ;
	backgroundpart->plotOn(frame1, DataError(RooAbsData::SumW2) ) ;
	cdata->cd(2);  frame1->Draw();
	cdata->Print("../plotdir/D2KSKpPiPi_LL/SPlot/splot_fit_LOG_D_IPCHI2_OWNPV_D2KSKpPiPi_LL.pdf");
	delete cdata; delete frame; delete frame1;

	//release pointers
	delete xframe;
	delete xframe_2;
	delete chain;
	delete c;
	return 0;
}














int main()
{
    /* code */
    int i1;
    int i2;
	int i3;
	int i4;
    cout << "Fit subtr_D2KSPiPiPi_DD? ('1' for yes, others for no)"<<"\n";
    cin >> i1;
    cout << "Fit subtr_D2KSKpPiPi_DD? ('1' for yes, others for no)"<<"\n";
    cin >> i2;
	cout << "Fit subtr_D2KSPiPiPi_LL? ('1' for yes, others for no)"<<"\n";
    cin >> i3;
    cout << "Fit subtr_D2KSKpPiPi_LL? ('1' for yes, others for no)"<<"\n";
    cin >> i4;
    
    if (i1==1) {
        fit_M_D2KSPiPiPi_DD();
    }
    
    if (i2==1) {
        fit_M_D2KSKpPiPi_DD();
    }

	if (i3==1) {
        fit_M_D2KSPiPiPi_LL();
    }
    
    if (i4==1) {
        fit_M_D2KSKpPiPi_LL();
    }
    return 0;
}
