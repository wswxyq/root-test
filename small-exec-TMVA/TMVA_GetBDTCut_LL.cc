#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TVector3.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
using namespace TMVA;
using namespace std;
using namespace RooFit ;
void TMVA_GetBDTCut_LL()
{

    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    //set origin file used to train
    TString fname = "../saved/TMVA_root_files/reduced_LL_1.root";
    TFile * input = new TFile( fname ); 


    Float_t var[8];
    //Define input variables
    reader->AddVariable( "D_PT", &var[0] );
    reader->AddVariable( "P1_ProbNNpi", &var[2] );
    reader->AddVariable( "P2_ProbNNpi", &var[3] );
    reader->AddVariable( "D_ReFit_chi2", &var[4] );
    reader->AddVariable( "KS_PT", &var[7] );
    reader->AddVariable( "var1:=D_ReFit_decayLength/D_ReFit_decayLengthErr", &var[5] );
	reader->AddVariable( "var2:=log(D_IPCHI2_OWNPV)", &var[6] );	


    // Book the MVA methods
    TString dir    = "D2KSPiPiPi_LL/weights/";
    TString prefix = "TMVAClassification";
    TString methodName = TString("BDT") + TString(" method");
    TString weightfile = dir + prefix + TString("_") + TString("BDT") + TString(".weights.xml");
    reader->BookMVA( methodName, weightfile );


    // Book output histograms
    UInt_t nbin = 100;
    TH1F   *histBdt(0);
    histBdt     = new TH1F( "MVA_BDT",  "MVA_BDT",  nbin, -0.8, 0.8 );


    // Event loop
    TTree* theTree = (TTree*)input->Get("ReducedTree");
    // output files
    TFile *target  = new TFile( "../saved/TMVA_root_files/TMVApp_LLx.root","RECREATE" );
    TTree *tree = theTree->CopyTree("");    //use this tree to calculate BDT cut
    
    Double_t BDT_response;
    TBranch *BDT_response_branch = tree->Branch("BDT_response", &BDT_response, "BDT_response/D");
    Double_t Ct;
    TBranch *Ct_branch = tree->Branch("Ct", &Ct, "Ct/D");

    Double_t userVard[9];
    Float_t userVarf[8];


    theTree->SetBranchAddress( "D_PT", &userVard[0] );
    theTree->SetBranchAddress( "P0_ProbNNpi", &userVard[1] );
    theTree->SetBranchAddress( "P1_ProbNNpi", &userVard[2] );
    theTree->SetBranchAddress( "P2_ProbNNpi", &userVard[3] );
    theTree->SetBranchAddress( "D_ReFit_chi2", &userVarf[0] );
    theTree->SetBranchAddress( "D_ReFit_decayLength", &userVarf[2] );
    theTree->SetBranchAddress( "D_ReFit_decayLengthErr", &userVarf[3] );
	theTree->SetBranchAddress( "D_IPCHI2_OWNPV", &userVard[4] );	
    theTree->SetBranchAddress( "KS_PT", &userVard[5] );




    TStopwatch sw;
    // calculate BDT for tree
    sw.Start();
    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
        //if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
        

        
        theTree->GetEntry(ievt);
        var[0]=userVard[0];
        var[1]=userVard[1];
        var[2]=userVard[2];
        var[3]=userVard[3];
        var[4]=userVarf[0];
        var[5]=userVarf[2]/userVarf[3];
        var[6]=TMath::Log(userVard[4]);
        var[7]=userVard[5];
        // Return the MVA outputs and fill into histograms
        histBdt    ->Fill( reader->EvaluateMVA( "BDT method" ) );

        BDT_response=reader->EvaluateMVA( "BDT method" );
        BDT_response_branch->Fill();


        if (ievt<100) {
            cout<<ievt<<"  "<<var[1]<<"  "<<var[2]<<"  "<<var[3]<<"  "<<var[4]<<endl;
            cout<<var[5]<<"  "<<var[6]<<"  "<<var[7]<<endl;

        }

        
    }

    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();


    target->cd();
    tree->Write("MVATree");
    histBdt ->Write();


    //calculate significance
    
    Double_t cut_up=tree->GetMaximum("BDT_response");
    Double_t cut_down=tree->GetMinimum("BDT_response");

    Long64_t num_cut=50;
    Double_t cut_value[ num_cut+1 ], significance_value[ num_cut+1 ];


    Double_t kk_signal_yield=1105.77;
    Double_t kk_background_yield=134890;
    Double_t kp_signal_yield=21980.1; //to fit
    Double_t kp_background_yield=400535.; //to fit
    Double_t temp_sig=0.;
    Double_t temp_bkg=0.;

	RooRealVar D_M("D_M", "D+ mass", 1779.65, 1959.65, "MeV") ;
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

	RooRealVar gauss_poly_frac("gauss_poly_frac","fraction of signal",0.,500000.) ;
	RooRealVar gauss_poly_frac_1("gauss_poly_frac_1","fraction of background",0.,500000.) ;
	RooAddPdf event("even","event",RooArgList(gauss, poly), RooArgList(gauss_poly_frac, gauss_poly_frac_1)) ;

    RooRealVar BDT_value("BDT_response", "BDT_response", -1000, 1000);
	RooArgSet variables(D_M, BDT_value); 
    

    for(Long64_t i = 0; i < num_cut; i++)
    {
        cut_value[i] = cut_down + (cut_up-cut_down)*i/num_cut;
        string bdtcut="BDT_response>";
        bdtcut+=to_string(cut_value[i]);
        
        RooDataSet *ds=new RooDataSet("ds_plus", "ds_plus", variables, Import(*tree), Cut(bdtcut.c_str()));
        event.fitTo(*ds, NumCPU(10));
        temp_sig=gauss_poly_frac.getVal()*kk_signal_yield/kp_signal_yield;
        temp_bkg=gauss_poly_frac_1.getVal()*kk_background_yield/kp_background_yield;

        significance_value[i] = temp_sig/sqrt(temp_sig+temp_bkg);

        //std::cout<<cut_value[i]<<"   "<<sum_sweight<<"   "<<sum_bweight
        //<<" "<<significance_value[i]<<std::endl;
        
    }


    
    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    c1->SetFillColor(42);
    c1->SetGrid();

    TGraph *gr = new TGraph(num_cut, cut_value, significance_value);
    gr->SetLineColor(2);
    gr->SetLineWidth(4);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->SetTitle("S vs Response");
    gr->GetXaxis()->SetTitle("BDT cut");
    gr->GetYaxis()->SetTitle("Significance");
    gr->Draw("ACP");

    // TCanvas::Update() draws the frame, after which one can change it
    c1->Update();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    c1->Modified();
    c1->Print("graph_LLx.eps");
    c1->Print("graph_LLx.png");

    target  ->Close();
    delete reader;

    Double_t bestcut;
    bestcut=cut_value[distance(significance_value, max_element(significance_value, significance_value + num_cut))];
    cout    << "The best BDT cut value: "<<  bestcut  << endl;
    
    for(int ii = 0; ii < num_cut; ii++)
    {
        cout<<cut_value[ii]<<"  "<<significance_value[ii]<<endl;
    }



    }

