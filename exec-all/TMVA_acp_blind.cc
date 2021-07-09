#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <fstream>
#include <map>
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

void TMVA_acp_blind() {

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
    const char* TChain_DD_file="./DD/TMVA_output_DD.root";
    const char* TChain_LL_file="./LL/TMVA_output_LL.root";
    TChain *chain_DD= new TChain(TChain_name);
    TChain *chain_LL= new TChain(TChain_name);
    chain_DD->Add(TChain_DD_file);
    chain_LL->Add(TChain_LL_file);


    //define variables in dataset
    RooRealVar D_M("D_M", "D+ mass", 1799.65, 1939.65, "MeV") ;
    RooRealVar BDT_response("BDT_response", "BDT_response value", 0., RooNumber::infinity());
    RooRealVar Ct("Ct", "Ct value", -RooNumber::infinity(), RooNumber::infinity());
    RooRealVar P0_ID("P0_ID", "P0_ID value", -1000., 1000.);
    RooArgSet variables(D_M, BDT_response, Ct, P0_ID);


    //creat dataset and write
    RooDataSet *DD_plus_plus=new RooDataSet("DD_plus_plus", "DD_plus_plus", variables, Import(*chain_DD),
                                            Cut(" P0_ID>0 && Ct>0 "));
    RooDataSet *DD_plus_minus=new RooDataSet("DD_plus_minus", "DD_plus_minus", variables, Import(*chain_DD),
            Cut(" P0_ID>0 && Ct<0 "));
    RooDataSet *DD_minus_plus=new RooDataSet("DD_minus_plus", "DD_minus_plus", variables, Import(*chain_DD),
            Cut(" P0_ID<0 && Ct>0 "));
    RooDataSet *DD_minus_minus=new RooDataSet("DD_minus_minus", "DD_minus_minus", variables, Import(*chain_DD),
            Cut(" P0_ID<0 && Ct<0 "));

    RooDataSet *LL_plus_plus=new RooDataSet("LL_plus_plus", "LL_plus_plus", variables, Import(*chain_LL),
                                            Cut(" P0_ID>0 && Ct>0 "));
    RooDataSet *LL_plus_minus=new RooDataSet("LL_plus_minus", "LL_plus_minus", variables, Import(*chain_LL),
            Cut(" P0_ID>0 && Ct<0 "));
    RooDataSet *LL_minus_plus=new RooDataSet("LL_minus_plus", "LL_minus_plus", variables, Import(*chain_LL),
            Cut(" P0_ID<0 && Ct>0 "));
    RooDataSet *LL_minus_minus=new RooDataSet("LL_minus_minus", "LL_minus_minus", variables, Import(*chain_LL),
            Cut(" P0_ID<0 && Ct<0 "));

    //pdf for D mass spectrum
    RooRealVar mean_DD_plus("mean_DD_plus","mean_DD_plus of gaussian", 1869.6) ;
    RooRealVar mean_1_DD_plus("mean_1_DD_plus","mean_1_DD_plus of gaussian(1)", 1870.03) ;
    RooRealVar sigma_DD_plus("sigma_DD_plus","width of gaussian", 12.6233) ;
    RooRealVar sigma_1_DD_plus("sigma_1_DD_plus","width of gaussian(1)", 7.30743) ;
    RooGaussian gauss_DD_plus("gauss_DD_plus","gaussian PDF", D_M, mean_DD_plus, sigma_DD_plus) ;
    RooGaussian gauss_1_DD_plus("gauss_1_DD_plus","gaussian PDF(1)", D_M, mean_1_DD_plus, sigma_1_DD_plus) ;
    RooRealVar gauss_frac_DD_plus("gauss_frac_DD_plus","fraction of component 1 in signal",0.322784) ;
    RooAddPdf gauss_all_DD_plus("gauss_all_DD_plus", "combination of gaussian function", RooArgList(gauss_DD_plus, gauss_1_DD_plus),
                                RooArgList(gauss_frac_DD_plus));
    RooRealVar x1_DD_plus("x1_DD_plus", "para1", 0.0050725);
    RooExponential poly_DD_plus("poly_DD_plus", "3-para-poly", D_M, x1_DD_plus);


    RooRealVar mean_DD_minus("mean_DD_minus","mean of gaussian", 1869.53) ;
    RooRealVar mean_1_DD_minus("mean_1_DD_minus","mean of gaussian(1)", 1870.12) ;
    RooRealVar sigma_DD_minus("sigma_DD_minus","width of gaussian", 11.6634) ;
    RooRealVar sigma_1_DD_minus("sigma_1_DD_minus","width of gaussian(1)", 6.96087) ;
    RooGaussian gauss_DD_minus("gauss_DD_minus","gaussian PDF", D_M, mean_DD_minus, sigma_DD_minus) ;
    RooGaussian gauss_1_DD_minus("gauss_1_DD_minus","gaussian PDF(1)", D_M, mean_1_DD_minus, sigma_1_DD_minus) ;
    RooRealVar gauss_frac_DD_minus("gauss_frac_DD_minus","fraction of component 1 in signal", 0.431087) ;
    RooAddPdf gauss_all_DD_minus("gauss_all_DD_minus", "combination of gaussian function", RooArgList(gauss_DD_minus, gauss_1_DD_minus),
                                 RooArgList(gauss_frac_DD_minus));
    RooRealVar x1_DD_minus("x1_DD_minus", "para1", 0.00519574);
    RooExponential poly_DD_minus("poly_DD_minus", "3-para-poly", D_M, x1_DD_minus);


    RooRealVar mean_LL_plus("mean_LL_plus","mean of gaussian", 1870.15) ;
    RooRealVar mean_1_LL_plus("mean_1_LL_plus","mean of gaussian(1)", 1863.72) ;
    RooRealVar sigma_LL_plus("sigma_LL_plus","width of gaussian", 6.78022) ;
    RooRealVar sigma_1_LL_plus("sigma_1_LL_plus","width of gaussian(1)", 29.5862) ;
    RooGaussian gauss_LL_plus("gauss_LL_plus","gaussian PDF", D_M, mean_LL_plus, sigma_LL_plus) ;
    RooGaussian gauss_1_LL_plus("gauss_1_LL_plus","gaussian PDF(1)", D_M, mean_1_LL_plus, sigma_1_LL_plus) ;
    RooRealVar gauss_frac_LL_plus("gauss_frac_LL_plus","fraction of component 1 in signal", 0.792908) ;
    RooAddPdf gauss_all_LL_plus("gauss_all_LL_plus", "combination of gaussian function", RooArgList(gauss_LL_plus, gauss_1_LL_plus),
                                RooArgList(gauss_frac_LL_plus));
    RooRealVar x1_LL_plus("x1_LL_plus", "para1", 0.00442008);
    RooExponential poly_LL_plus("poly_LL_plus", "3-para-poly", D_M, x1_LL_plus);


    RooRealVar mean_LL_minus("mean_LL_minus","mean of gaussian", 1870.04) ;
    RooRealVar mean_1_LL_minus("mean_1_LL_minus","mean of gaussian(1)", 1868.37) ;
    RooRealVar sigma_LL_minus("sigma_LL_minus","width of gaussian", 6.42809) ;
    RooRealVar sigma_1_LL_minus("sigma_1_LL_minus","width of gaussian(1)", 17.8949) ;
    RooGaussian gauss_LL_minus("gauss_LL_minus","gaussian PDF", D_M, mean_LL_minus, sigma_LL_minus) ;
    RooGaussian gauss_1_LL_minus("gauss_1_LL_minus","gaussian PDF(1)", D_M, mean_1_LL_minus, sigma_1_LL_minus) ;
    RooRealVar gauss_frac_LL_minus("gauss_frac_LL_minus","fraction of component 1 in signal", 0.781868) ;
    RooAddPdf gauss_all_LL_minus("gauss_all_LL_minus", "combination of gaussian function", RooArgList(gauss_LL_minus, gauss_1_LL_minus),
                                 RooArgList(gauss_frac_LL_minus));
    RooRealVar x1_LL_minus("x1_LL_minus", "para1", 0.00428315);
    RooExponential poly_LL_minus("poly_LL_minus", "3-para-poly", D_M, x1_LL_minus);


    //global--------------------------------------------------------------
    RooRealVar AT_plus("AT_plus", "AT_plus", 0, -0.1, 0.1);
    RooUnblindPrecision AT_plus_unblind("AT_plus_unblind", "AT_plus_unblind", "Atblind", 0., 0.1, AT_plus, blindCat, 1);
    RooRealVar AT_minus("AT_minus", "AT_minus", 0, -0.1, 0.1);
    RooUnblindPrecision AT_minus_unblind("AT_minus_unblind", "AT_minus_unblind", "Atbarblind", 0., 0.1, AT_minus, blindCat, 1);

    //DD------------------------------------------------------------------
    RooRealVar Ntot_DD_plus("Ntot_DD_plus", "Ntot_DD_plus", 0., 500000.); //number of D+
    RooRealVar Ntot_DD_minus("Ntot_DD_minus", "Ntot_DD_minus", 0., 500000.); //number od D-

    RooRealVar Ntot_DD_bkg_plus("Ntot_DD_bkg_plus", "Ntot_DD_bkg_plus", 0., 500000.);
    RooRealVar Ntot_DD_bkg_minus("Ntot_DD_bkg_minus", "Ntot_DD_bkg_minus", 0., 500000.);

    RooRealVar AT_DD_bkg_plus("AT_DD_bkg_plus", "AT_DD_bkg_plus", -0.1, 0.1);
    RooRealVar AT_DD_bkg_minus("AT_DD_bkg_minus", "AT_DD_bkg_minus", -0.1, 0.1);

    RooFormulaVar Nsig_DD_plus_plus("Nsig_DD_plus_plus", "0.5*@0*(1+@1)", RooArgList(Ntot_DD_plus, AT_plus_unblind));
    RooFormulaVar Nsig_DD_plus_minus("Nsig_DD_plus_minus", "0.5*@0*(1-@1)", RooArgList(Ntot_DD_plus, AT_plus_unblind));
    RooFormulaVar Nsig_DD_minus_minus("Nsig_DD_minus_minus", "0.5*@0*(1+@1)", RooArgList(Ntot_DD_minus, AT_minus_unblind));
    RooFormulaVar Nsig_DD_minus_plus("Nsig_DD_minus_plus", "0.5*@0*(1-@1)", RooArgList(Ntot_DD_minus, AT_minus_unblind));

    RooFormulaVar Nbkg_DD_plus_plus("Nbkg_DD_plus_plus", "0.5*@0*(1+@1)", RooArgList(Ntot_DD_bkg_plus, AT_DD_bkg_plus));
    RooFormulaVar Nbkg_DD_plus_minus("Nbkg_DD_plus_minus", "0.5*@0*(1-@1)", RooArgList(Ntot_DD_bkg_plus, AT_DD_bkg_plus));
    RooFormulaVar Nbkg_DD_minus_minus("Nbkg_DD_minus_minus", "0.5*@0*(1+@1)", RooArgList(Ntot_DD_bkg_minus, AT_DD_bkg_minus));
    RooFormulaVar Nbkg_DD_minus_plus("Nbkg_DD_minus_plus", "0.5*@0*(1-@1)", RooArgList(Ntot_DD_bkg_minus, AT_DD_bkg_minus));

    RooAddPdf  model_DD_plus_plus("model_DD_plus_plus","model_DD_plus_plus",
                                  RooArgList(gauss_all_DD_plus, poly_DD_plus), RooArgList(Nsig_DD_plus_plus, Nbkg_DD_plus_plus)) ;
    RooAddPdf  model_DD_minus_plus("model_DD_minus_plus","model_DD_minus_plus",
                                   RooArgList(gauss_all_DD_minus, poly_DD_minus), RooArgList(Nsig_DD_minus_plus, Nbkg_DD_minus_plus)) ;
    RooAddPdf  model_DD_plus_minus("model_DD_plus_minus","model_DD_plus_minus",
                                   RooArgList(gauss_all_DD_plus, poly_DD_plus), RooArgList(Nsig_DD_plus_minus, Nbkg_DD_plus_minus)) ;
    RooAddPdf  model_DD_minus_minus("model_DD_minus_minus","model_DD_minus_minus",
                                    RooArgList(gauss_all_DD_minus, poly_DD_minus), RooArgList(Nsig_DD_minus_minus, Nbkg_DD_minus_minus)) ;



    //LL------------------------------------------------------------------
    RooRealVar Ntot_LL_plus("Ntot_LL_plus", "Ntot_LL_plus", 0., 1000000.); //number of D+
    RooRealVar Ntot_LL_minus("Ntot_LL_minus", "Ntot_LL_minus", 0., 1000000.); //number od D-

    RooRealVar Ntot_LL_bkg_plus("Ntot_LL_bkg_plus", "Ntot_LL_bkg_plus", 0., 1000000.);
    RooRealVar Ntot_LL_bkg_minus("Ntot_LL_bkg_minus", "Ntot_LL_bkg_minus", 0., 1000000.);

    RooRealVar AT_LL_bkg_plus("AT_LL_bkg_plus", "AT_LL_bkg_plus", -0.1, 0.1);
    RooRealVar AT_LL_bkg_minus("AT_LL_bkg_minus", "AT_LL_bkg_minus", -0.1, 0.1);

    RooFormulaVar Nsig_LL_plus_plus("Nsig_LL_plus_plus", "0.5*@0*(1+@1)", RooArgList(Ntot_LL_plus, AT_plus_unblind));
    RooFormulaVar Nsig_LL_plus_minus("Nsig_LL_plus_minus", "0.5*@0*(1-@1)", RooArgList(Ntot_LL_plus, AT_plus_unblind));
    RooFormulaVar Nsig_LL_minus_minus("Nsig_LL_minus_minus", "0.5*@0*(1+@1)", RooArgList(Ntot_LL_minus, AT_minus_unblind));
    RooFormulaVar Nsig_LL_minus_plus("Nsig_LL_minus_plus", "0.5*@0*(1-@1)", RooArgList(Ntot_LL_minus, AT_minus_unblind));

    RooFormulaVar Nbkg_LL_plus_plus("Nbkg_LL_plus_plus", "0.5*@0*(1+@1)", RooArgList(Ntot_LL_bkg_plus, AT_LL_bkg_plus));
    RooFormulaVar Nbkg_LL_plus_minus("Nbkg_LL_plus_minus", "0.5*@0*(1-@1)", RooArgList(Ntot_LL_bkg_plus, AT_LL_bkg_plus));
    RooFormulaVar Nbkg_LL_minus_minus("Nbkg_LL_minus_minus", "0.5*@0*(1+@1)", RooArgList(Ntot_LL_bkg_minus, AT_LL_bkg_minus));
    RooFormulaVar Nbkg_LL_minus_plus("Nbkg_LL_minus_plus", "0.5*@0*(1-@1)", RooArgList(Ntot_LL_bkg_minus, AT_LL_bkg_minus));

    RooAddPdf  model_LL_plus_plus("model_LL_plus_plus","model_LL_plus_plus",
                                  RooArgList(gauss_all_LL_plus, poly_LL_plus), RooArgList(Nsig_LL_plus_plus, Nbkg_LL_plus_plus)) ;
    RooAddPdf  model_LL_minus_plus("model_LL_minus_plus","model_LL_minus_plus",
                                   RooArgList(gauss_all_LL_minus, poly_LL_minus), RooArgList(Nsig_LL_minus_plus, Nbkg_LL_minus_plus)) ;
    RooAddPdf  model_LL_plus_minus("model_LL_plus_minus","model_LL_plus_minus",
                                   RooArgList(gauss_all_LL_plus, poly_LL_plus), RooArgList(Nsig_LL_plus_minus, Nbkg_LL_plus_minus)) ;
    RooAddPdf  model_LL_minus_minus("model_LL_minus_minus","model_LL_minus_minus",
                                    RooArgList(gauss_all_LL_minus, poly_LL_minus), RooArgList(Nsig_LL_minus_minus, Nbkg_LL_minus_minus)) ;



    //---create combined datasample for simultaneous fit----
    RooCategory alldata("alldata","alldata") ;
    alldata.defineType("DD_plus_plus") ;
    alldata.defineType("DD_plus_minus") ;
    alldata.defineType("DD_minus_plus") ;
    alldata.defineType("DD_minus_minus") ;
    alldata.defineType("LL_plus_plus") ;
    alldata.defineType("LL_plus_minus") ;
    alldata.defineType("LL_minus_plus") ;
    alldata.defineType("LL_minus_minus") ;

    map<string, RooDataSet*> m= {
        {"DD_plus_plus",DD_plus_plus}, {"DD_plus_minus",DD_plus_minus},
        {"DD_minus_plus",DD_minus_plus}, {"DD_minus_minus",DD_minus_minus},
        {"LL_plus_plus",LL_plus_plus}, {"LL_plus_minus",LL_plus_minus},
        {"LL_minus_plus",LL_minus_plus}, {"LL_minus_minus",LL_minus_minus}
    };

    RooDataSet combData_all("combData", "combined data", variables, Index(alldata), Import(m));


    //here, dh1 and dh1m are the two datasets and mass is the variable in which you fit

    //---define simultaneous PDF----
    RooSimultaneous simPdf("simPdf","simultaneous pdf", alldata) ;
    simPdf.addPdf(model_DD_plus_plus,"DD_plus_plus") ;
    simPdf.addPdf(model_DD_plus_minus,"DD_plus_minus") ;
    simPdf.addPdf(model_DD_minus_plus,"DD_minus_plus") ;
    simPdf.addPdf(model_DD_minus_minus,"DD_minus_minus") ;
    simPdf.addPdf(model_LL_plus_plus,"LL_plus_plus") ;
    simPdf.addPdf(model_LL_plus_minus,"LL_plus_minus") ;
    simPdf.addPdf(model_LL_minus_plus,"LL_minus_plus") ;
    simPdf.addPdf(model_LL_minus_minus,"LL_minus_minus") ;



    cout<<"begin fit"<<endl;
    //---do the sim. fit---
    auto result= simPdf.fitTo(combData_all, RooFit::NumCPU(10), RooFit::Save(kTRUE), RooFit::Minos(kTRUE)) ;
    Int_t fit_index = result->status();
    //model_plus_plus.fitTo(*ds_plus_plus);

    cout<<"end fit"<<endl;



    //plot
    RooFormulaVar Acp_blind("Acp", "0.5*(@0-@1)", RooArgList(AT_plus_unblind, AT_minus_unblind));

    //---draw a component----
    TCanvas* c = new TCanvas("total_plot","total_plot", 1800, 2000) ;
    c->Divide(2,4) ;

    //=====================================================
    RooPlot* frame1 = D_M.frame(Title("D+_DD(Ct>0)"));
    combData_all.plotOn(frame1,Cut("alldata==alldata::DD_plus_plus")) ;
    simPdf.plotOn(frame1,Slice(alldata,"DD_plus_plus"),ProjWData(alldata, combData_all)) ;
    //simPdf.paramOn(frame1);
    combData_all.statOn(frame1);
    c->cd(1) ;
    frame1->Draw();

    RooPlot* frame2 = D_M.frame(Title("D+_DD(Ct><0)"));
    combData_all.plotOn(frame2,Cut("alldata==alldata::DD_plus_minus")) ;
    simPdf.plotOn(frame2,Slice(alldata,"DD_plus_minus"),ProjWData(alldata, combData_all)) ;
    //simPdf.paramOn(frame2);
    combData_all.statOn(frame2);
    c->cd(2) ;
    frame2->Draw();

    RooPlot* frame3 = D_M.frame(Title("D-_DD(Ct>0)"));
    combData_all.plotOn(frame3,Cut("alldata==alldata::DD_minus_plus")) ;
    simPdf.plotOn(frame3,Slice(alldata,"DD_minus_plus"),ProjWData(alldata, combData_all)) ;
    //simPdf.paramOn(frame3);
    combData_all.statOn(frame3);
    c->cd(3) ;
    frame3->Draw();

    RooPlot* frame4 = D_M.frame(Title("D-_DD(Ct<0)"));
    combData_all.plotOn(frame4,Cut("alldata==alldata::DD_minus_minus")) ;
    simPdf.plotOn(frame4,Slice(alldata,"DD_minus_minus"),ProjWData(alldata, combData_all)) ;
    //simPdf.paramOn(frame4);
    combData_all.statOn(frame4);
    c->cd(4) ;
    frame4->Draw();

    RooPlot* frame5 = D_M.frame(Title("D+_LL(Ct>0)"));
    combData_all.plotOn(frame5,Cut("alldata==alldata::LL_plus_plus")) ;
    simPdf.plotOn(frame5,Slice(alldata,"LL_plus_plus"),ProjWData(alldata, combData_all)) ;
    //simPdf.paramOn(frame5);
    combData_all.statOn(frame5);
    c->cd(5) ;
    frame5->Draw();

    RooPlot* frame6 = D_M.frame(Title("D+_LL(Ct<0)"));
    combData_all.plotOn(frame6,Cut("alldata==alldata::LL_plus_minus")) ;
    simPdf.plotOn(frame6,Slice(alldata,"LL_plus_minus"),ProjWData(alldata, combData_all)) ;
    //simPdf.paramOn(frame6);
    combData_all.statOn(frame6);
    c->cd(6) ;
    frame6->Draw();

    RooPlot* frame7 = D_M.frame(Title("D-_LL(Ct>0)"));
    combData_all.plotOn(frame7,Cut("alldata==alldata::LL_minus_plus")) ;
    simPdf.plotOn(frame7,Slice(alldata,"LL_minus_plus"),ProjWData(alldata, combData_all)) ;
    //simPdf.paramOn(frame7);
    combData_all.statOn(frame7);
    c->cd(7) ;
    frame7->Draw();

    RooPlot* frame8 = D_M.frame(Title("D-_LL(Ct<0)"));
    combData_all.plotOn(frame8,Cut("alldata==alldata::LL_minus_minus")) ;
    simPdf.plotOn(frame8,Slice(alldata,"LL_minus_minus"),ProjWData(alldata, combData_all)) ;
    //simPdf.paramOn(frame8);
    combData_all.statOn(frame8);
    c->cd(8) ;
    frame8->Draw();


    c->Print("./all/At.eps");
    c->Print("./all/At.pdf");

    Ntot_DD_plus.Print();
    Ntot_DD_minus.Print();
    Ntot_DD_bkg_plus.Print();
    Ntot_DD_bkg_minus.Print();
    Ntot_LL_plus.Print();
    Ntot_LL_minus.Print();
    Ntot_LL_bkg_plus.Print();
    Ntot_LL_bkg_minus.Print();

    AT_plus.Print();
    AT_minus.Print();
    AT_DD_bkg_plus.Print();
    AT_DD_bkg_minus.Print();
    AT_LL_bkg_plus.Print();
    AT_LL_bkg_minus.Print();
    AT_plus_unblind.Print();
    AT_minus_unblind.Print();
    //result->Print("V");




    //cout<< "Acp_unblind = " << Acp_unblind.getVal() <<" +/- "
    //	<< Acp_unblind.getPropagatedError(*result)<<endl;
    //Acp.getPropagatedError(*result);


    std::ofstream myfile;
    myfile.open ("./all/fitresult.txt");
    myfile<<"Fitted number of events(DD): "<<
          Ntot_DD_plus.getVal()+Ntot_DD_minus.getVal()+Ntot_DD_bkg_plus.getVal()+Ntot_DD_bkg_minus.getVal() <<endl;
    myfile<<"Fitted number of events(LL): "<<
          Ntot_LL_plus.getVal()+Ntot_LL_minus.getVal()+Ntot_LL_bkg_plus.getVal()+Ntot_LL_bkg_minus.getVal() <<endl;

    myfile<<"Number of total entries in dataset(DD): "<<DD_plus_plus->numEntries()+
          DD_plus_minus->numEntries()+DD_minus_plus->numEntries()+DD_minus_minus->numEntries()<<endl;
    myfile<<"Number of total entries in dataset(LL): "<<LL_plus_plus->numEntries()+
          LL_plus_minus->numEntries()+LL_minus_plus->numEntries()+LL_minus_minus->numEntries()<<endl;

    myfile<<"signal yield(DD): "<< Ntot_DD_plus.getVal()+Ntot_DD_minus.getVal()<<endl;
    myfile<<"signal yield(LL): "<< Ntot_LL_plus.getVal()+Ntot_LL_minus.getVal()<<endl;
    myfile<<"background yield(DD): "<< Ntot_DD_bkg_plus.getVal()+Ntot_DD_bkg_minus.getVal()<<endl;
    myfile<<"background yield(LL): "<< Ntot_LL_bkg_plus.getVal()+Ntot_LL_bkg_minus.getVal()<<endl;
    myfile<<"============================="<<endl;
    myfile<<"At(unblind): "<<AT_plus.getVal()<<endl;
    myfile<<"At_bar(unblind): "<<AT_minus.getVal()<<endl;
    myfile<<"At(blind): "<<AT_plus_unblind.getVal()<<endl;
    myfile<<"At_bar(blind): "<<AT_minus_unblind.getVal()<<endl;

    myfile<<"Acp_blind: "<< Acp_blind.getVal() << " +/- "
          <<Acp_blind.getPropagatedError(*result)<<endl;
    myfile<<"fit status: "<<fit_index<<endl;
    myfile.close();

}

