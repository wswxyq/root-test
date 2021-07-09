#include <iostream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

#include "TLegend.h"
#include "TCut.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooStats/SPlot.h"
using namespace RooStats ;

using namespace std;
using namespace RooFit;

int fitPipi_wSplot() {

    TChain *chain= new TChain("chain");
    chain->Add("D02pipi.root");

    int nbins=100;

    //////////////////////////////////add variable and cut
    RooRealVar mass("D0_M", "mass D^{0}", 1810, 1920, "MeV");
    RooRealVar pione_ProbNNk("pione_ProbNNk", "pione_ProbNNk", 0, 1);
    RooRealVar pitwo_ProbNNk("pitwo_ProbNNk", "pitwo_ProbNNk", 0, 1);
    RooRealVar pione_ProbNNpi("pione_ProbNNpi", "pione_ProbNNpi", 0, 1);
    RooRealVar pitwo_isMuon("pitwo_isMuon", "pitwo_isMuon",0,1);
    RooRealVar pitwo_ProbNNmu("pitwo_ProbNNmu", "pitwo_ProbNNmu", 0, 1);


    RooArgSet variables(mass, pione_ProbNNk,pitwo_ProbNNk, pione_ProbNNpi);
    variables.add(RooArgSet(pitwo_isMuon, pitwo_ProbNNmu));
    //////////////////////////////////
    //Peak component
    RooRealVar mean1("mean1","mean1",1860.5);
    RooRealVar sigma1("sigma1","sigma1",23.5) ;
    RooRealVar a1("a1", "a1", 1.61);
    RooRealVar n1("n1", "n1", 0.33);
    RooCBShape gauss1("gauss1", "gauss1", mass, mean1, sigma1, a1, n1);

    RooRealVar mean2("mean2","mean2",1865.8);
    RooRealVar sigma2("sigma2","sigma2",8.33) ;
    RooGaussian gauss2("gauss2", "gauss2", mass, mean2, sigma2);

    RooRealVar frac1("frac1", "frac1", 0.157);
    RooAddPdf peak("peak", "peak", RooArgSet(gauss1, gauss2),RooArgSet(frac1));

    //////////////////////////////////
    //Combinatoric component
    RooRealVar c1("c1", "c1", -0.1114); //-0.05, -1, 1
    RooRealVar c2("c2", "c2", 0.039); //-0.07, -1, 1);
    RooChebychev poly1("poly1", "poly1", mass, RooArgSet(c1, c2));

    //////////////////////////////////
    //Kpi component
    RooRealVar mean3("mean3","mean3",1700);
    RooRealVar sigma3("sigma3","sigma3",28.3) ;
    RooGaussian Kpi_PDF("Kpi_PDF", "Kpi_PDF", mass, mean3, sigma3);

    /////////////////////////////////
    //Total PDF
    RooRealVar Npipi("Npipi", "Npipi", 1500, 3000);
    RooRealVar Nkomb("Nkomb", "Nkomb", 7000, 11000);
    RooRealVar Nkpi("Nkpi", "Nkpi", 40, 1000);
    RooAddPdf PDF("PDF", "PDF", RooArgSet(peak, poly1, Kpi_PDF), RooArgSet(Npipi, Nkomb, Nkpi));


    //////////////////////////////////
    RooDataSet *data=new RooDataSet("data", "data", variables, Import(*chain), Cut("pione_ProbNNpi>0.85 && pione_ProbNNk<0.05 && pitwo_ProbNNk<0.3 && pitwo_isMuon==1"));
    cout<<"entries dataset: "<<data->numEntries()<<endl;

    PDF.fitTo(*data);

    RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, &PDF, RooArgList(Npipi,Nkomb,Nkpi) );

    std::cout << "Check SWeights:" << std::endl;
    std::cout << std::endl <<  "Yield of pipi is "
              << Npipi.getVal() << ".  From sWeights it is "
              << sData->GetYieldFromSWeight("Npipi") << std::endl;
    std::cout << std::endl <<  "Yield of kpi is "
              << Nkpi.getVal() << ".  From sWeights it is "
              << sData->GetYieldFromSWeight("Nkpi") << std::endl;
    std::cout << std::endl <<  "Yield of komb is "
              << Nkomb.getVal() << ".  From sWeights it is "
              << sData->GetYieldFromSWeight("Nkomb") << std::endl;



    TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
    cdata->Divide(1,3);

    // create weightfed data set
    RooDataSet * dataw_z = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"Npipi_sw") ;


    RooPlot *frame = pitwo_ProbNNmu.frame( Title(" "));
    dataw_z->plotOn(frame, DataError(RooAbsData::SumW2) ) ;
    cdata->cd(1);
    frame->Draw();

    RooDataSet * dataw_komb = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"Nkomb_sw") ;
    RooPlot* frame3 = pitwo_ProbNNmu.frame() ;
    dataw_komb->plotOn(frame3,DataError(RooAbsData::SumW2) ) ;
    cdata->cd(2);
    frame3->Draw();

    RooDataSet * dataw_kpi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"Nkpi_sw") ;
    RooPlot* frame4 = pitwo_ProbNNmu.frame() ;
    dataw_kpi->plotOn(frame4,DataError(RooAbsData::SumW2) ) ;
    cdata->cd(3);
    frame4->Draw();

    cdata->Print("/eos/user/t/tnanut/D02mumu/plots/sPlot_ProbNNmu.pdf");

    return 0;
}
