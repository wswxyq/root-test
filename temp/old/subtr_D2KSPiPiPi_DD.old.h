#include <iostream> 
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"

using namespace std;

int subtr_D2KSPiPiPi_DD(){

///////////////////////////////////

	//Load a file 

//////////////////////////////////
	Double_t ymax=0;
	TChain *chain= new TChain("D2KSPiPiPi_DD/DecayTree");
	chain->Add("/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root");
	//you can add as many files as needed
	cout<<"Number of entries: "<< chain->GetEntries()<<endl; 
	//prints in the command line the number of entries

    TCut trigger = "( (D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
(D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1) ) && ( (D_Hlt1TrackMVADecision_TOS==1) || \
(D_Hlt1TwoTrackMVADecision_TOS==1) )";
    TCut SignalWindow = "(D_M > 1839.65) && (D_M < 1899.65)";
    TCut SideBand = "( (D_M > 1779.65) && (D_M < 1809.65) ) || ( (D_M > 1929.65) && (D_M < 1959.65) )";
//////////////////////////////////
	//Step 1: Draw a histogram of the D0 mass and save it to a pdf

///////////D_M

	TH1F *h1_D_M = new TH1F("h1_D_M", "D_M", 100, 1700, 2100);
	TH1F *h2_D_M = new TH1F("h2_D_M", "D_M", 100, 1700, 2100);
	TH1F *h3_D_M = new TH1F("h3_D_M", "D_M", 100, 1700, 2100);
	chain->Project("h1_D_M", "D_M", trigger) ;
	chain->Project("h2_D_M", "D_M", trigger&&SignalWindow) ;
	chain->Project("h3_D_M", "D_M", trigger&&SideBand) ;

	//this fills the histogram h1 with variable D0_M
	TCanvas *c1_D_M = new TCanvas( "c1", "c1", 200, 10, 700, 500 );
	ymax=h1_D_M->GetMaximum();
	TLine *line1 =new TLine(1839.65,0,1839.65,ymax);
	TLine *line2 =new TLine(1899.65,0,1899.65,ymax);


	h2_D_M->SetFillColor(8);
	h2_D_M->SetFillStyle(3005);
	h3_D_M->SetFillColor(4);
	h3_D_M->SetFillStyle(3016);

	h1_D_M->Draw();
	h2_D_M->Draw("same");
	h3_D_M->Draw("same");
	line1->Draw("same");
	line2->Draw("same");
	c1_D_M->Print("h_M_D.pdf");


///////D_PT

	//TH1F *h1_D_PT = new TH1F("h1_D_PT", "D_PT", 100, 0, 20000);
	TH1F *h2_D_PT = new TH1F("h2_D_PT", "D_PT", 100, 0, 20000);
	TH1F *h3_D_PT = new TH1F("h3_D_PT", "D_PT", 100, 0, 20000);

	chain->Project("h1_D_PT", "D_PT", trigger&&SignalWindow) ;
	chain->Project("h2_D_PT", "D_PT", trigger&&SideBand) ;
	chain->Project("h3_D_PT", "D_PT", trigger&&SignalWindow) ;
	//this fills the histogram h*_D_PT with variable D_PT

	TCanvas *c1_D_PT = new TCanvas( "c1_D_PT", "c1_D_PT", 200, 10, 700, 500 );
	h2_D_PT->SetFillColor(8);
	h2_D_PT->SetFillStyle(3005);
	h2_D_PT->Draw();
	c1_D_PT->Print("h_D_PT_bkg.pdf");

	TCanvas *c2_D_PT = new TCanvas( "c2_D_PT", "c2_D_PT", 200, 10, 700, 500 );
	h3_D_PT->Add(h2_D_PT, -1);
	h3_D_PT->SetFillColor(4);
	h3_D_PT->SetFillStyle(3016);
	h3_D_PT->Draw();
	c2_D_PT->Print("h_D_PT_sgnl.pdf");

	TCanvas *c3_D_PT = new TCanvas( "c3_D_PT", "c3_D_PT", 200, 10, 700, 500 );
	h2_D_PT->DrawNormalized();
	h3_D_PT->DrawNormalized("same");
	c3_D_PT->Print("h_D_PT_bkg&sgnl.pdf");




///////D_ReFit_P

	//TH1F *h1_D_ReFit_P = new TH1F("h1_D_ReFit_P", "D_ReFit_P", 100, 0, 500000);
	TH1F *h2_D_ReFit_P = new TH1F("h2_D_ReFit_P", "D_ReFit_P", 100, 0, 300000);
	TH1F *h3_D_ReFit_P = new TH1F("h3_D_ReFit_P", "D_ReFit_P", 100, 0, 300000);

	chain->Project("h1_D_ReFit_P", "D_ReFit_P", trigger&&SignalWindow) ;
	chain->Project("h2_D_ReFit_P", "D_ReFit_P", trigger&&SideBand) ;
	chain->Project("h3_D_ReFit_P", "D_ReFit_P", trigger&&SignalWindow) ;
	//this fills the histogram h*_D_ReFit_P with variable D_ReFit_P

	TCanvas *c1_D_ReFit_P = new TCanvas( "c1_D_ReFit_P", "c1_D_ReFit_P", 200, 10, 700, 500 );
	h2_D_ReFit_P->SetFillColor(8);
	h2_D_ReFit_P->SetFillStyle(3005);
	h2_D_ReFit_P->Draw();
	c1_D_ReFit_P->Print("h_D_ReFit_P_bkg.pdf");

	TCanvas *c2_D_ReFit_P = new TCanvas( "c2_D_ReFit_P", "c2_D_ReFit_P", 200, 10, 700, 500 );
	h3_D_ReFit_P->Add(h2_D_ReFit_P, -1);
	h3_D_ReFit_P->SetFillColor(4);
	h3_D_ReFit_P->SetFillStyle(3016);
	h3_D_ReFit_P->Draw();
	c2_D_ReFit_P->Print("h_D_ReFit_P_sgnl.pdf");

	TCanvas *c3_D_ReFit_P = new TCanvas( "c3_D_ReFit_P", "c3_D_ReFit_P", 200, 10, 700, 500 );
	h2_D_ReFit_P->DrawNormalized();
	h3_D_ReFit_P->DrawNormalized("same");
	c3_D_ReFit_P->Print("h_D_ReFit_P_bkg&sgnl.pdf");



///////D_ReFit_chi2

	//TH1F *h1_D_ReFit_chi2 = new TH1F("h1_D_ReFit_chi2", "D_ReFit_chi2", 100, 0, 500000);
	TH1F *h2_D_ReFit_chi2 = new TH1F("h2_D_ReFit_chi2", "D_ReFit_chi2", 100, 0, 100);
	TH1F *h3_D_ReFit_chi2 = new TH1F("h3_D_ReFit_chi2", "D_ReFit_chi2", 100, 0, 100);

	chain->Project("h1_D_ReFit_chi2", "D_ReFit_chi2", trigger&&SignalWindow) ;
	chain->Project("h2_D_ReFit_chi2", "D_ReFit_chi2", trigger&&SideBand) ;
	chain->Project("h3_D_ReFit_chi2", "D_ReFit_chi2", trigger&&SignalWindow) ;
	//this fills the histogram h*_D_ReFit_chi2 with variable D_ReFit_chi2

	TCanvas *c1_D_ReFit_chi2 = new TCanvas( "c1_D_ReFit_chi2", "c1_D_ReFit_chi2", 200, 10, 700, 500 );
	h2_D_ReFit_chi2->SetFillColor(8);
	h2_D_ReFit_chi2->SetFillStyle(3005);
	h2_D_ReFit_chi2->Draw();
	c1_D_ReFit_chi2->Print("h_D_ReFit_chi2_bkg.pdf");

	TCanvas *c2_D_ReFit_chi2 = new TCanvas( "c2_D_ReFit_chi2", "c2_D_ReFit_chi2", 200, 10, 700, 500 );
	h3_D_ReFit_chi2->Add(h2_D_ReFit_chi2, -1);
	h3_D_ReFit_chi2->SetFillColor(4);
	h3_D_ReFit_chi2->SetFillStyle(3016);
	h3_D_ReFit_chi2->Draw();
	c2_D_ReFit_chi2->Print("h_D_ReFit_chi2_sgnl.pdf");

	TCanvas *c3_D_ReFit_chi2 = new TCanvas( "c3_D_ReFit_chi2", "c3_D_ReFit_chi2", 200, 10, 700, 500 );
	h2_D_ReFit_chi2->DrawNormalized();
	h3_D_ReFit_chi2->DrawNormalized("same");
	c3_D_ReFit_chi2->Print("h_D_ReFit_chi2_bkg&sgnl.pdf");




///////KS_PT

	//TH1F *h1_KS_PT = new TH1F("h1_KS_PT", "KS_PT", 100, 0, 500000);
	TH1F *h2_KS_PT = new TH1F("h2_KS_PT", "KS_PT", 100, 400, 600);
	TH1F *h3_KS_PT = new TH1F("h3_KS_PT", "KS_PT", 100, 400, 600);

	chain->Project("h1_KS_PT", "KS_PT", trigger&&SignalWindow) ;
	chain->Project("h2_KS_PT", "KS_PT", trigger&&SideBand) ;
	chain->Project("h3_KS_PT", "KS_PT", trigger&&SignalWindow) ;
	//this fills the histogram h*_KS_PT with variable KS_PT

	TCanvas *c1_KS_PT = new TCanvas( "c1_KS_PT", "c1_KS_PT", 200, 10, 700, 500 );
	h2_KS_PT->SetFillColor(8);
	h2_KS_PT->SetFillStyle(3005);
	h2_KS_PT->Draw();
	c1_KS_PT->Print("h_KS_PT_bkg.pdf");

	TCanvas *c2_KS_PT = new TCanvas( "c2_KS_PT", "c2_KS_PT", 200, 10, 700, 500 );
	h3_KS_PT->Add(h2_KS_PT, -1);
	h3_KS_PT->SetFillColor(4);
	h3_KS_PT->SetFillStyle(3016);
	h3_KS_PT->Draw();
	c2_KS_PT->Print("h_KS_PT_sgnl.pdf");

	TCanvas *c3_KS_PT = new TCanvas( "c3_KS_PT", "c3_KS_PT", 200, 10, 700, 500 );
	h2_KS_PT->DrawNormalized();
	h3_KS_PT->DrawNormalized("same");
	c3_KS_PT->Print("h_KS_PT_bkg&sgnl.pdf");




///////D_ReFit_KS0_P

	//TH1F *h1_D_ReFit_KS0_P = new TH1F("h1_D_ReFit_KS0_P", "D_ReFit_KS0_P", 100, 0, 500000);
	TH1F *h2_D_ReFit_KS0_P = new TH1F("h2_D_ReFit_KS0_P", "D_ReFit_KS0_P", 100, 0, 100000);
	TH1F *h3_D_ReFit_KS0_P = new TH1F("h3_D_ReFit_KS0_P", "D_ReFit_KS0_P", 100, 0, 100000);

	chain->Project("h1_D_ReFit_KS0_P", "D_ReFit_KS0_P", trigger&&SignalWindow) ;
	chain->Project("h2_D_ReFit_KS0_P", "D_ReFit_KS0_P", trigger&&SideBand) ;
	chain->Project("h3_D_ReFit_KS0_P", "D_ReFit_KS0_P", trigger&&SignalWindow) ;
	//this fills the histogram h*_D_ReFit_KS0_P with variable D_ReFit_KS0_P

	TCanvas *c1_D_ReFit_KS0_P = new TCanvas( "c1_D_ReFit_KS0_P", "c1_D_ReFit_KS0_P", 200, 10, 700, 500 );
	h2_D_ReFit_KS0_P->SetFillColor(8);
	h2_D_ReFit_KS0_P->SetFillStyle(3005);
	h2_D_ReFit_KS0_P->Draw();
	c1_D_ReFit_KS0_P->Print("h_D_ReFit_KS0_P_bkg.pdf");

	TCanvas *c2_D_ReFit_KS0_P = new TCanvas( "c2_D_ReFit_KS0_P", "c2_D_ReFit_KS0_P", 200, 10, 700, 500 );
	h3_D_ReFit_KS0_P->Add(h2_D_ReFit_KS0_P, -1);
	h3_D_ReFit_KS0_P->SetFillColor(4);
	h3_D_ReFit_KS0_P->SetFillStyle(3016);
	h3_D_ReFit_KS0_P->Draw();
	c2_D_ReFit_KS0_P->Print("h_D_ReFit_KS0_P_sgnl.pdf");

	TCanvas *c3_D_ReFit_KS0_P = new TCanvas( "c3_D_ReFit_KS0_P", "c3_D_ReFit_KS0_P", 200, 10, 700, 500 );
	h2_D_ReFit_KS0_P->DrawNormalized();
	h3_D_ReFit_KS0_P->DrawNormalized("same");
	c3_D_ReFit_KS0_P->Print("h_D_ReFit_KS0_P_bkg&sgnl.pdf");




///////P1_ProbNNpi

	//TH1F *h1_P1_ProbNNpi = new TH1F("h1_P1_ProbNNpi", "P1_ProbNNpi", 100, 0, 500000);
	TH1F *h2_P1_ProbNNpi = new TH1F("h2_P1_ProbNNpi", "P1_ProbNNpi", 100, 0.8, 1.01);
	TH1F *h3_P1_ProbNNpi = new TH1F("h3_P1_ProbNNpi", "P1_ProbNNpi", 100, 0.8, 1.01);

	chain->Project("h1_P1_ProbNNpi", "P1_ProbNNpi", trigger&&SignalWindow) ;
	chain->Project("h2_P1_ProbNNpi", "P1_ProbNNpi", trigger&&SideBand) ;
	chain->Project("h3_P1_ProbNNpi", "P1_ProbNNpi", trigger&&SignalWindow) ;
	//this fills the histogram h*_P1_ProbNNpi with variable P1_ProbNNpi

	TCanvas *c1_P1_ProbNNpi = new TCanvas( "c1_P1_ProbNNpi", "c1_P1_ProbNNpi", 200, 10, 700, 500 );
	h2_P1_ProbNNpi->SetFillColor(8);
	h2_P1_ProbNNpi->SetFillStyle(3005);
	h2_P1_ProbNNpi->Draw();
	c1_P1_ProbNNpi->Print("h_P1_ProbNNpi_bkg.pdf");

	TCanvas *c2_P1_ProbNNpi = new TCanvas( "c2_P1_ProbNNpi", "c2_P1_ProbNNpi", 200, 10, 700, 500 );
	h3_P1_ProbNNpi->Add(h2_P1_ProbNNpi, -1);
	h3_P1_ProbNNpi->SetFillColor(4);
	h3_P1_ProbNNpi->SetFillStyle(3016);
	h3_P1_ProbNNpi->Draw();
	c2_P1_ProbNNpi->Print("h_P1_ProbNNpi_sgnl.pdf");

	TCanvas *c3_P1_ProbNNpi = new TCanvas( "c3_P1_ProbNNpi", "c3_P1_ProbNNpi", 200, 10, 700, 500 );
	h2_P1_ProbNNpi->DrawNormalized();
	h3_P1_ProbNNpi->DrawNormalized("same");
	c3_P1_ProbNNpi->Print("h_P1_ProbNNpi_bkg&sgnl.pdf");



///////P2_ProbNNpi

	//TH1F *h1_P2_ProbNNpi = new TH1F("h1_P2_ProbNNpi", "P2_ProbNNpi", 100, 0, 500000);
	TH1F *h2_P2_ProbNNpi = new TH1F("h2_P2_ProbNNpi", "P2_ProbNNpi", 100, 0.8, 1.01);
	TH1F *h3_P2_ProbNNpi = new TH1F("h3_P2_ProbNNpi", "P2_ProbNNpi", 100, 0.8, 1.01);

	chain->Project("h1_P2_ProbNNpi", "P2_ProbNNpi", trigger&&SignalWindow) ;
	chain->Project("h2_P2_ProbNNpi", "P2_ProbNNpi", trigger&&SideBand) ;
	chain->Project("h3_P2_ProbNNpi", "P2_ProbNNpi", trigger&&SignalWindow) ;
	//this fills the histogram h*_P2_ProbNNpi with variable P2_ProbNNpi

	TCanvas *c1_P2_ProbNNpi = new TCanvas( "c1_P2_ProbNNpi", "c1_P2_ProbNNpi", 200, 10, 700, 500 );
	h2_P2_ProbNNpi->SetFillColor(8);
	h2_P2_ProbNNpi->SetFillStyle(3005);
	h2_P2_ProbNNpi->Draw();
	c1_P2_ProbNNpi->Print("h_P2_ProbNNpi_bkg.pdf");

	TCanvas *c2_P2_ProbNNpi = new TCanvas( "c2_P2_ProbNNpi", "c2_P2_ProbNNpi", 200, 10, 700, 500 );
	h3_P2_ProbNNpi->Add(h2_P2_ProbNNpi, -1);
	h3_P2_ProbNNpi->SetFillColor(4);
	h3_P2_ProbNNpi->SetFillStyle(3016);
	h3_P2_ProbNNpi->Draw();
	c2_P2_ProbNNpi->Print("h_P2_ProbNNpi_sgnl.pdf");

	TCanvas *c3_P2_ProbNNpi = new TCanvas( "c3_P2_ProbNNpi", "c3_P2_ProbNNpi", 200, 10, 700, 500 );
	h2_P2_ProbNNpi->DrawNormalized();
	h3_P2_ProbNNpi->DrawNormalized("same");
	c3_P2_ProbNNpi->Print("h_P2_ProbNNpi_bkg&sgnl.pdf");


///////P0_ProbNNpi

	//TH1F *h1_P0_ProbNNpi = new TH1F("h1_P0_ProbNNpi", "P0_ProbNNpi", 100, 0, 500000);
	TH1F *h2_P0_ProbNNpi = new TH1F("h2_P0_ProbNNpi", "P0_ProbNNpi", 100, 0.8, 1.01);
	TH1F *h3_P0_ProbNNpi = new TH1F("h3_P0_ProbNNpi", "P0_ProbNNpi", 100, 0.8, 1.01);

	chain->Project("h1_P0_ProbNNpi", "P0_ProbNNpi", trigger&&SignalWindow) ;
	chain->Project("h2_P0_ProbNNpi", "P0_ProbNNpi", trigger&&SideBand) ;
	chain->Project("h3_P0_ProbNNpi", "P0_ProbNNpi", trigger&&SignalWindow) ;
	//this fills the histogram h*_P0_ProbNNpi with variable P0_ProbNNpi

	TCanvas *c1_P0_ProbNNpi = new TCanvas( "c1_P0_ProbNNpi", "c1_P0_ProbNNpi", 200, 10, 700, 500 );
	h2_P0_ProbNNpi->SetFillColor(8);
	h2_P0_ProbNNpi->SetFillStyle(3005);
	h2_P0_ProbNNpi->Draw();
	c1_P0_ProbNNpi->Print("h_P0_ProbNNpi_bkg.pdf");

	TCanvas *c2_P0_ProbNNpi = new TCanvas( "c2_P0_ProbNNpi", "c2_P0_ProbNNpi", 200, 10, 700, 500 );
	h3_P0_ProbNNpi->Add(h2_P0_ProbNNpi, -1);
	h3_P0_ProbNNpi->SetFillColor(4);
	h3_P0_ProbNNpi->SetFillStyle(3016);
	h3_P0_ProbNNpi->Draw();
	c2_P0_ProbNNpi->Print("h_P0_ProbNNpi_sgnl.pdf");

	TCanvas *c3_P0_ProbNNpi = new TCanvas( "c3_P0_ProbNNpi", "c3_P0_ProbNNpi", 200, 10, 700, 500 );
	h2_P0_ProbNNpi->DrawNormalized();
	h3_P0_ProbNNpi->DrawNormalized("same");
	c3_P0_ProbNNpi->Print("h_P0_ProbNNpi_bkg&sgnl.pdf");



///////D_ReFit_decayLength/D_ReFit_decayLengthErr(D_ReFit_decayLE)

	//TH1F *h1_D_ReFit_decayLE = new TH1F("h1_D_ReFit_decayLE", "D_ReFit_decayLE", 100, 0, 500000);
	TH1F *h2_D_ReFit_decayLE = new TH1F("h2_D_ReFit_decayLE", "D_ReFit_decayLE", 100, 0, 300);
	TH1F *h3_D_ReFit_decayLE = new TH1F("h3_D_ReFit_decayLE", "D_ReFit_decayLE", 100, 0, 300);

	chain->Project("h1_D_ReFit_decayLE", "D_ReFit_decayLength/D_ReFit_decayLengthErr", trigger&&SignalWindow) ;
	chain->Project("h2_D_ReFit_decayLE", "D_ReFit_decayLength/D_ReFit_decayLengthErr", trigger&&SideBand) ;
	chain->Project("h3_D_ReFit_decayLE", "D_ReFit_decayLength/D_ReFit_decayLengthErr", trigger&&SignalWindow) ;
	//this fills the histogram h*_D_ReFit_decayLE with variable D_ReFit_decayLE

	TCanvas *c1_D_ReFit_decayLE = new TCanvas( "c1_D_ReFit_decayLE", "c1_D_ReFit_decayLE", 200, 10, 700, 500 );
	h2_D_ReFit_decayLE->SetFillColor(8);
	h2_D_ReFit_decayLE->SetFillStyle(3005);
	h2_D_ReFit_decayLE->Draw();
	c1_D_ReFit_decayLE->Print("h_D_ReFit_decayLE_bkg.pdf");

	TCanvas *c2_D_ReFit_decayLE = new TCanvas( "c2_D_ReFit_decayLE", "c2_D_ReFit_decayLE", 200, 10, 700, 500 );
	h3_D_ReFit_decayLE->Add(h2_D_ReFit_decayLE, -1);
	h3_D_ReFit_decayLE->SetFillColor(4);
	h3_D_ReFit_decayLE->SetFillStyle(3016);
	h3_D_ReFit_decayLE->Draw();
	c2_D_ReFit_decayLE->Print("h_D_ReFit_decayLE_sgnl.pdf");

	TCanvas *c3_D_ReFit_decayLE = new TCanvas( "c3_D_ReFit_decayLE", "c3_D_ReFit_decayLE", 200, 10, 700, 500 );
	h2_D_ReFit_decayLE->DrawNormalized();
	h3_D_ReFit_decayLE->DrawNormalized("same");
	c3_D_ReFit_decayLE->Print("h_D_ReFit_decayLE_bkg&sgnl.pdf");


return 0;
}