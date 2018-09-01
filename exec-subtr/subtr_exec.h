#include <iostream> 
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "subtr_common.h"

using namespace std;

int subtr_D2KSPiPiPi_DD(){

	///////////////////////////////////
    //subtr_D2KSPiPiPi_DD();
    const char* TChain_name="D2KSPiPiPi_DD/DecayTree";
    const char* TChain_file="/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root";
    const char* trigger_content="( (D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
(D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1) ) && ( (D_Hlt1TrackMVADecision_TOS==1) || \
(D_Hlt1TwoTrackMVADecision_TOS==1) )";
    const char* SignalWindow_content="(D_M > 1839.65) && (D_M < 1899.65)";
    const char* SideBand_content="( (D_M > 1779.65) && (D_M < 1809.65) ) || ( (D_M > 1929.65) && (D_M < 1959.65) )";


	//////////////////////////////////
	Double_t ymax=0;
	TChain *chain= new TChain(TChain_name);
	chain->Add(TChain_file);
	//you can add as many files as needed
	cout<<"Number of entries: "<< chain->GetEntries()<<endl; 
	//prints in the command line the number of entries

    TCut trigger = trigger_content;
    TCut SignalWindow = SignalWindow_content;
    TCut SideBand = SideBand_content;
	//----------D_M--------------------------------------------//

	TH1F *h1_D_M = new TH1F("h1_D_M", "D_M", 100, 1700, 2100);
	TH1F *h2_D_M = new TH1F("h2_D_M", "D_M", 100, 1700, 2100);
	TH1F *h3_D_M = new TH1F("h3_D_M", "D_M", 100, 1700, 2100);
	chain->Project("h1_D_M", "D_M", trigger) ;
	chain->Project("h2_D_M", "D_M", trigger&&SignalWindow) ;
	chain->Project("h3_D_M", "D_M", trigger&&SideBand) ;

	//this fills the histogram h1 with variable D_M
	TCanvas *c1_D_M = new TCanvas( "c1", "c1", 200, 10, 700, 500 );
	ymax=h1_D_M->GetMaximum();
	TLine *line1 =new TLine(1839.65,0,1839.65,ymax);
	TLine *line2 =new TLine(1899.65,0,1899.65,ymax);
	//set fill content
	h2_D_M->SetFillColor(8);
	h2_D_M->SetFillStyle(3005);
	h3_D_M->SetFillColor(4);
	h3_D_M->SetFillStyle(3016);
	//draw and print
	h1_D_M->Draw();
	h2_D_M->Draw("same");
	h3_D_M->Draw("same");
	line1->Draw("same");
	line2->Draw("same");
	c1_D_M->Print("../plotdir/D2KSPiPiPi_DD/h_M_D.pdf");
	delete chain;
    delete h2_D_M;
    delete h3_D_M;
    delete c1_D_M;



	/////// D_PT 
    subtr_common("D_PT", "D_PT", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 20000, 0);

	/////// D_ReFit_P
    subtr_common("D_ReFit_P", "D_ReFit_P", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 300000, 1);

	///////D_ReFit_chi2
    subtr_common("D_ReFit_chi2", "D_ReFit_chi2", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 100, 1);

	///////KS_PT
    subtr_common("KS_PT", "KS_PT", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 10000, 0);

	///////D_ReFit_KS0_P
    subtr_common("D_ReFit_KS0_P", "D_ReFit_KS0_P", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 100000, 1);

	///////P1_ProbNNpi
    subtr_common("P1_ProbNNpi", "P1_ProbNNpi", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0.8, 1.01, 0);

	///////P2_ProbNNpi
    subtr_common("P2_ProbNNpi", "P2_ProbNNpi", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0.8, 1.01, 0);

	///////P0_ProbNNpi
    subtr_common("P0_ProbNNpi", "P0_ProbNNpi", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0.8, 1.01, 0);

	///////D_ReFit_decayLength/D_ReFit_decayLengthErr(D_ReFit_decayLE)
    subtr_common("D_ReFit_decayLength/D_ReFit_decayLengthErr", "D_ReFit_decayLength_vs_Err", "D2KSPiPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 300, 0);


return 0;
}



int subtr_D2KSKpPiPi_DD(){

	///////////////////////////////////
    //subtr_D2KSKpPiPi_DD();
    const char* TChain_name="D2KSKpPiPi_DD/DecayTree";
    const char* TChain_file="/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root";
    const char* trigger_content="( (D_L0HadronDecision_TOS==1) || (D_L0MuonDecision_TIS==1) || (D_L0DiMuonDecision_TIS==1) || \
(D_L0ElectronDecision_TIS==1) || (D_L0PhotonDecision_TIS==1) ) && ( (D_Hlt1TrackMVADecision_TOS==1) || \
(D_Hlt1TwoTrackMVADecision_TOS==1) )";
    const char* SignalWindow_content="(D_M > 1839.65) && (D_M < 1899.65)";
    const char* SideBand_content="( (D_M > 1809.65) && (D_M < 1839.65) ) || ( (D_M > 1899.65) && (D_M < 1929.65) )";


	//////////////////////////////////
	Double_t ymax=0;
	TChain *chain= new TChain(TChain_name);
	chain->Add(TChain_file);
	//you can add as many files as needed
	cout<<"Number of entries: "<< chain->GetEntries()<<endl; 
	//prints in the command line the number of entries

    TCut trigger = trigger_content;
    TCut SignalWindow = SignalWindow_content;
    TCut SideBand = SideBand_content;
	//----------D_M--------------------------------------------//

	TH1F *h1_D_M = new TH1F("h1_D_M", "D_M", 100, 1700, 2100);
	TH1F *h2_D_M = new TH1F("h2_D_M", "D_M", 100, 1700, 2100);
	TH1F *h3_D_M = new TH1F("h3_D_M", "D_M", 100, 1700, 2100);
	chain->Project("h1_D_M", "D_M", trigger) ;
	chain->Project("h2_D_M", "D_M", trigger&&SignalWindow) ;
	chain->Project("h3_D_M", "D_M", trigger&&SideBand) ;

	//this fills the histogram h1 with variable D_M
	TCanvas *c1_D_M = new TCanvas( "c1", "c1", 200, 10, 700, 500 );
	ymax=h1_D_M->GetMaximum();
	TLine *line1 =new TLine(1839.65,0,1839.65,ymax);
	TLine *line2 =new TLine(1899.65,0,1899.65,ymax);
	//set fill content
	h2_D_M->SetFillColor(8);
	h2_D_M->SetFillStyle(3005);
	h3_D_M->SetFillColor(4);
	h3_D_M->SetFillStyle(3016);
	//draw and print
	h1_D_M->Draw();
	h2_D_M->Draw("same");
	h3_D_M->Draw("same");
	line1->Draw("same");
	line2->Draw("same");
	c1_D_M->Print("../plotdir/D2KSKpPiPi_DD/h_M_D.pdf");
	delete chain;
    delete h2_D_M;
    delete h3_D_M;
    delete c1_D_M;



	/////// D_PT 
    subtr_common("D_PT", "D_PT", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 20000, 1);

	/////// D_ReFit_P
    subtr_common("D_ReFit_P", "D_ReFit_P", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 300000, 1);

	///////D_ReFit_chi2
    subtr_common("D_ReFit_chi2", "D_ReFit_chi2", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 100, 1);

	///////KS_PT
    subtr_common("KS_PT", "KS_PT", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 10000, 1);

	///////D_ReFit_KS0_P
    subtr_common("D_ReFit_KS0_P", "D_ReFit_KS0_P", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 100000, 1);

	///////P1_ProbNNpi
    subtr_common("P1_ProbNNpi", "P1_ProbNNpi", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0.8, 1.01, 1);

	///////P2_ProbNNpi
    subtr_common("P2_ProbNNpi", "P2_ProbNNpi", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0.8, 1.01, 1);

	///////P0_ProbNNk
    subtr_common("P0_ProbNNk", "P0_ProbNNk", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0.8, 1.01, 1);

	///////D_ReFit_decayLength/D_ReFit_decayLengthErr(D_ReFit_decayLE)
    subtr_common("D_ReFit_decayLength/D_ReFit_decayLengthErr", "D_ReFit_decayLength_vs_Err", "D2KSKpPiPi_DD", TChain_name, TChain_file, trigger_content,
SignalWindow_content, SideBand_content, 0, 70, 1);


return 0;
}