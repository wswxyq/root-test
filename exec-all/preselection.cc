//This script generate the Reduced tree after cut

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooDataSet.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"

using namespace RooFit ;
using namespace std;
using namespace TMath;

const char* eosdir="/eos/user/s/shaowei/";

int all_D2KSKpPiPi_DD()
{

	//define TChain, TFile 
	const char* TChain_name="D2KSKpPiPi_DD/DecayTree";
	const char* TChain_file1="/eos/lhcb/user/m/mamartin/D2KsKPiPi/2016/MU/*.root";
	const char* TChain_file2="/eos/lhcb/user/m/mamartin/D2KsKPiPi/2016/MD/*.root";

	TChain *chainx= new TChain(TChain_name);
	chainx->Add(TChain_file1);
	chainx->Add(TChain_file2);
    TFile *newfile = new TFile("./DD/D2KSKpPiPi_DD.root","recreate");


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
        (D_Hlt1TwoTrackMVADecision_TOS==1) )   )&& (D_M>1799.65&&D_M<1939.65) &&(D_IPCHI2_OWNPV>0)");

	newfile->cd();
    newtree->Write("ReducedTree");
	return 0;
}




int all_D2KSKpPiPi_LL()
{


	const char* TChain_name="D2KSKpPiPi_LL/DecayTree";
	const char* TChain_file1="/eos/lhcb/user/m/mamartin/D2KsKPiPi/2016/MU/*.root";
	const char* TChain_file2="/eos/lhcb/user/m/mamartin/D2KsKPiPi/2016/MD/*.root";

	TChain *chainx= new TChain(TChain_name);
	chainx->Add(TChain_file1);
	chainx->Add(TChain_file2);
    TFile *newfile = new TFile("./LL/D2KSKpPiPi_LL.root","recreate");

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
        (D_Hlt1TwoTrackMVADecision_TOS==1) )   )&& (D_M>1799.65&&D_M<1939.65) &&(D_IPCHI2_OWNPV>0)&&( (KS_PT>0.)||(KS_PT<=0.) )\
        &&( (D_PT>0.)||(D_PT<=0.) )\
        &&( (D_ReFit_P>0.)||(D_ReFit_P<=0.) )\
        &&( (D_ReFit_KS0_P>0.)||(D_ReFit_KS0_P<=0.) )\
        &&( (P0_ProbNNpi>0.)||(P0_ProbNNpi<=0.) )\
        &&( (P1_ProbNNpi>0.)||(P1_ProbNNpi<=0.) )\
        &&( (P2_ProbNNpi>0.)||(P2_ProbNNpi<=0.) )\
        &&( (D_ReFit_decayLength>0.)||(D_ReFit_decayLength<=0.) )\
        &&( (D_ReFit_decayLengthErr>0.)||(D_ReFit_decayLengthErr<0.) )\
        ");
	newfile->cd();
    newtree->Write("ReducedTree");

	return 0;
}





int preselection()
{
    /* code */
    int i1;
    int i2;
	int i3;
	int i4;
    cout << "Cut D2KSKpPiPi_DD? ('1' for yes, others for no)"<<"\n";
    cin >> i1;
    cout << "Cut D2KSKpPiPi_LL? ('1' for yes, others for no)"<<"\n";
    cin >> i2;
	/*cout << "Fit subtr_D2KSPiPiPi_LL? ('1' for yes, others for no)"<<"\n";
    cin >> i3;
    cout << "Fit subtr_D2KSKpPiPi_LL? ('1' for yes, others for no)"<<"\n";
    cin >> i4;
    */
    if (i1==1) {
        all_D2KSKpPiPi_DD();
    }
    
    if (i2==1) {
        all_D2KSKpPiPi_LL();
    }
	/*
	if (i3==1) {
        fit_M_D2KSPiPiPi_LL();
    }
    
    if (i4==1) {
        fit_M_D2KSKpPiPi_LL();
    }
	*/
    return 0;
}
