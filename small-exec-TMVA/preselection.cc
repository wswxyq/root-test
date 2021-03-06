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

int preselection()
{

    //define TChain, TFile
    const char* TChain_name="D2KSPiPiPi_LL/DecayTree";
    const char* TChain_file="/eos/lhcb/user/m/mamartin/D2Kshhh/Test/MD/*.root";

    TChain *chainx= new TChain(TChain_name);
    chainx->Add(TChain_file);
    TFile *newfile = new TFile("../saved/temp.root","recreate");


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
    newfile->cd();
    newtree->Write("ReducedTree");
    return 0;
}


