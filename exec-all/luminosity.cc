#include <iostream>
#include <fstream>
#include "TChain.h"
//#include "ofstream.h"


void luminosity() {
    const char* TChain_name="GetIntegratedLuminosity/LumiTuple";
    const char* TChain_file1="/eos/lhcb/user/m/mamartin/D2KsKPiPi/2016/MU/*.root";
    const char* TChain_file2="/eos/lhcb/user/m/mamartin/D2KsKPiPi/2016/MU/*.root";
    //const char* TChain_file="../../root_data/DV*.root";
    TChain *chain= new TChain(TChain_name);
    chain->Add(TChain_file1);
    chain->Add(TChain_file2);
    Long64_t numentries= chain->GetEntries();
    Double_t lmn, sum_lmn=0;
    chain->SetBranchAddress("IntegratedLuminosity", &lmn);
    for(Long64_t i = 0; i < numentries; i++)
    {
        chain->GetEntry(i);
        sum_lmn+=lmn;

    }
    std::cout<<sum_lmn<<std::endl;

    std::ofstream myfile;
    myfile.open ("luminosity.txt");
    myfile << "luminosity:"<<sum_lmn<<"\n";
    myfile.close();
}