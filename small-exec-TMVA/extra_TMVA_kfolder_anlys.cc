#include <iostream>
#include "TMVA/mvas.h"
void extra_TMVA_kfolder_anlys()
{
    TString fileName = "../plotdir/D2KSPiPiPi_DD/TMVA/TMVA_refine_1_1.root";
    TString dataset = "dataset";

    TMVA::mvas(dataset,fileName,TMVA::kMVAType);
    TMVA::mvas(dataset,fileName,TMVA::kProbaType);
    TMVA::mvas(dataset,fileName,TMVA::kRarityType);
    TMVA::mvas(dataset,fileName,TMVA::kCompareType);
}