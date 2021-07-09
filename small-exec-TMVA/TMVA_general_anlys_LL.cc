//analyse result from TMVA_general*.cc on overtraining etc.
#include <iostream>
#include "TMVA/mvas.h"
void TMVA_general_anlys_LL()
{
    TString fileName = "../plotdir/D2KSPiPiPi_LL/TMVA/TMVA_refine_LL_1.root";
    TString dataset = "D2KSPiPiPi_LL";

    TMVA::mvas(dataset,fileName,TMVA::kMVAType);
    TMVA::mvas(dataset,fileName,TMVA::kProbaType);
    TMVA::mvas(dataset,fileName,TMVA::kRarityType);
    TMVA::mvas(dataset,fileName,TMVA::kCompareType);
}