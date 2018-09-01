//analyse result from TMVA_general*.cc on overtraining etc.
#include <iostream> 
#include "TMVA/mvas.h"
void TMVA_general_anlys_DD()
{
    TString fileName = "../plotdir/D2KSPiPiPi_DD/TMVA/TMVA_refine_DD_1.root";
    TString dataset = "D2KSPiPiPi_DD";

    TMVA::mvas(dataset,fileName,TMVA::kMVAType);
    TMVA::mvas(dataset,fileName,TMVA::kProbaType);
    TMVA::mvas(dataset,fileName,TMVA::kRarityType);
    TMVA::mvas(dataset,fileName,TMVA::kCompareType);
}