#include <iostream> 
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMVA/Event.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"
#include "TMVA/CrossValidation.h"
#include "TMVA/Results.h"
void extra_TMVA_kfolder()
{
    //Setup input Dataset
    TString inputFileName = "../saved/TMVA_root_files/reduced_LL_1.root";
    auto inputFile = TFile::Open( inputFileName );
    //Declare Factory
    TMVA::Tools::Instance();
    auto outputFile = TFile::Open("../plotdir/D2KSPiPiPi_DD/TMVA/TMVA_refine_1_1.root", "RECREATE");

    //Define input variables
    TMVA::DataLoader * loader = new TMVA::DataLoader("dataset");
    loader->AddVariable( "D_PT", 'F' );
    loader->AddVariable( "P0_ProbNNpi", 'F' );
    loader->AddVariable( "P1_ProbNNpi", 'F' );
    //loader->AddVariable( "P2_ProbNNpi", 'F' );
    loader->AddVariable( "D_ReFit_chi2", 'F' );
    //loader->AddVariable( "KS_PT", 'F' );
    loader->AddVariable( "var1:=D_ReFit_decayLength/D_ReFit_decayLengthErr", 'F' );
	loader->AddVariable( "var2:=log(D_IPCHI2_OWNPV)",'F');	

    TTree *signalTree     = (TTree*)inputFile->Get("ReducedTree");
    TTree *backgroundTree = (TTree*)inputFile->Get("ReducedTree");

    Double_t signalWeight     = 1.;
    Double_t backgroundWeight = 1.;

    loader->AddSignalTree    ( signalTree,     signalWeight     );
    loader->AddBackgroundTree( backgroundTree, backgroundWeight );

    loader->SetSignalWeightExpression("swei");
    loader->SetBackgroundWeightExpression("bwei");

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_signal = "";
    TCut mycut_background = "";

    loader->PrepareTrainingAndTestTree( mycut_signal, mycut_background,
        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

    TString cvOptions = "!V:!Silent:ModelPersistence:AnalysisType=Classification:NumFolds=3";
        ":SplitExpr=""";
    auto cv = new TMVA::CrossValidation("TMVACrossValidation", loader, outputFile, cvOptions);

    // Use a kernel density estimator to approximate the PDFs
    //TMVA::CrossValidation CV(jobname, loader, options);

    cv->BookMethod(TMVA::Types::kBDT, "BDT","!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:\
        AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
    cv->Evaluate();
    TMVA::CrossValidationResult & result = (TMVA::CrossValidationResult &) cv->GetResults()[0];
    result.Print();
    auto c = new TCanvas();
    result.GetROCCurves()->Draw("AL");
    c->BuildLegend();
    c->Draw();
	c->Print("../plotdir/D2KSPiPiPi_DD/TMVA_out_DD_kfold.pdf");
    std::cout << "Average ROC Integral = " << result.GetROCAverage() 
          << " +/- " << result.GetROCStandardDeviation()/sqrt(cv->GetNumFolds()) << std::endl;
    outputFile->Close();



}





