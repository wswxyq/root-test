//train model with s-weight using NDT method
//result stored in file TMVA_refine_*.root and folder named after tree name
#include <iostream> 
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMVA/Event.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"
void TMVA_general_DD()
{
    //Declare Factory
    TMVA::Tools::Instance();
    auto outputFile = TFile::Open("../plotdir/D2KSPiPiPi_DD/TMVA/TMVA_refine_DD_1.root", "RECREATE");
    TMVA::Factory factory("TMVAClassification", outputFile,
"!V:ROC:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" );

    //Define input variables
    TMVA::DataLoader * loader = new TMVA::DataLoader("D2KSPiPiPi_DD");
    loader->AddVariable( "D_PT", 'F' );
    //loader->AddVariable( "P0_ProbNNpi", 'F' );
    loader->AddVariable( "P1_ProbNNpi", 'F' );
    loader->AddVariable( "P2_ProbNNpi", 'F' );
    loader->AddVariable( "D_ReFit_chi2", 'F' );
    loader->AddVariable( "KS_PT", 'F' );
    loader->AddVariable( "var1:=D_ReFit_decayLength/D_ReFit_decayLengthErr", 'F' );
	loader->AddVariable( "var2:=log(D_IPCHI2_OWNPV)",'F');	


    //Setup input Dataset
    TString inputFileName = "../saved/TMVA_root_files/reduced_DD_1.root";
    auto inputFile = TFile::Open( inputFileName );

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

    /*
    // Use a kernel density estimator to approximate the PDFs
    factory.BookMethod(loader, TMVA::Types::kLikelihood, "LikelihoodKDE",
"!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:\
KDEborder=None:NAvEvtPerBin=50" ); 

    // Fisher discriminant (same as LD)
    factory.BookMethod(loader, TMVA::Types::kFisher, "Fisher", "!H:!V:Fisher:VarTransform=None:\
CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

    */
    //Multi-Layer Perceptron (Neural Network)
    /*factory.BookMethod(loader, TMVA::Types::kMLP, "MLP",
"!H:!V:NeuronType=tanh:VarTransform=N:NCycles=100:HiddenLayers=N+5:TestRate=5:!UseRegulator" ); */

    
    factory.BookMethod(loader,TMVA::Types::kBDT, "BDT","!V:NTrees=150:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:\
AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    auto c1 = factory.GetROCCurve(loader);
    c1->Draw();
	c1->Print("../plotdir/D2KSPiPiPi_DD/TMVA_out_DD.pdf");

    outputFile->Close();
}