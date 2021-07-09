#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TVector3.h"
using namespace TMVA;
using namespace std;
void TMVA_GetResponse_DD_x()
{

    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    //set origin file used to train
    TString fname = "../saved/TMVA_root_files/reduced_DD_2.root";
    TFile * input = new TFile( fname );


    Float_t var[8];
    //Define input variables
    reader->AddVariable( "D_PT", &var[0] );
    reader->AddVariable( "P0_ProbNNpi", &var[1] );
    reader->AddVariable( "P1_ProbNNpi", &var[2] );
    reader->AddVariable( "P2_ProbNNpi", &var[3] );
    reader->AddVariable( "D_ReFit_chi2", &var[4] );
    reader->AddVariable( "KS_PT", &var[5] );
    reader->AddVariable( "var1:=D_ReFit_decayLength/D_ReFit_decayLengthErr", &var[6] );
    reader->AddVariable( "var2:=log(D_IPCHI2_OWNPV)", &var[7] );


    // Book the MVA methods
    TString dir    = "D2KSPiPiPi_DD/weights/";
    TString prefix = "TMVAClassification";
    TString methodName = TString("BDT") + TString(" method");
    TString weightfile = dir + prefix + TString("_") + TString("BDT") + TString(".weights.xml");
    reader->BookMVA( methodName, weightfile );


    // Book output histograms
    UInt_t nbin = 100;
    TH1F   *histBdt(0);
    histBdt     = new TH1F( "MVA_BDT",  "MVA_BDT",  nbin, -0.8, 0.8 );


    // Event loop
    TTree* theTree = (TTree*)input->Get("ReducedTree");
    // output files
    TFile *target  = new TFile( "../saved/TMVA_root_files/TMVApp_DDx.root","RECREATE" );
    TTree *tree = theTree->CopyTree("");

    Double_t BDT_response;
    TBranch *BDT_response_branch = tree->Branch("BDT_response", &BDT_response, "BDT_response/D");
    Double_t Ct;
    TBranch *Ct_branch = tree->Branch("Ct", &Ct, "Ct/D");

    Double_t userVard[9];
    Float_t userVarf[8];

    Double_t P0_array[3], P1_array[3], P2_array[3];

    Int_t id_temp;
    TVector3 v1, v2, v3;

    theTree->SetBranchAddress( "D_PT", &userVard[0] );
    theTree->SetBranchAddress( "P0_ProbNNpi", &userVard[1] );
    theTree->SetBranchAddress( "P1_ProbNNpi", &userVard[2] );
    theTree->SetBranchAddress( "P2_ProbNNpi", &userVard[3] );
    theTree->SetBranchAddress( "D_ReFit_chi2", &userVarf[0] );
    theTree->SetBranchAddress( "KS_PT", &userVarf[1] );
    theTree->SetBranchAddress( "D_ReFit_decayLength", &userVarf[2] );
    theTree->SetBranchAddress( "D_ReFit_decayLengthErr", &userVarf[3] );
    theTree->SetBranchAddress( "D_IPCHI2_OWNPV", &userVard[4] );

    theTree->SetBranchAddress( "P0_ID", &id_temp );
    theTree->SetBranchAddress( "P0_PX", &P0_array[0] );
    theTree->SetBranchAddress( "P0_PY", &P0_array[1] );
    theTree->SetBranchAddress( "P0_PZ", &P0_array[2] );
    theTree->SetBranchAddress( "P1_PX", &P1_array[0] );
    theTree->SetBranchAddress( "P1_PY", &P1_array[1] );
    theTree->SetBranchAddress( "P1_PZ", &P1_array[2] );
    theTree->SetBranchAddress( "P2_PX", &P2_array[0] );
    theTree->SetBranchAddress( "P2_PY", &P2_array[1] );
    theTree->SetBranchAddress( "P2_PZ", &P2_array[2] );


    TStopwatch sw;
    sw.Start();
    for (Long64_t ievt=0; ievt<theTree->GetEntries(); ievt++) {
        //if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;



        theTree->GetEntry(ievt);
        var[0]=userVard[0];
        var[1]=userVard[1];
        var[2]=userVard[2];
        var[3]=userVard[3];
        var[4]=userVarf[0];
        var[5]=userVarf[1];
        var[6]=userVarf[2]/userVarf[3];
        var[7]=TMath::Log(userVard[4]);
        // Return the MVA outputs and fill into histograms
        histBdt    ->Fill( reader->EvaluateMVA( "BDT method" ) );

        BDT_response=reader->EvaluateMVA( "BDT method" );
        BDT_response_branch->Fill();

        v1.SetXYZ(P0_array[0],P0_array[1],P0_array[2]);
        v2.SetXYZ(P1_array[0],P1_array[1],P1_array[2]);
        v3.SetXYZ(P2_array[0],P2_array[1],P2_array[2]);
        Ct=v1.Dot(v2.Cross(v3));
        Ct_branch->Fill();
        if (ievt<100) {
            //cout<<ievt<<"  "<<id_temp<<"  "<<P0_array[0]<<"  "<<P0_array[1]<<"  "<<P0_array[2]<<endl;
            //cout<<v1(0)<<"  "<<v1(1)<<"  "<<v1(2)<<endl;

        }


    }
    cout    << "Number of event after cut "<<  tree->GetEntries("BDT_response>0.2433")  << endl;

    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: ";
    sw.Print();


    target->cd();
    tree->Write("MVATree");
    histBdt ->Write();


    //calculate significance

    Double_t cut_up=tree->GetMaximum("BDT_response");
    Double_t cut_down=tree->GetMinimum("BDT_response");

    Long64_t num_cut=100;
    Double_t cut_value[ num_cut+1 ], significance_value[ num_cut+1 ];

    Double_t readsweight, readbweight;
    Double_t sum_sweight, sum_bweight;

    tree->SetBranchAddress("swei", &readsweight );
    tree->SetBranchAddress("bwei", &readbweight );


    for(Long64_t i = 0; i <= num_cut; i++)
    {
        cut_value[i] = cut_down + (cut_up-cut_down)*i/num_cut;

        sum_bweight=0.;
        sum_sweight=0.;

        for(Long64_t ii = 0; ii < tree->GetEntries(); ii++)
        {
            tree->GetEntry(ii);
            if(BDT_response>=cut_value[i]) {
                sum_sweight+=readsweight;
                sum_bweight+=readbweight;
            }
        }

        significance_value[i] = sum_sweight/sqrt(sum_sweight+sum_bweight);

        //std::cout<<cut_value[i]<<"   "<<sum_sweight<<"   "<<sum_bweight
        //<<" "<<significance_value[i]<<std::endl;

    }

    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    c1->SetFillColor(42);
    c1->SetGrid();

    TGraph *gr = new TGraph(num_cut+1, cut_value, significance_value);
    gr->SetLineColor(2);
    gr->SetLineWidth(4);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->SetTitle("S vs Response");
    gr->GetXaxis()->SetTitle("BDT cut");
    gr->GetYaxis()->SetTitle("Significance");
    gr->Draw("ACP");

    // TCanvas::Update() draws the frame, after which one can change it
    c1->Update();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    c1->Modified();
    c1->Print("graph_DDx.eps");
    c1->Print("graph_DDx.png");

    target  ->Close();
    delete reader;

    Double_t bestcut;
    bestcut=cut_value[distance(significance_value, max_element(significance_value, significance_value + num_cut))];
    cout    << "The best BDT cut value: "<<  bestcut  << endl;



}

