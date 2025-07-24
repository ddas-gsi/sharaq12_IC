#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <stdio.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"

void sharaq12_ICAnalysis_v2()
{
    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    TTree *tree = (TTree *)file->Get("tree_new");

    TFile *graphFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/graphCuts/FE9_X_PID_cut.root", "READ");
    TCutG *cut = (TCutG *)graphFile->Get("myCut");

    Double_t F3_T_0;  // Time at F3 Diamond (array).
    Double_t PID_T_0; // Time of Flight between F3 and FE9.
    Double_t FE9_X_0; // Horizontal position at FE9 focal plane.
    Double_t FE9_Y_0; // Vertical position at FE9 focal plane.
    Double_t FE9_A_0; // Horizontal angle at FE9 focal plane.
    Double_t FE9_B_0; // Vertical angle at FE9 focal plane.

    Double_t Beam_T_0; // Time of Flight between FE9 and FE12.
    Double_t S0_X_0;   // Horizontal position at FE12/S0 focal plane.
    Double_t S0_Y_0;   // Vertical position at FE12/S0 focal plane.
    Double_t S0_A_0;   // Horizontal angle at FE12/S0 focal plane.
    Double_t S0_B_0;   // Vertical angle at FE12/S0 focal plane.

    Double_t S1PID_T_0; // Time at S1 (average between anodes).
    Double_t S1_X_0;    // Horizontal position at S1 focal plane.
    Double_t S1_Y_0;    // Vertical position at S1 focal plane.
    Double_t S1_A_0;    // Horizontal angle at S1 focal plane.
    Double_t S1_B_0;    // Vertical angle at S1 focal plane.

    tree->SetBranchAddress("F3_T_0", &F3_T_0);
    tree->SetBranchAddress("PID_T_0", &PID_T_0);
    tree->SetBranchAddress("FE9_X_0", &FE9_X_0);
    tree->SetBranchAddress("FE9_Y_0", &FE9_Y_0);
    tree->SetBranchAddress("FE9_A_0", &FE9_A_0);
    tree->SetBranchAddress("FE9_B_0", &FE9_B_0);

    tree->SetBranchAddress("Beam_T_0", &Beam_T_0);
    tree->SetBranchAddress("S0_X_0", &S0_X_0);
    tree->SetBranchAddress("S0_Y_0", &S0_Y_0);
    tree->SetBranchAddress("S0_A_0", &S0_A_0);
    tree->SetBranchAddress("S0_B_0", &S0_B_0);

    tree->SetBranchAddress("S1PID_T_0", &S1PID_T_0);
    tree->SetBranchAddress("S1_X_0", &S1_X_0);
    tree->SetBranchAddress("S1_Y_0", &S1_Y_0);
    tree->SetBranchAddress("S1_A_0", &S1_A_0);
    tree->SetBranchAddress("S1_B_0", &S1_B_0);

    TCanvas *c1 = new TCanvas("c1", "IC Canvas", 800, 600);

    TList *hlist_E_Cal = new TList();

    TH2F *h_IC_E_Cal = new TH2F("IC_E_Cal", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 50); // IC Energy Data
    hlist_E_Cal->Add(h_IC_E_Cal);

    vector<vector<double>> IC_E_Cal_vec;
    Double_t IC_E_Cal[30];

    for (int k = 0; k < 30; k++)
    {
        string IC_E_SetBranch = "IC_E_Cal_" + to_string(k);
        tree->SetBranchAddress(IC_E_SetBranch.c_str(), &IC_E_Cal[k]); // .c_str() converts any string to CONSTANT
    }

    int entries = tree->GetEntries();
    cout << "Number of Entries: " << entries << endl;

    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        for (int ix = 0; ix < 30; ix++)
        {
            h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);
        }
    }
    h_IC_E_Cal->Draw("colz");
    // h_IC_E_Cal_2Dcut->Draw("colz");
    // hlist_E_Cal->Draw("colz");

    TFile *icAnalysisFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/902_icHistos.root", "recreate");

    c1->Update();

    cout << "Do you want to continue more Analysis (Y/N)?: " << endl;
    std::string moreAnalysis;
    std::getline(std::cin, moreAnalysis);

    int k = 1;

    while (moreAnalysis == "Y" || moreAnalysis == "y")
    {
        string hist2Dcut_iName = "IC_E_Cal_2Dcut_" + to_string(k);

        TH2F *h_IC_E_Cal_2Dcut = new TH2F(hist2Dcut_iName.c_str(), "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 50); // IC Energy Data
        hlist_E_Cal->Add(h_IC_E_Cal_2Dcut);                                                                                                           // maybe need to make h_IC_E_Cal_2Dcut with dynamic naming scheme

        k = k + 1;
        // cout << k <<endl;

        cout << "***** Give the Gate Condition *****" << endl;

        cout << "Is it a Graph Cut (use flag 'graph' or 'g') or 1D Cut (use flag '1d')?" << endl;
        std::string cutType;
        std::getline(std::cin, cutType);

        if (cutType == "graph" || cutType == "g")
        {
            cout << "Give the graphCut ROOT file with Path:" << endl;
            std::string graphCutFileName;
            std::getline(std::cin, graphCutFileName); // Get file name from terminal input
            std::string graphCutFileWithPath = "/u/ddas/software/work/artemis-oedo/output/Analysis/graphCuts/" + graphCutFileName;

            TFile *graphCutFile = new TFile(graphCutFileWithPath.c_str(), "READ"); // Open the ROOT file using the provided file name
            if (!graphCutFile || graphCutFile->IsZombie())
            {
                std::cerr << "Error: Could not open file " << graphCutFile << std::endl;
                // return 1;
            }
            TCutG *cut = (TCutG *)graphCutFile->Get("myCut");

            cout << "Give the X-axis Branch Name:" << endl;
            std::string branchNameX;
            std::getline(std::cin, branchNameX);

            cout << "Give the Y-axis Branch Name:" << endl;
            std::string branchNameY;
            std::getline(std::cin, branchNameY);
        }

        else if (cutType == "1d")
        {
            cout << "Give the Branch Name:" << endl;
            std::string branchName;
            std::getline(std::cin, branchName);
            cout << "Give the Branch Value:" << endl;
            std::string branchValue;
            std::getline(std::cin, branchValue);
        }

        cout << "Gate Condition loaded" << endl;

        for (int i = 0; i < entries; i++)
        {
            tree->GetEntry(i);
            for (int ix = 0; ix < 30; ix++)
            {
                if (cutType == "graph" || cutType == "g")
                {
                    // if (cut->IsInside(PID_T_0, FE9_X_0) && fabs(Beam_T_0 - 245.0) < 1)
                    if (cut->IsInside(PID_T_0, FE9_X_0))
                        h_IC_E_Cal_2Dcut->Fill(ix, IC_E_Cal[ix]);
                }
                // else if (cutType == "1d")
                // {
                //     if (branchName == "Beam_T_0")
                //     {
                //         if (fabs(Beam_T_0 - stof(branchValue)) < 5)
                //         {
                //             h_IC_E_Cal_2Dcut->Fill(ix, IC_E_Cal[ix]);
                //         }
                //     }
                // }
            }
        }
        h_IC_E_Cal_2Dcut->Draw("colz");

        c1->Update();

        cout << "sdfgsd" << endl;
        cout << "Do you want to continue more Analysis (Y/N)?: " << endl;
        std::getline(std::cin, moreAnalysis);

        // graphCutFile->Close();
    }

    icAnalysisFile->cd();
    hlist_E_Cal->Write();
    c1->Write();
    icAnalysisFile->Close();
}