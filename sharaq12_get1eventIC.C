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

// Energy Loss for BeamEnergy 15.5 MeV
vector<double> EnergyLoss15_5MeV = {19.307, 19.777, 20.276, 20.82, 21.413, 22.06, 22.77,
                                    23.554, 24.425, 25.4, 26.502, 27.759, 29.211, 30.905,
                                    32.892, 35.247, 37.76, 35.479, 12.255, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0};

vector<double> Segment = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                          17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

const int n = 30;

void sharaq12_get1eventIC()
{
    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    TTree *tree = (TTree *)file->Get("tree_new");

    // TFile *graphFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/graphCuts/FE9_X_PID_cut.root");
    // TCutG *cut = (TCutG *)graphFile->Get("myCut");

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

    int count = 0;

    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);

        if (IC_E_Cal[0] > 19.30 && IC_E_Cal[0] < 19.31)
            // if (i == 419661)  // for the case if (IC_E_Cal[0] > 15.0 && IC_E_Cal[0] < 16.0)
            // if (i == 419762) // for the 15.5 MeV/u beam energy seems reasonable
            if (i == 363409)
            {
                cout << "Entry No.: " << i << endl;

                for (int ix = 0; ix < 30; ix++)
                {
                    cout << "IC Segment: " << ix << " >> EnergyLoss: " << IC_E_Cal[ix] << endl;

                    h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);
                }

                count = count + 1;
            }
    }
    cout << count << endl;

    // gStyle->SetPalette(kRed);
    // h_IC_E_Cal->SetFillColor(kRed);
    // h_IC_E_Cal->SetLineColor(kBlue);
    // h_IC_E_Cal->Draw("colz");
    h_IC_E_Cal->Draw("BOX");

    // Use .data() to pass raw pointers (like const Double_t*) to TGraph
    // TGraph *graph = new TGraph(n, Segment.data(), EnergyLoss15MeV.data());
    TGraph *graph = new TGraph(n, Segment.data(), EnergyLoss15_5MeV.data());
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kRed);
    graph->Draw("P SAME");

    c1->Update();
}