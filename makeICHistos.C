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

vector<double> scalingFactor = {0.577467478021203, 0.581816361603791, 0.564955317464165,
                                0.562139191789609, 0.561005004269991, 0.561218801267859, 0.567053744111581,
                                0.569979158925311, 0.560227642018658, 0.55741623032533, 0.558140994953815,
                                0.557823129251701, 0.555961278493404, 0.556697353011889, 0.555589846462464,
                                0.554531914566923, 0.556387512959337, 0.551068847585898, 0.556926626884374,
                                0.555737144160237, 0.546561722045824, 0.545918783911169, 0.550806012050982,
                                0.549184647472296, 0.540979901878595, 0.546395536631919, 0.551659412203756,
                                0.542413675610234, 0.54024933472174, 0.530419340460305};

vector<double> fastBeamCalibGain = {1.602698506, 1.389071469, 1.377926764, 1.360651156, 1.331318079, 1.321526422, 1.330078099,
                                    1.407405634, 1.315865412, 1.243962221, 1.454810263, 1.205535935, 1.421234305, 1.176428169,
                                    1.173315345, 0.9011144476, 0.8718705323, 0.9236897765, 0.8994656762, 0.9463624396,
                                    0.8797988209, 0.9286721693, 0.8944084016, 0.8927418722, 0.8460740386, 0.8595045522,
                                    0.8507238845, 0.8110059205, 0.7798595606, 0.7066028773};

vector<double> fastBeamCalibOffset = {-12.34962539, -9.676244073, -10.06149722, -9.959041221, -9.651483604, -9.546878672,
                                      -9.507970187, -10.36683731, -9.533352802, -8.732453079, -11.42562881, -8.283599067,
                                      -11.13727523, -7.990747165, -8.006340186, -4.514376126, -4.108062137, -4.913789452,
                                      -4.483116823, -5.138871074, -4.47090886, -5.156721535, -4.601901164, -4.628643422,
                                      -4.185189839, -4.265171009, -4.046910533, -3.704934567, -3.313521891, -2.490195728};

void makeICHistos()
{
    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    // TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All_newCalib.root");
    TTree *tree = (TTree *)file->Get("tree_new");

    // for uncalibrated data
    vector<vector<double>> IC_E_vec;
    Double_t IC_E[30];
    // for calibrated data
    vector<vector<double>> IC_E_Cal_vec;
    Double_t IC_E_Cal[30];

    for (int k = 0; k < 30; k++)
    {
        // Uncalibrated Energy branch
        string IC_E_SetBranch = "IC_E_" + to_string(k);
        tree->SetBranchAddress(IC_E_SetBranch.c_str(), &IC_E[k]); // .c_str() converts any string to CONSTANT

        // calibrated Energy branch
        string IC_E_Cal_SetBranch = "IC_E_Cal_" + to_string(k);
        tree->SetBranchAddress(IC_E_Cal_SetBranch.c_str(), &IC_E_Cal[k]); // .c_str() converts any string to CONSTANT
    }
    int entries = tree->GetEntries();
    cout << "Number of Entries: " << entries << endl;

    TCanvas *c1 = new TCanvas("c1", "IC Canvas", 1500, 600);
    c1->Divide(2, 1);

    TList *hlist_E = new TList();
    // for Uncalibrated Energy
    TH2F *h_IC_E = new TH2F("IC_E", "IC Uncalibrated Energy Data; Segment; IC UnCal Energy", 31, -0.5, 30.5, 1000, 0, 100); // IC UnCal Energy Data
    hlist_E->Add(h_IC_E);
    // for calibrated Energy
    TH2F *h_IC_E_Cal = new TH2F("IC_E_Cal", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 100); // IC Cal Energy Data
    hlist_E->Add(h_IC_E_Cal);

    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        for (int ix = 0; ix < 30; ix++)
        {
            // uncalibrated energy
            h_IC_E->Fill(ix, IC_E[ix]);
            // calibrated energy
            // h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);
            // h_IC_E_Cal->Fill(ix, IC_E[ix] * fastBeamCalibGain[ix] + fastBeamCalibOffset[ix]);
            h_IC_E_Cal->Fill(ix, IC_E[ix] * scalingFactor[ix]);
        }
    }
    c1->cd(1);
    h_IC_E->Draw("colz");

    c1->cd(2);
    h_IC_E_Cal->Draw("colz");

    c1->Update();
}