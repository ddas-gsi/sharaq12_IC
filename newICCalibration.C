#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cmath>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"

// ==================================================================================

double ICSegment[30] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                        20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

double ElossLISE_78MeV[30] = {4.429, 4.433, 4.437, 4.441, 4.445, 4.449, 4.453, 4.457, 4.461, 4.465,
                              4.469, 4.473, 4.477, 4.481, 4.485, 4.489, 4.493, 4.493, 4.48, 4.485,
                              4.49, 4.495, 4.499, 4.504, 4.509, 4.514, 4.519, 4.523, 4.528, 4.533}; // for 2024 FastBeam run 0089

double scalingFactor2024[30] = {0.328100808214003, 0.354331023347641, 0.353376871615164, 0.341221667306953,
                                0.299871820818998, 0.2978110984671, 0.356931018451723, 0.329294421869228,
                                0.355599840573934, 0.351940599678406, 0.391313865417451, 0.408803019640458,
                                0.359500216808261, 0.409018301309844, 0.352190096273146, 0.357800432006759,
                                0.363726148939098, 0.352455737113362, 0.357382175563994, 0.354226230906535,
                                0.349125630797702, 0.345915579668321, 0.354422202790317, 0.35380711856152,
                                0.346675482839218, 0.352375450812634, 0.357493196633124, 0.353796092051126,
                                0.352195387547155, 0.353249222664682};

// From charge to voltage:

double IC_V_M[30] = {0.000168032, 0.000155687, 0.000148696, 0.000148154, 0.000149055, 0.000156954,
                     0.000162803, 0.000174383, 0.000156753, 0.000152314, 0.000163803, 0.000156695,
                     0.000171786, 0.000165882, 0.000163502, 0.000157136, 0.000159309, 0.000151592,
                     0.000173288, 0.000160501, 0.00015524, 0.000164984, 0.000158545, 0.000144887,
                     0.000162171, 0.00016246, 0.000177322, 0.00015614, 0.000162964, 0.000161043};

double IC_V_N[30] = {0.040754, 0.0463789, 0.0584141, 0.0502596, 0.0565706, 0.054101, 0.0711771,
                     0.054378, 0.0681855, 0.049819, 0.109847, 0.00333231, 0.116113, 0.00149862,
                     0.0468207, 0.0291755, 0.0266396, 0.0293459, 0.028094, 0.0401911, 0.0266147,
                     0.0313477, 0.0269507, 0.0282567, 0.029073, 0.0322126, 0.0310689, 0.0282302,
                     0.0214893, 0.0223735};

// From voltage to energy:

double IC_E_M[30] = {185.76, 211.80, 204.22, 208.80, 209.99, 203.18, 193.55, 193.99, 193.04,
                     207.48, 163.56, 195.58, 158.65, 177.20, 189.61, 205.84, 214.86, 213.22,
                     196.01, 198.81, 216.87, 199.13, 213.05, 223.90, 213.76, 204.70, 188.67,
                     212.10, 211.80, 188.89};

double IC_E_N[30] = {-7.57, -9.82, -11.93, -10.49, -11.88, -10.99, -13.78, -10.55, -13.16,
                     -10.34, -17.97, -0.65, -18.42, -0.27, -8.88, -6.01, -5.72, -6.26, -5.51,
                     -7.99, -5.77, -6.24, -5.74, -6.33, -6.21, -6.59, -5.86, -5.99, -4.55, -4.23};

// ==================================================================================

void newICCalibration()
{
    // TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    // TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/opticsOutput/FastBeam_0089_2024.root");
    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_1053new.root");
    TTree *tree = (TTree *)file->Get("tree_new");

    TCanvas *c1 = new TCanvas("c1", "IC Canvas", 1500, 1200);
    c1->Divide(2, 2);

    TList *hlist_E_Cal = new TList();

    TH2F *h_IC_C_UnCal = new TH2F("IC_C_UnCal", "IC UnCalibrated Charge Data; Segment; IC Uncalibrated Charge", 31, -0.5, 30.5, 1000, 0, 500); // IC Charge Data
    hlist_E_Cal->Add(h_IC_C_UnCal);

    TH2F *h_IC_E_UnCal = new TH2F("IC_E_UnCal", "IC UnCalibrated Energy Data; Segment; IC Uncalibrated Energy", 31, -0.5, 30.5, 1000, 0, 100); // IC Energy Data
    hlist_E_Cal->Add(h_IC_E_UnCal);

    TH2F *h_IC_E_Cal = new TH2F("IC_E_Cal", "IC Calibrated/Scaled Energy Data; Segment; IC Calibrated Energy", 31, -0.5, 30.5, 1000, 0, 100); // IC Energy Data
    hlist_E_Cal->Add(h_IC_E_Cal);

    Double_t IC_C[30];
    Double_t IC_E[30];
    Double_t IC_E_Cal[30];
    Double_t IC_E_Cal_Carlos[30];
    // vector<vector<double>> IC_E_Cal_vec;

    for (int k = 0; k < 30; k++)
    {
        string IC_C_SetBranch = "IC_C_" + to_string(k);
        tree->SetBranchAddress(IC_C_SetBranch.c_str(), &IC_C[k]); // .c_str() converts any string to CONSTANT

        string IC_E_SetBranch = "IC_E_" + to_string(k);
        tree->SetBranchAddress(IC_E_SetBranch.c_str(), &IC_E[k]); // .c_str() converts any string to CONSTANT

        string IC_E_Cal_SetBranch = "IC_E_Cal_" + to_string(k);
        tree->SetBranchAddress(IC_E_Cal_SetBranch.c_str(), &IC_E_Cal[k]); // .c_str() converts any string to CONSTANT
    }

    int entries = tree->GetEntries();
    cout << "Number of Entries: " << entries << endl;

    int count = 0;

    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);

        for (int ix = 0; ix < 30; ix++)
        {
            h_IC_C_UnCal->Fill(ix, IC_C[ix]);

            h_IC_E_UnCal->Fill(ix, IC_E[ix]);

            // h_IC_E_Cal->Fill(ix, IC_E[ix] * scalingFactor2024[ix]);
            h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);

            // // Using Carlos Params
            // IC_E_Cal_Carlos[ix] = (IC_C[ix] * IC_V_M[ix] + IC_V_N[ix]) * IC_E_M[ix] + IC_E_N[ix];
            // h_IC_E_Cal_Carlos->Fill(ix, IC_E_Cal_Carlos[ix]);
        }
    }

    c1->cd(1);
    h_IC_C_UnCal->Draw("colz");

    c1->cd(2);
    h_IC_E_UnCal->Draw("colz");

    c1->cd(3);
    h_IC_E_Cal->Draw("colz");

    // Overlay TGraph for ElossLISE_78MeV
    c1->cd(3); // Draw on top of h_IC_E_Cal
    TGraph *g_ElossLISE = new TGraph(30, ICSegment, ElossLISE_78MeV);
    g_ElossLISE->SetMarkerStyle(20);   // Solid circle marker
    g_ElossLISE->SetMarkerSize(1.0);   // Adjust size
    g_ElossLISE->SetMarkerColor(kRed); // Set color to red
    g_ElossLISE->Draw("P SAME");

    // c1->cd(4);
    // h_IC_E_Cal->Draw("colz");

    c1->Update();
}