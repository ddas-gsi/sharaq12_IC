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
// LISE++ EnergyLoss Data for 15.5 MeV/u 50Ca20

std::vector<double> ICSegment = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                                 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

std::vector<double> EnergyLoss15_5MeV = {19.307, 19.777, 20.276, 20.82, 21.413, 22.06, 22.77, 23.554,
                                         24.425, 25.4, 26.502, 27.759, 29.211, 30.905, 32.892, 35.247,
                                         37.76, 35.479, 12.255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// ==================================================================================

// Define the constants
double c = 300000000;   // velocity in m/s
double amu_E = 931.494; // in MeV/u

// Distance
double distance_F3_FE9 = 68.54307; // distance in meter
// double distance_FE9_FE12 = 14.87168; // distance in meter
// double distance_FE9_FE12 = 14.68793; // distance in meter   // Updated distance (2023.05.01)
// double distance_FE12_S1 = 9.163;     // distance in meter    // Updated distance (2023.05.01) >> FE12-S0 = 1.350m , S0-S1 = 7.813m
double distance_FE9_FE12 = 14.87168; // distance in meter   // Carlos used this distance
double distance_FE12_S1 = 9.4892;    // distance in meter   // Carlos used this distance

double calculateBeamEnergy(double TOF, double pathLength)
{
    /*
    Returns kinetic Energy in MeV/u for every event from the TOF information
    */
    double kineticE;
    if (TOF > 0)
    {
        double velocity = (pathLength / TOF) * pow(10, 9);
        double beta = velocity / c;
        double gamma = 1 / sqrt(1 - pow(beta, 2));
        double totalE = gamma * amu_E;
        kineticE = totalE - amu_E;
    }
    else if (TOF <= 0)
    {
        kineticE = 0;
    }
    return kineticE;
}

void sharaq12_ICAnalysis()
{

    // ==================================================================================
    // Input File

    // TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    // TTree *tree = (TTree *)file->Get("tree_new");

    // TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/sharaq12phys_1003new.root");
    // TTree *tree = (TTree *)file->Get("tree_new");

    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/sharaq12phys_1005new.root");
    TTree *tree = (TTree *)file->Get("tree_new");

    // ==================================================================================

    // ==================================================================================
    // Graph Cut Files

    // TFile *graphFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/graphCuts/FE9_X_PID_cut.root", "READ");
    // TCutG *cut = (TCutG *)graphFile->Get("myCut");

    // TFile *graphFileFE9 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/graphCuts/FE9_X_PID_cut_50Ca_051124.root", "READ");
    // TCutG *FE9cut = (TCutG *)graphFileFE9->Get("newCut0511");

    // TFile *graphFileFE9_50Ca20 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/graphCuts/FE9_X_PID_cut_50Ca20_170125.root", "READ");
    // TCutG *cut50Ca20FE9 = (TCutG *)graphFileFE9_50Ca20->Get("cut50Ca20FE9");

    // TFile *graphFileS1 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/graphCuts/S1_PID_cut_50Ca_150125.root", "READ");
    // TCutG *S1cut = (TCutG *)graphFileS1->Get("s1pidCut150125");

    TFile *graphFileFE9_50Ca20 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/1005_FE9_PID_50Ca20_2401.root", "READ");
    TCutG *cut50Ca20FE9 = (TCutG *)graphFileFE9_50Ca20->Get("FE9pidCut2401_1005");

    TFile *graphFileS1 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/1005_S1_PID_50Ca20_2401.root", "READ");
    TCutG *S1cut = (TCutG *)graphFileS1->Get("s1pidCut2401_1005");

    // ==================================================================================

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

    Double_t FE9_E_0;  // Energy in MeV/u at FE9. This comes from PID_T_0 TOF.
    Double_t FE12_E_0; // Energy in MeV/u at FE12. This comes from Beam_T_0 TOF.
    Double_t S1_E_0;   // Energy in MeV/u at S1. This comes from S1PID_T_0 TOF.

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

    tree->SetBranchAddress("FE9_E_0", &FE9_E_0);
    tree->SetBranchAddress("FE12_E_0", &FE12_E_0);
    tree->SetBranchAddress("S1_E_0", &S1_E_0);

    TCanvas *c1 = new TCanvas("c1", "IC Canvas", 1800, 1350);
    c1->Divide(3, 3);

    TList *hlist_E_Cal = new TList();

    // TH2F *h_IC_E_Cal = new TH2F("IC_E_Cal", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60); // IC Energy Data
    // hlist_E_Cal->Add(h_IC_E_Cal);
    TH2F *h_IC_E_2Dcut = new TH2F("IC_E_2Dcut", "IC Uncalibrated Energy Data; Segment; IC UnCal Energy", 31, -0.5, 30.5, 1000, 0, 100); // IC Uncalibrated Energy Data
    hlist_E_Cal->Add(h_IC_E_2Dcut);
    TH2F *h_IC_E_Cal_2Dcut = new TH2F("IC_E_Cal_2Dcut", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60); // IC Energy Data
    hlist_E_Cal->Add(h_IC_E_Cal_2Dcut);

    TH2F *h_IC_E_Cal_UpdatedC_2Dcut = new TH2F("IC_E_Cal_UpdatedC_2Dcut", "IC Calibrated Energy Data Updated Params; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60); // IC Energy Data
    hlist_E_Cal->Add(h_IC_E_Cal_UpdatedC_2Dcut);

    Double_t IC_E[30];
    Double_t IC_E_Cal[30];
    Double_t IC_C[30];
    Double_t IC_E_Cal_UpdatedC[30];
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

    TH2F *h_FE9_PID = new TH2F("FE9_PID", "FE9_X_0:PID_T_0; TOF (ns); X Position", 1000, 1050, 1150, 1000, -50, 50);
    TH2F *h_S1_PID = new TH2F("S1_PID", "S1_X_0:S1PID_T_0; TOF (ns); X Position", 1000, 600, 700, 1000, -200, 200);

    TH1F *h_FE9_E_0 = new TH1F("FE9_E_0", "Beam Energy at FE9; Energy (MeV/u); Counts", 1000, 15, 30);
    TH1F *h_FE12_E_0 = new TH1F("FE12_E_0", "Beam Energy at FE12; Energy (MeV/u); Counts", 1000, 1, 40);
    TH1F *h_S1_E_0 = new TH1F("S1_E_0", "Beam Energy at S1; Energy (MeV/u); Counts", 1000, 1, 30);
    TH1F *h_S1_A_0 = new TH1F("S1_A_0", "S1 horizontal Angle; Horizontal Angle; Counts", 1000, -0.5, 0.5);
    TH1F *h_S1_B_0 = new TH1F("S1_B_0", "S1 Vertical Angle; Vertical Angle; Counts", 1000, -0.5, 0.5);
    TH1F *h_S1_X_0 = new TH1F("S1_X_0", "S1 X Position; Horizontal Position; Counts", 1000, -200, 200);
    TH1F *h_S1_Y_0 = new TH1F("S1_Y_0", "S1 Y Position; Vertical Position; Counts", 1000, -200, 200);

    TH1F *h_S1_E_0_cut = new TH1F("S1_E_0_cut", "Beam Energy at S1 Cut; Energy (MeV/u); Counts", 1000, 1, 30);
    TH1F *h_S1_A_0_cut = new TH1F("S1_A_0_cut", "S1 horizontal Angle Cut; Horizontal Angle; Counts", 1000, -0.5, 0.5);
    TH1F *h_S1_B_0_cut = new TH1F("S1_B_0_cut", "S1 Vertical Angle Cut; Vertical Angle; Counts", 1000, -0.5, 0.5);
    TH1F *h_S1_X_0_cut = new TH1F("S1_X_0_cut", "S1 X Position Cut; Horizontal Position; Counts", 1000, -200, 200);
    TH1F *h_S1_Y_0_cut = new TH1F("S1_Y_0_cut", "S1 Y Position Cut; Vertical Position; Counts", 1000, -200, 200);

    vector<double> PID_T_0_vec;   // Time of Flight between F3 and FE9. Time in nanosecond.
    vector<double> Beam_T_0_vec;  // Time of Flight between FE9 and FE12. Time in nanosecond.
    vector<double> S1PID_T_0_vec; // Time at S1 (average between anodes). Time in nanosecond.

    vector<double> FE9_E_0_vec;  // Energy in MeV/u at FE9. This comes from PID_T_0 TOF.
    vector<double> FE12_E_0_vec; // Energy in MeV/u at FE12. This comes from Beam_T_0 TOF.
    vector<double> S1_E_0_vec;   // Energy in MeV/u at S1. This comes from S1PID_T_0 TOF.

    int entries = tree->GetEntries();
    cout << "Number of Entries: " << entries << endl;

    int count = 0;

    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        // count++;

        // double FE9_E_0 = calculateBeamEnergy(PID_T_0, distance_F3_FE9);
        // h_FE9_E_0->Fill(FE9_E_0);
        // FE9_E_0_vec.push_back(FE9_E_0);

        // double FE12_E_0 = calculateBeamEnergy(Beam_T_0 + 22.2, distance_FE9_FE12); // 22.2 ns is the TOF offset for Beam_T_0
        // h_FE12_E_0->Fill(FE12_E_0);
        // FE12_E_0_vec.push_back(FE12_E_0);

        // double S1_E_0 = calculateBeamEnergy(S1PID_T_0 - 467.679, distance_FE12_S1); // 467.679 ns is the TOF offset for S1PID_T_0
        // h_S1_E_0->Fill(S1_E_0);
        // S1_E_0_vec.push_back(S1_E_0);

        h_FE9_PID->Fill(PID_T_0, FE9_X_0);
        h_S1_PID->Fill(S1PID_T_0, S1_X_0);

        h_FE9_E_0->Fill(FE9_E_0);
        h_FE12_E_0->Fill(FE12_E_0);
        h_S1_E_0->Fill(S1_E_0);
        h_S1_X_0->Fill(S1_X_0);
        h_S1_Y_0->Fill(S1_Y_0);
        h_S1_A_0->Fill(S1_A_0);
        h_S1_B_0->Fill(S1_B_0);

        // if (fabs(S1_E_0 - 15.50) < 0.5)
        // if (FE9cut->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (S1_E_0 > 15.50 && S1_E_0 < 15.55) && (fabs(S1_Y_0 - 6.115) < 2))
        // if (FE9cut->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.1) && (fabs(S1_X_0 - 11.86) < 2))
        if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0))
        // if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && not S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.5) && (fabs(S1_Y_0 - 6.115) < 3) && (fabs(S1_X_0 - 11.86) < 15) && (fabs(S1_A_0 < -0.01747) < 0.04))
        // if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.1) && (fabs(S1_Y_0 - (50)) < 5) && (fabs(S1_X_0 - 20) < 15) && (fabs(S1_A_0 - (-0.01747)) < 0.04) && (fabs(S1_B_0 - 0.003003) < 0.012))
        // if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.1) && (fabs(S1_Y_0 - (6.792)) < 5) && (fabs(S1_X_0 - 2.798) < 15) && (fabs(S1_A_0 - (0.013)) < 0.04) && (fabs(S1_B_0 - 0.0033) < 0.012))
        // if (FE9cut->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.5) && (fabs(S1_Y_0 - 6.115) < 2))
        {

            h_S1_E_0_cut->Fill(S1_E_0);
            h_S1_X_0_cut->Fill(S1_X_0);
            h_S1_Y_0_cut->Fill(S1_Y_0);
            h_S1_A_0_cut->Fill(S1_A_0);
            h_S1_B_0_cut->Fill(S1_B_0);

            cout << "Event: " << i << "  " << S1_E_0 << " " << S1_X_0 << " " << S1_Y_0 << " " << S1_A_0 << " " << S1_B_0 << endl;
            // count++;

            for (int ix = 0; ix < 30; ix++)
            {
                // h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);
                // Gates are applied here
                // if (cut->IsInside(PID_T_0, FE9_X_0) && fabs(Beam_T_0 - 245.0) < 1)
                //     h_IC_E_Cal_2Dcut->Fill(ix, IC_E_Cal[ix]);

                // if ((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3]) / 4 > 10)
                if (((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3] + IC_E[4] + IC_E[5] + IC_E[6] + IC_E[7]) / 8.0) > 7)
                {
                    h_IC_E_2Dcut->Fill(ix, IC_E[ix]);
                    h_IC_E_Cal_2Dcut->Fill(ix, IC_E_Cal[ix]);
                    count++;

                    // Using Updated Params
                    IC_E_Cal_UpdatedC[ix] = (IC_C[ix] * IC_V_M[ix] + IC_V_N[ix]) * IC_E_M[ix] + IC_E_N[ix];
                    h_IC_E_Cal_UpdatedC_2Dcut->Fill(ix, IC_E_Cal_UpdatedC[ix]);
                }
            }
        }
    }

    // cout << count << endl;
    cout << count / 30 << endl;

    // h_IC_E_Cal->Draw("colz");
    // h_IC_E_Cal_2Dcut->Draw("colz");
    // hlist_E_Cal->Draw("colz");

    // TFile *icAnalysisFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/902_icHistos.root", "recreate");
    // hlist_E_Cal->Write();
    // c1->Write();
    // icAnalysisFile->Close();

    const int n = 30;
    // Use .data() to pass raw pointers (like const Double_t*) to TGraph
    TGraph *graph = new TGraph(n, ICSegment.data(), EnergyLoss15_5MeV.data());
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kRed);

    // c1->SetLogy(); // Set log scale for Y-axis

    c1->cd(1);
    h_FE9_PID->Draw("colz");
    cut50Ca20FE9->Draw("SAME");
    cut50Ca20FE9->SetLineColor(kRed); // Set line color to red
    cut50Ca20FE9->SetLineWidth(2);    // Set line width to 2

    c1->cd(2);
    h_S1_PID->Draw("colz");
    S1cut->Draw("SAME");
    S1cut->SetLineColor(kRed); // Set line color to red
    S1cut->SetLineWidth(2);    // Set line width to 2

    c1->cd(3);
    // h_FE9_E_0->Draw();
    // h_FE12_E_0->Draw();
    h_S1_E_0->Draw();
    h_S1_E_0_cut->SetLineColor(kRed); // Set line color to red
    h_S1_E_0_cut->Draw("SAME");

    c1->cd(4);
    // h_S1_E_0->Draw();
    // h_IC_E_2Dcut->Draw("colz");
    h_S1_X_0->Draw();
    h_S1_X_0_cut->SetLineColor(kRed); // Set line color to red
    h_S1_X_0_cut->Draw("SAME");

    c1->cd(5);
    // h_IC_E_2Dcut->Draw("colz");
    // h_IC_E_Cal_2Dcut->Draw("colz");
    // graph->Draw("P SAME");
    h_S1_Y_0->Draw();
    h_S1_Y_0_cut->SetLineColor(kRed); // Set line color to red
    h_S1_Y_0_cut->Draw("SAME");

    c1->cd(6);
    // h_IC_E_Cal_2Dcut->Draw("colz");
    // h_IC_E_Cal_UpdatedC_2Dcut->Draw("colz");
    // graph->Draw("P SAME");
    h_S1_A_0->Draw("colz");
    h_S1_A_0_cut->SetLineColor(kRed); // Set line color to red
    h_S1_A_0_cut->Draw("SAME");

    c1->cd(7);
    h_S1_B_0->Draw("colz");
    h_S1_B_0_cut->SetLineColor(kRed); // Set line color to red
    h_S1_B_0_cut->Draw("SAME");

    c1->cd(8);
    // h_IC_E_2Dcut->Draw("colz");
    h_IC_E_Cal_2Dcut->Draw("colz");
    graph->Draw("P SAME");

    c1->cd(9);
    // h_IC_E_Cal_2Dcut->Draw("colz");
    h_IC_E_Cal_UpdatedC_2Dcut->Draw("colz");
    graph->Draw("P SAME");

    c1->Update();
}