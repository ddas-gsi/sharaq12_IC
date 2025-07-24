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
#include "TCutG.h"

// ==================================================================================
double scalingFactor2024[30] = {0.328100808214003, 0.354331023347641, 0.353376871615164, 0.341221667306953,
                                0.299871820818998, 0.2978110984671, 0.356931018451723, 0.329294421869228,
                                0.355599840573934, 0.351940599678406, 0.391313865417451, 0.408803019640458,
                                0.359500216808261, 0.409018301309844, 0.352190096273146, 0.357800432006759,
                                0.363726148939098, 0.352455737113362, 0.357382175563994, 0.354226230906535,
                                0.349125630797702, 0.345915579668321, 0.354422202790317, 0.35380711856152,
                                0.346675482839218, 0.352375450812634, 0.357493196633124, 0.353796092051126,
                                0.352195387547155, 0.353249222664682};


double slowBeamScalingFactor2024[30] = {1.26338526912181, 1.193116359447, 1.1918085943246, 1.17401841379863, 
                                        1.06834068916363, 1.05059082252175, 1.11467740180611, 1.18675880416648,
                                        1.09826208630759, 1.05018386329223, 1.12629919107201, 0.966855412860375,
                                        1.1792201616738, 0.981287811996473, 0.988117252869986, 0.94271901552886,
                                        0.884537994920774, 0.920783898428363, 0.875873047148435, 0.954646910516312,
                                        0.969882242924157, 1, 1, 1, 1, 1, 1, 1, 1, 1};


double slowBeamScalingFactor2024v2[30] = {1.24356918787034, 1.19642915750435, 1.18418962285133, 1.17846597626726, 
                                          1.06987163297353, 1.04709269095848, 1.11999037978005, 1.18380290005004, 
                                          1.09117010023674, 1.04917247243532, 1.12979603829485, 0.96579216547196, 
                                          1.17231337956594, 0.980110832763987, 0.99499625733323, 0.920454212978164, 
                                          0.876854381481155, 0.911032622619079, 0.877540683403249, 0.947080149187554, 
                                          0.989079619782111, 1.08265621886245, 1, 1, 1, 1, 1, 1, 1, 1 };

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

std::vector<double> EnergyLoss16_5MeV = {17.839, 18.227, 18.635, 19.038, 19.486, 19.969, 20.49, 21.055, 
                                         21.663, 22.333, 23.071, 23.888, 24.799, 25.822, 26.984, 28.314, 
                                         29.849, 31.649, 33.784, 36.312, 38.167, 28.794, 3.757, 0, 0, 0, 0, 0, 0, 0};
    

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

void reverse_sharaq12_ICAnalysis()
{

    // ==================================================================================
    // Input File

    // TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    // TTree *tree = (TTree *)file->Get("tree_new");

    // TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/sharaq12phys_1003new.root");
    // TTree *tree = (TTree *)file->Get("tree_new");

    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_1053new.root");
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

    // 50Ca20 Cut Files ---------------------------------------

    // TFile *graphFileFE9_50Ca20 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/1005_FE9_PID_50Ca20_2401.root", "READ");
    // TCutG *FE9cut50Ca20 = (TCutG *)graphFileFE9_50Ca20->Get("FE9pidCut2401_1005");
    // TFile *graphFileS1_50Ca20 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/1005_S1_PID_50Ca20_2401.root", "READ");
    // TCutG *S1cut50Ca20 = (TCutG *)graphFileS1_50Ca20->Get("s1pidCut2401_1005");

    // TFile *graphFileFE9_50Ca20 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_50Ca20.root", "READ");
    // TCutG *FE9cut50Ca20 = (TCutG *)graphFileFE9_50Ca20->Get("FE9pidCut_50Ca_1053_2024");
    // TFile *graphFileS1_50Ca20 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_50Ca20.root", "READ");
    // TCutG *S1cut50Ca20 = (TCutG *)graphFileS1_50Ca20->Get("s1pidCut_50Ca20_1053_2024");

    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_50Ca20.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_50Ca20.cxx");
    TCutG *FE9cut50Ca20 = (TCutG *)gROOT->FindObject("FE9pidCut_50Ca_1053_2024");
    TCutG *S1cut50Ca20 = (TCutG *)gROOT->FindObject("s1pidCut_50Ca20_1053_2024");

    // 51Sc21 Cut Files ---------------------------------------

    // TFile *graphFileFE9_51Sc21 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/1005_FE9_PID_51Sc21_0502.root", "READ");
    // TCutG *FE9cut51Sc21 = (TCutG *)graphFileFE9_51Sc21->Get("FE9pidCut_51Sc_1005");

    // TFile *graphFileS1_51Sc21 = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/1005_S1_PID_51Sc21_0502.root", "READ");
    // TCutG *S1cut51Sc21 = (TCutG *)graphFileS1_51Sc21->Get("s1pidCut_51Sc_1005");
    // --------------------------------------------------------

    // Check the cut files:
    if (!FE9cut50Ca20)
    // if (!FE9cut51Sc21)
    {
        std::cerr << "Error: Could not load the FE9cut50Ca20 graphical cut!" << std::endl;
        return; // Exit with error code
    }
    if (!S1cut50Ca20)
    // if (!S1cut51Sc21)
    {
        std::cerr << "Error: Could not load the S1cut50Ca20 graphical cut!" << std::endl;
        return; // Exit with error code
    }

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
    // TCanvas *c1 = new TCanvas("c1", "IC Canvas", 1500, 850);
    c1->Divide(3, 3);

    TList *hlist_E_Cal = new TList();

    // TH2F *h_IC_E_Cal = new TH2F("IC_E_Cal", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60); // IC Energy Data
    // hlist_E_Cal->Add(h_IC_E_Cal);
    TH2F *h_IC_E_2Dcut = new TH2F("IC_E_2Dcut", "IC Uncalibrated Energy Data; Segment; IC UnCal Energy", 31, -0.5, 30.5, 1000, 0, 200); // IC Uncalibrated Energy Data
    hlist_E_Cal->Add(h_IC_E_2Dcut);

    TH2F *h_IC_E_Cal_2Dcut = new TH2F("IC_E_Cal_2Dcut", "IC Calibrated Energy Data With Fast Beam; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60); // IC Energy Data
    hlist_E_Cal->Add(h_IC_E_Cal_2Dcut);

    TH2F *h_IC_E_FastSlow_Calib = new TH2F("IC_E_FastSlow_Calib", "IC Calibrated Energy Data With Fast-Slow Beam; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60); // IC Energy Data
    hlist_E_Cal->Add(h_IC_E_FastSlow_Calib);

    TH2F *h_IC_E_Cal_UpdatedC_2Dcut = new TH2F("IC_E_Cal_UpdatedC_2Dcut", "IC Calibrated Energy Data Updated Params; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 100); // IC Energy Data
    hlist_E_Cal->Add(h_IC_E_Cal_UpdatedC_2Dcut);

    
    Double_t IC_E[30];
    Double_t IC_E_Cal[30];
    Double_t IC_C[30];
    Double_t IC_E_Cal_UpdatedC[30];
    Double_t IC_E_Cal_UpdatedC_S1X[30];
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
    TH2F *h_S1_XY = new TH2F("S1_XY", "S1_Y_0:S1_X_0; X Position; Y Position", 1000, -200, 200, 1000, -200, 200);
    TH2F *h_S1_XY_cut = new TH2F("S1_XY_cut", "S1_Y_0:S1_X_0; X Position; Y Position", 1000, -200, 200, 1000, -200, 200);

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
    int countBad = 0;

    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        // tree->GetEntry(1545609);
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

        // h_FE9_PID->Fill(PID_T_0, FE9_X_0);
        // h_S1_PID->Fill(S1PID_T_0, S1_X_0);

        h_FE9_E_0->Fill(FE9_E_0);
        h_FE12_E_0->Fill(FE12_E_0);
        h_S1_E_0->Fill(S1_E_0);
        h_S1_X_0->Fill(S1_X_0);
        h_S1_Y_0->Fill(S1_Y_0);
        h_S1_A_0->Fill(S1_A_0);
        h_S1_B_0->Fill(S1_B_0);
        h_S1_XY->Fill(S1_X_0, S1_Y_0);

        // if (fabs(S1_E_0 - 15.50) < 0.5)
        // if (FE9cut->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (S1_E_0 > 15.50 && S1_E_0 < 15.55) && (fabs(S1_Y_0 - 6.115) < 2))
        // if (FE9cut->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.1) && (fabs(S1_X_0 - 11.86) < 2))
        // if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.5) && (fabs(S1_Y_0 - 6.115) < 3) && (fabs(S1_X_0 - 11.86) < 15))
        // if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && not S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.5) && (fabs(S1_Y_0 - 6.115) < 3) && (fabs(S1_X_0 - 11.86) < 15) && (fabs(S1_A_0 < -0.01747) < 0.04))
        // if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.1) && (fabs(S1_Y_0 - (50)) < 5) && (fabs(S1_X_0 - 20) < 15) && (fabs(S1_A_0 - (-0.01747)) < 0.04) && (fabs(S1_B_0 - 0.003003) < 0.012))
        // if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && (fabs(S1_E_0 - 15.5) < 0.1) && (fabs(S1_Y_0 - (6.792)) < 5) && (fabs(S1_X_0 - 2.798) < 15) && (fabs(S1_A_0 - (0.013)) < 0.04) && (fabs(S1_B_0 - 0.0033) < 0.012))
        // if (cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && S1cut->IsInside(S1PID_T_0, S1_X_0) && S1_E_0 >= 0)
        if (S1_E_0 >= 0.1)
        // if (S1_X_0 < 0 && S1_X_0 > -150)
        // if ((fabs(S1_X_0 - (-15)) < 5) && (fabs(S1_E_0 - 15.5) < 0.5) && S1cut->IsInside(S1PID_T_0, S1_X_0) && cut50Ca20FE9->IsInside(PID_T_0, FE9_X_0) && (fabs(S1_Y_0 - (6.792)) < 1))
        // if (FE9cut51Sc21->IsInside(PID_T_0, FE9_X_0) && S1cut51Sc21->IsInside(S1PID_T_0, S1_X_0))
        // if (FE9cut50Ca20->IsInside(PID_T_0, FE9_X_0) && S1cut50Ca20->IsInside(S1PID_T_0, S1_X_0) && abs(S1_E_0 - 13.5) < 0.2 && abs(S1_Y_0 - 0) < 10)
        // if (FE9cut50Ca20->IsInside(PID_T_0, FE9_X_0) && S1cut50Ca20->IsInside(S1PID_T_0, S1_X_0))
        {
            // h_S1_E_0_cut->Fill(S1_E_0);
            // h_S1_X_0_cut->Fill(S1_X_0);
            // h_S1_Y_0_cut->Fill(S1_Y_0);
            // h_S1_A_0_cut->Fill(S1_A_0);
            // h_S1_B_0_cut->Fill(S1_B_0);

            // cout << "Event: " << i << "  " << S1_E_0 << " " << S1_X_0 << " " << S1_Y_0 << " " << S1_A_0 << " " << S1_B_0 << endl;
            // count++;

            int eventNbr;

            for (int ix = 0; ix < 30; ix++)
            {
                // h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);
                // Gates are applied here
                // if (cut->IsInside(PID_T_0, FE9_X_0) && fabs(Beam_T_0 - 245.0) < 1)
                //     h_IC_E_Cal_2Dcut->Fill(ix, IC_E_Cal[ix]);

                // if (((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3] + IC_E[4] + IC_E[5] + IC_E[6] + IC_E[7]) / 4) < 10)
                // if (((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3] + IC_E[4] + IC_E[5] + IC_E[6] + IC_E[7]) / 4) > 10)
                // if (((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3] + IC_E[4] + IC_E[5] + IC_E[6] + IC_E[7]) / 8) < 10)
                // if (((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3] + IC_E[4] + IC_E[5] + IC_E[6] + IC_E[7]) / 8) > 10)
                // if (((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3] + IC_E[4] + IC_E[5] + IC_E[6] + IC_E[7]) / 8.0) > 7)
                // if (((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3] + IC_E[4] + IC_E[5] + IC_E[6] + IC_E[7]) / 8.0) < 7)
                // if (((IC_E[0] + IC_E[1] + IC_E[2] + IC_E[3] + IC_E[4] + IC_E[5] + IC_E[6] + IC_E[7] + IC_E[8] + IC_E[9]) / 10) > 6.5)
                if (((IC_E_Cal[0] + IC_E_Cal[1] + IC_E_Cal[2] + IC_E_Cal[3] + IC_E_Cal[4] + IC_E_Cal[5] + IC_E_Cal[6] + IC_E_Cal[7]) / 8.0) > 6) // for 2024 Data
                {
                    h_IC_E_2Dcut->Fill(ix, IC_E[ix]);
                    // h_IC_E_Cal_2Dcut->Fill(ix, IC_E_Cal[ix]);
                    h_IC_E_Cal_2Dcut->Fill(ix, IC_E[ix] * scalingFactor2024[ix]);

                    // h_IC_E_FastSlow_Calib->Fill(ix, IC_E[ix] * scalingFactor2024[ix] * slowBeamScalingFactor2024[ix]);
                    h_IC_E_FastSlow_Calib->Fill(ix, IC_E[ix] * scalingFactor2024[ix] * slowBeamScalingFactor2024v2[ix]);
                    // count++;

                    // // Using Updated 2022 Params
                    // IC_E_Cal_UpdatedC[ix] = (IC_C[ix] * IC_V_M[ix] + IC_V_N[ix]) * IC_E_M[ix] + IC_E_N[ix];
                    // h_IC_E_Cal_UpdatedC_2Dcut->Fill(ix, IC_E_Cal_UpdatedC[ix]);

                    eventNbr = i;
                }
            }

            if (i == eventNbr)
            {
                tree->GetEntry(eventNbr);
                // tree->GetEntry(1545609);

                // if (S1_E_0 >= 0 && S1_X_0 >= 0)
                // if (S1_E_0 >= 0 && S1_X_0 >= -150)
                if (S1_E_0 >= 0.1)
                {
                    cout << "Good Event: " << eventNbr << "  " << S1_E_0 << " " << S1_X_0 << " " << S1_Y_0 << " " << S1_A_0 << " " << S1_B_0 << endl;
                    count++;

                    h_FE9_PID->Fill(PID_T_0, FE9_X_0);
                    h_S1_PID->Fill(S1PID_T_0, S1_X_0);
                    // h_S1_XY->Fill(S1_X_0, S1_Y_0);
                    h_S1_XY_cut->Fill(S1_X_0, S1_Y_0);

                    h_S1_E_0_cut->Fill(S1_E_0);
                    h_S1_X_0_cut->Fill(S1_X_0);
                    h_S1_Y_0_cut->Fill(S1_Y_0);
                    h_S1_A_0_cut->Fill(S1_A_0);
                    h_S1_B_0_cut->Fill(S1_B_0);
                }
                else
                {
                    countBad++;
                }

                // h_FE9_PID->Fill(PID_T_0, FE9_X_0);
                // h_S1_PID->Fill(S1PID_T_0, S1_X_0);

                // h_S1_E_0_cut->Fill(S1_E_0);
                // h_S1_X_0_cut->Fill(S1_X_0);
                // h_S1_Y_0_cut->Fill(S1_Y_0);
                // h_S1_A_0_cut->Fill(S1_A_0);
                // h_S1_B_0_cut->Fill(S1_B_0);

                // cout << "Bad Event: " << eventNbr << "  " << S1_E_0 << " " << S1_X_0 << " " << S1_Y_0 << " " << S1_A_0 << " " << S1_B_0 << endl;
            }
        }
    }

    cout << "Good Count: " << count << endl;
    cout << "Bad Count: " << countBad << endl;
    // cout << count / 30 << endl;

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

    TGraph *graph2 = new TGraph(n, ICSegment.data(), EnergyLoss16_5MeV.data());
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerColor(kOrange);

    // c1->SetLogy(); // Set log scale for Y-axis

    c1->cd(1);
    h_FE9_PID->Draw("colz");
    FE9cut50Ca20->Draw("SAME");
    FE9cut50Ca20->SetLineColor(kRed);
    FE9cut50Ca20->SetLineWidth(2);
    // FE9cut51Sc21->Draw("SAME");
    // FE9cut51Sc21->SetLineColor(kRed); // Set line color to red
    // FE9cut51Sc21->SetLineWidth(2);    // Set line width to 2

    c1->cd(2);
    h_S1_PID->Draw("colz");
    S1cut50Ca20->Draw("SAME");
    S1cut50Ca20->SetLineColor(kRed);
    S1cut50Ca20->SetLineWidth(2);
    // S1cut51Sc21->Draw("SAME");
    // S1cut51Sc21->SetLineColor(kRed); // Set line color to red
    // S1cut51Sc21->SetLineWidth(2);    // Set line width to 2

    c1->cd(3);
    // h_FE9_E_0->Draw();
    // h_FE12_E_0->Draw();
    h_S1_E_0->Draw();
    h_S1_E_0_cut->SetLineColor(kRed); // Set line color to red
    h_S1_E_0_cut->Draw("SAME");
    // h_S1_XY->Draw("colz");
    // h_S1_XY_cut->SetMarkerColor(kRed); // Set line color to red
    // h_S1_XY_cut->SetMarkerStyle(10);   // Set line width to 10
    // // h_S1_XY_cut->Draw("SAME");
    // h_S1_XY_cut->Draw("colz");

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

    // c1->cd(8);
    // h_IC_E_2Dcut->Draw("colz");
    // // h_IC_E_Cal_2Dcut->Draw("colz");
    // graph->Draw("P SAME");

    // c1->cd(9);
    c1->cd(8);
    h_IC_E_Cal_2Dcut->Draw("colz");
    // h_IC_E_Cal_UpdatedC_2Dcut->Draw("colz");
    graph->Draw("P SAME");
    graph2->Draw("P SAME");

    c1->cd(9);
    h_IC_E_FastSlow_Calib->Draw("colz");
    graph->Draw("P SAME");
    graph2->Draw("P SAME");

   

    // c1->cd(8);
    // h_IC_E_2Dcut->SetLineColor(kBlack);
    // h_IC_E_2Dcut->SetMarkerColor(kBlack);
    // h_IC_E_2Dcut->SetFillColor(0);
    // h_IC_E_2Dcut->Draw("BOX");
    // graph->Draw("P SAME");

    // c1->cd(9);
    // h_IC_E_Cal_2Dcut->SetLineColor(kBlack);
    // h_IC_E_Cal_2Dcut->SetMarkerColor(kBlack);
    // h_IC_E_Cal_2Dcut->SetFillColor(0);
    // h_IC_E_Cal_2Dcut->Draw("BOX");
    // graph->Draw("P SAME");

    c1->Update();
}