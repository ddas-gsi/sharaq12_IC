#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <stdio.h>
#include <string>

#define HEIGHT 512
#define WIDTH 10

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"

// void analyzer(int fEvent = 0, int lEvent = 1)
void analyzer_spline(int RUN_Nbr_int)
{
  // string RUN_Nbr = "1052";
  string RUN_Nbr = std::to_string(RUN_Nbr_int);

  // std::cout << "Copying the root file " << RUN_Nbr << "_Spline_2024.root to the AoQ folder ..." << std::endl;
  // std::string copy_command = "cp /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/" + RUN_Nbr +
  //                            "_Spline_2024.root /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/AoQ/" + RUN_Nbr + "_Spline_2024.root";
  // std::system(copy_command.c_str());

  string inputfile = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/" + RUN_Nbr + "_Spline_2024.root";
  TFile *file = new TFile(inputfile.c_str(), "");
  TTree *tree = (TTree *)file->Get("tree_new");

  // Variables required for the SHARAQ12 analysis in units of [mm], [mrad], [MeV] and [ns].

  double F3_T[10];  // Time at F3.
  double PID_T[10]; // Time of Flight between F3 and FE9.
  double FE9_X[10]; // Horizontal position at FE9.
  double FE9_Y[10]; // Vertical position at FE9.
  double FE9_A[10]; // Horizontal angle at FE9.
  double FE9_B[10]; // Vertical angle at FE9.

  double Beam_T[10]; // Time of Flight between FE9 and FE12.
  double S0_X[10];   // Horizontal position at S0.
  double S0_Y[10];   // Vertical position at S0.
  double S0_A[10];   // Horizontal angle at S0.
  double S0_B[10];   // Vertical angle at S0.
  // double S0_T[10];   // Time at S0.

  double S1PID_T[10]; // Time of Flight between FE12 and S1.
  double S1_X[10];    // Horizontal position at S1.
  double S1_Y[10];    // Vertical position at S1.
  double S1_A[10];    // Horizontal angle at S1.
  double S1_B[10];    // Vertical angle at S1.

  double IC_X[10]; // Horizontal position at IC Entrance.
  double IC_Y[10]; // Vertical position at IC Entrance.
  // double S1_T[10]; // Time at S1.

  double FE9_E[10];  // Energy in MeV/u at FE9.
  double FE12_E[10]; // Energy in MeV/u at FE12.
  double S1_E[10];   // Energy in MeV/u at S1.

  double betaFE9[10];  // Beta at FE9.
  double betaFE12[10]; // Beta at FE12.
  double betaS1[10];   // Beta at S1.

  // for my Analysis

  double splineRange;
  double splineTotalE;
  double splinePeakE;
  double x_scaling;
  double y_scaling;
  double chi2;
  double avg3padE;
  double a3EBetaS1;
  double sPeakES1X;
  double sRangeBetaS1;

  // double range_expt;
  // double range_simCoarse;
  // double peak_sim;
  // double rangeAll;
  // double peakAll;
  // double rangeCa;
  // double peakCa;
  // double deltaE;
  // double totalE;
  // double peakS1YCa;
  // double peakS1YXCa;
  // double deltaCa_st;
  // double deltaCa_lise;
  // double peakS1YAll;
  // double peakS1YXAll;
  // double deltaAll_st;
  // double deltaAll_lise;
  // double deltaCa_stS1Y;
  // double deltaAll_stS1Y;

  double IC_E[30];     // Energy of each segment of the IC chamber.
  double IC_E_Cal[30]; // Calibrated energy of each segment of the IC chamber.
  double IC_C[30];

  double TTT_F_E[5]; // Highest TTT front energies.
  double TTT_B_E[5]; // Highest TTT back energies.
  int TTT_F_ID[5];   // Strip IDs of the highest TTT front energies.
  int TTT_B_ID[5];   // Strip IDs of the highest TTT back energies.

  double CsI_E[16]; // Energy of each TTT-backing CsI crystal.
  double CsI_Time;
  double CsI_Energy;

  double YY1_E[5]; // Highest YY1 energies.
  int YY1_ID[5];   // Strip IDs of the highest YY1 energies.
  double YY1_Time;

  double YY1_T_RAW[6];  // RAW TINA2 YY1 time of each sector.
  double TTT_T_RAW;     // RAW TINA2 TTT time.
  double CsI_T_RAW[16]; // RAW TINA2 CsI time of each sector.

  // Functions useful for the analysis.

  vector<double> f_Coordinates(int ID_Front, int ID_Back);                                                              // Provides Coordinates of TTT pixel.
  double f_Theta(double X, double Y, double Z, double X_0, double Y_0, double Z_0);                                     // Calculates the angle between two vectors.
  double f_Phi(double X, double Y);                                                                                     // Azimuthal angle of a point with respect to the beam axis.
  double f_Excitation_Energy(double Proton_Energy, double Scattering_Angle, double Beam_Energy, double Incident_Angle); // Missing mass spectrum calculator.
  double f_list(int a, int b);                                                                                          // Interstrip Cross-Talk correction for TTT back strips.
  double f_AQ_2_S1_Y(int a);                                                                                            // F3-FE9 ToF per event number (1250 sections)).
  double f_AQ_2_S1_X(int a);                                                                                            // F3-FE9 ToF per event number (1250 sections)).

  double f_Energy(double TTT_E, double Theta); // Calculates Total Energy from TTT Energy and angle.

  double list[3][512]; // First index = 0 (1) means vertical (horizontal)strips, second index indicates SRPPAC index starting from FE91 and third index strip number starting from 0.
  for (int i = 7; i < 10; i++)
  {
    for (int j = 0; j < 512; j++)
    {
      list[i - 7][j] = f_list(i, j);
    }
  }

  double AQ_2_S1_Y[319];
  for (int i = 0; i < 319; i++)
  {
    AQ_2_S1_Y[i] = f_AQ_2_S1_Y(i);
  }

  double AQ_2_S1_X[559];
  for (int i = 0; i < 559; i++)
  {
    AQ_2_S1_X[i] = f_AQ_2_S1_X(i);
  }

  // Tree reading.
  tree->SetBranchStatus("*", 0);

  for (int k = 0; k < 7; k++)
  {
    string F3_T_SetBranch = "F3_T_" + to_string(k);
    tree->SetBranchStatus(F3_T_SetBranch.c_str(), 1);
    tree->SetBranchAddress(F3_T_SetBranch.c_str(), &F3_T[k]);

    string PID_T_SetBranch = "PID_T_" + to_string(k);
    tree->SetBranchStatus(PID_T_SetBranch.c_str(), 1);
    tree->SetBranchAddress(PID_T_SetBranch.c_str(), &PID_T[k]);

    string FE9_X_SetBranch = "FE9_X_" + to_string(k);
    tree->SetBranchStatus(FE9_X_SetBranch.c_str(), 1);
    tree->SetBranchAddress(FE9_X_SetBranch.c_str(), &FE9_X[k]);

    string FE9_Y_SetBranch = "FE9_Y_" + to_string(k);
    tree->SetBranchStatus(FE9_Y_SetBranch.c_str(), 1);
    tree->SetBranchAddress(FE9_Y_SetBranch.c_str(), &FE9_Y[k]);

    string FE9_A_SetBranch = "FE9_A_" + to_string(k);
    tree->SetBranchStatus(FE9_A_SetBranch.c_str(), 1);
    tree->SetBranchAddress(FE9_A_SetBranch.c_str(), &FE9_A[k]);

    string FE9_B_SetBranch = "FE9_B_" + to_string(k);
    tree->SetBranchStatus(FE9_B_SetBranch.c_str(), 1);
    tree->SetBranchAddress(FE9_B_SetBranch.c_str(), &FE9_B[k]);
    string Beam_T_SetBranch = "Beam_T_" + to_string(k);
    tree->SetBranchStatus(Beam_T_SetBranch.c_str(), 1);
    tree->SetBranchAddress(Beam_T_SetBranch.c_str(), &Beam_T[k]);

    string S0_X_SetBranch = "S0_X_" + to_string(k);
    tree->SetBranchStatus(S0_X_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S0_X_SetBranch.c_str(), &S0_X[k]);

    string S0_Y_SetBranch = "S0_Y_" + to_string(k);
    tree->SetBranchStatus(S0_Y_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S0_Y_SetBranch.c_str(), &S0_Y[k]);

    string S0_A_SetBranch = "S0_A_" + to_string(k);
    tree->SetBranchStatus(S0_A_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S0_A_SetBranch.c_str(), &S0_A[k]);

    string S0_B_SetBranch = "S0_B_" + to_string(k);
    tree->SetBranchStatus(S0_B_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S0_B_SetBranch.c_str(), &S0_B[k]);

    string S1PID_T_SetBranch = "S1PID_T_" + to_string(k);
    tree->SetBranchStatus(S1PID_T_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S1PID_T_SetBranch.c_str(), &S1PID_T[k]);

    string S1_X_SetBranch = "S1_X_" + to_string(k);
    tree->SetBranchStatus(S1_X_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S1_X_SetBranch.c_str(), &S1_X[k]);

    string S1_Y_SetBranch = "S1_Y_" + to_string(k);
    tree->SetBranchStatus(S1_Y_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S1_Y_SetBranch.c_str(), &S1_Y[k]);

    string S1_A_SetBranch = "S1_A_" + to_string(k);
    tree->SetBranchStatus(S1_A_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S1_A_SetBranch.c_str(), &S1_A[k]);

    string S1_B_SetBranch = "S1_B_" + to_string(k);
    tree->SetBranchStatus(S1_B_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S1_B_SetBranch.c_str(), &S1_B[k]);

    // for my analysis
    string FE9_E_SetBranch = "FE9_E_" + to_string(k);
    tree->SetBranchStatus(FE9_E_SetBranch.c_str(), 1);
    tree->SetBranchAddress(FE9_E_SetBranch.c_str(), &FE9_E[k]);

    string FE12_E_SetBranch = "FE12_E_" + to_string(k);
    tree->SetBranchStatus(FE12_E_SetBranch.c_str(), 1);
    tree->SetBranchAddress(FE12_E_SetBranch.c_str(), &FE12_E[k]);

    string S1_E_SetBranch = "S1_E_" + to_string(k);
    tree->SetBranchStatus(S1_E_SetBranch.c_str(), 1);
    tree->SetBranchAddress(S1_E_SetBranch.c_str(), &S1_E[k]);

    string betaFE9_SetBranch = "betaFE9_" + to_string(k);
    tree->SetBranchStatus(betaFE9_SetBranch.c_str(), 1);
    tree->SetBranchAddress(betaFE9_SetBranch.c_str(), &betaFE9[k]);

    string betaFE12_SetBranch = "betaFE12_" + to_string(k);
    tree->SetBranchStatus(betaFE12_SetBranch.c_str(), 1);
    tree->SetBranchAddress(betaFE12_SetBranch.c_str(), &betaFE12[k]);

    string betaS1_SetBranch = "betaS1_" + to_string(k);
    tree->SetBranchStatus(betaS1_SetBranch.c_str(), 1);
    tree->SetBranchAddress(betaS1_SetBranch.c_str(), &betaS1[k]);
  }

  for (int k = 0; k < 30; k++)
  {
    string IC_E_SetBranch = "IC_E_" + to_string(k);
    tree->SetBranchStatus(IC_E_SetBranch.c_str(), 1);
    tree->SetBranchAddress(IC_E_SetBranch.c_str(), &IC_E[k]);

    string IC_E_Cal_SetBranch = "IC_E_Cal_" + to_string(k);
    tree->SetBranchStatus(IC_E_Cal_SetBranch.c_str(), 1);
    tree->SetBranchAddress(IC_E_Cal_SetBranch.c_str(), &IC_E_Cal[k]);

    string IC_C_SetBranch = "IC_C_" + to_string(k);
    tree->SetBranchStatus(IC_C_SetBranch.c_str(), 1);
    tree->SetBranchAddress(IC_C_SetBranch.c_str(), &IC_C[k]);
  }

  // for my Analysis

  tree->SetBranchStatus("splineRange", 1);
  tree->SetBranchAddress("splineRange", &splineRange);
  tree->SetBranchStatus("splineTotalE", 1);
  tree->SetBranchAddress("splineTotalE", &splineTotalE);
  tree->SetBranchStatus("splinePeakE", 1);
  tree->SetBranchAddress("splinePeakE", &splinePeakE);
  tree->SetBranchStatus("x_scaling", 1);
  tree->SetBranchAddress("x_scaling", &x_scaling);
  tree->SetBranchStatus("y_scaling", 1);
  tree->SetBranchAddress("y_scaling", &y_scaling);
  tree->SetBranchStatus("chi2", 1);
  tree->SetBranchAddress("chi2", &chi2);
  tree->SetBranchStatus("avg3padE", 1);
  tree->SetBranchAddress("avg3padE", &avg3padE);
  tree->SetBranchStatus("a3EBetaS1", 1);
  tree->SetBranchAddress("a3EBetaS1", &a3EBetaS1);
  tree->SetBranchStatus("sPeakES1X", 1);
  tree->SetBranchAddress("sPeakES1X", &sPeakES1X);
  tree->SetBranchStatus("sRangeBetaS1", 1);
  tree->SetBranchAddress("sRangeBetaS1", &sRangeBetaS1);

  // tree->SetBranchStatus("range_expt", 1);
  // tree->SetBranchAddress("range_expt", &range_expt);

  // tree->SetBranchStatus("range_simCoarse", 1);
  // tree->SetBranchAddress("range_simCoarse", &range_simCoarse);

  // tree->SetBranchStatus("peak_sim", 1);
  // tree->SetBranchAddress("peak_sim", &peak_sim);

  // tree->SetBranchStatus("rangeAll", 1);
  // tree->SetBranchAddress("rangeAll", &rangeAll);

  // tree->SetBranchStatus("peakAll", 1);
  // tree->SetBranchAddress("peakAll", &peakAll);

  // tree->SetBranchStatus("rangeCa", 1);
  // tree->SetBranchAddress("rangeCa", &rangeCa);

  // tree->SetBranchStatus("peakCa", 1);
  // tree->SetBranchAddress("peakCa", &peakCa);

  // tree->SetBranchStatus("deltaE", 1);
  // tree->SetBranchAddress("deltaE", &deltaE);

  // tree->SetBranchStatus("totalE", 1);
  // tree->SetBranchAddress("totalE", &totalE);

  // tree->SetBranchStatus("peakS1YCa", 1);
  // tree->SetBranchAddress("peakS1YCa", &peakS1YCa);

  // tree->SetBranchStatus("peakS1YXCa", 1);
  // tree->SetBranchAddress("peakS1YXCa", &peakS1YXCa);

  // tree->SetBranchStatus("deltaCa_st", 1);
  // tree->SetBranchAddress("deltaCa_st", &deltaCa_st);

  // tree->SetBranchStatus("deltaCa_lise", 1);
  // tree->SetBranchAddress("deltaCa_lise", &deltaCa_lise);

  // tree->SetBranchStatus("peakS1YAll", 1);
  // tree->SetBranchAddress("peakS1YAll", &peakS1YAll);

  // tree->SetBranchStatus("peakS1YXAll", 1);
  // tree->SetBranchAddress("peakS1YXAll", &peakS1YXAll);

  // tree->SetBranchStatus("deltaAll_st", 1);
  // tree->SetBranchAddress("deltaAll_st", &deltaAll_st);

  // tree->SetBranchStatus("deltaAll_lise", 1);
  // tree->SetBranchAddress("deltaAll_lise", &deltaAll_lise);

  // tree->SetBranchStatus("deltaCa_stS1Y", 1);
  // tree->SetBranchAddress("deltaCa_stS1Y", &deltaCa_stS1Y);

  // tree->SetBranchStatus("deltaAll_stS1Y", 1);
  // tree->SetBranchAddress("deltaAll_stS1Y", &deltaAll_stS1Y);

  // Variables calculated during the analysis.

  double PI = 3.14159265358979323846;

  double M_49Ca = 45601.606;  // Mass of 49Ca (MeV).
  double M_50Ca = 46535.0915; // Mass of 50Ca (MeV).
  double M_51Sc = 47462.924;  // Mass of 51Sc (MeV).
  double M_49K = 45613.5757;  // Mass of 49K (MeV).

  // F3-FE9 PID Gates.
  double FS_50Ca_X_A = 1083.2;
  double FS_50Ca_X_B = 1098.4;
  double FS_50Ca_X_C = 1087.4;
  double FS_50Ca_X_D = 1073.1;
  double FS_50Ca_Y_A = 37.6;
  double FS_50Ca_Y_B = -21.50;
  double FS_50Ca_Y_C = -45.18;
  double FS_50Ca_Y_D = 26.07;

  double FS_50Ca_M_AB = (FS_50Ca_Y_A - FS_50Ca_Y_B) / (FS_50Ca_X_A - FS_50Ca_X_B);
  double FS_50Ca_M_BC = (FS_50Ca_Y_B - FS_50Ca_Y_C) / (FS_50Ca_X_B - FS_50Ca_X_C);
  double FS_50Ca_M_CD = (FS_50Ca_Y_C - FS_50Ca_Y_D) / (FS_50Ca_X_C - FS_50Ca_X_D);
  double FS_50Ca_M_DA = (FS_50Ca_Y_D - FS_50Ca_Y_A) / (FS_50Ca_X_D - FS_50Ca_X_A);
  double FS_50Ca_N_AB = FS_50Ca_Y_A - FS_50Ca_M_AB * FS_50Ca_X_A;
  double FS_50Ca_N_BC = FS_50Ca_Y_B - FS_50Ca_M_BC * FS_50Ca_X_B;
  double FS_50Ca_N_CD = FS_50Ca_Y_C - FS_50Ca_M_CD * FS_50Ca_X_C;
  double FS_50Ca_N_DA = FS_50Ca_Y_D - FS_50Ca_M_DA * FS_50Ca_X_D;

  // S0 Beam Energy Calculation.
  double beta;
  double gamma;

  double Z_ToF = 14871.68;
  double E_Beam;
  double E_Beam_C[10];
  double E_Beam_2[10];
  double E_Loss;
  // double Beam_param[9] = {4.1620*pow(10,-23), -5.4000*pow(10,-19), 2.9192*pow(10,-15), -8.5111*pow(10,-12), 1.4455*pow(10,-8), -1.4362*pow(10,-5), 7.8612*pow(10,-3), -9.3805*pow(10,-1), 1.2268*pow(10,1)};
  // double Beam_param_2[9] = {9.5015*pow(10,-24), -1.0846*pow(10,-19), 5.0308*pow(10,-16), -1.2190*pow(10,-12), 1.6543*pow(10,-9), -1.2658*pow(10,-6), 5.4521*pow(10,-4), 8.6821*pow(10,-1), -6.5966*pow(10,-1)};
  double Beam_param[6] = {1.0449 * pow(10, -6), -1.0025 * pow(10, -4), 3.9070 * pow(10, -3), -7.9621 * pow(10, -2), 1.9126 * pow(10, 0), -6.0741 * pow(10, 0)};
  double Beam_param_2[6] = {2.1635 * pow(10, -7), -2.1897 * pow(10, -6), 9.1640 * pow(10, -5), -2.0621 * pow(10, -3), 1.0273 * pow(10, 0), -2.2604 * pow(10, -1)};
  double Light_param[9] = {-1.0457 * pow(10, -8), 8.5944 * pow(10, -7), -2.8971 * pow(10, -5), 5.1546 * pow(10, -4), -5.1807 * pow(10, -3), 2.9085 * pow(10, -2), -8.3397 * pow(10, -2), 1.093 * pow(10, 0), 8.4547 * pow(10, -3)};
  double S0_X_Shift = 1.;
  double S0_Y_Shift = 6.;
  double S0_Gate = 20.;

  double Z_ToF_S1 = 9172.;
  double E_Beam_S1;
  double E_Beam_C_S1;
  double E_Loss_S1;
  double Beam_param_S1[9] = {4.1620 * pow(10, -23), -5.4000 * pow(10, -19), 2.9192 * pow(10, -15), -8.5111 * pow(10, -12), 1.4455 * pow(10, -8), -1.4362 * pow(10, -5), 7.8612 * pow(10, -3), -9.3805 * pow(10, -1), 1.2268 * pow(10, 1)};

  // S1 Ionization Chamber.
  double IC_segment = 25.233; // Length of each individual segment [mm].
  double IC_Brho;             // Brho.

  // For My analysis.
  double IC_Brho_array[10];            // IC_Brho array for branch fillup
  double beta_staraq_array[10];        // for beta_sharaq branch fillup
  double beta_sharaq_f_oedo_array[10]; // for beta_sharaq_f_oedo branch fillup. Beta SHARAQ from OEDO beam Energy

  double AQ[10];   // A/Q at IC.
  double AQ_2[10]; // AQ at S1 (2).

  double IC_Brho_zero = 1.3671;
  double IC_X0 = -0.4;
  double IC_delta = 2222.;
  double IC_delta_a = 50.;

  /***************************
          Output Tree
  ***************************/

  // We create a new File with a  smaller Tree in which we fill the good events only.
  // string output_new = "/home/jupiter/test/work/artemis-oedo/output/Analysis/testanalyzed_" + to_string(fEvent / 1000000) + ".root";
  // TFile *ofile_new = new TFile(output_new.c_str(), "recreate");
  // TTree *tree_new = new TTree("tree_new","Good events sharaq12");

  // TFile *ofile_new = new TFile("/home/jupiter/test/work/artemis-oedo/output/Analysis/testanalyzed_try2.root","recreate");
  // TFile *ofile_new = new TFile("/Users/cfergon12/geant4_workdir/Files/F3DSgood.root","recreate");
  // TFile *ofile_new = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/1056_Analysis_2024.root", "recreate");

  string outputfile = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/AoQ/" + RUN_Nbr + "_Spline_2024.root";
  TFile *ofile_new = new TFile(outputfile.c_str(), "recreate");
  TTree *tree_new = new TTree("tree_new", "Good events sharaq12");

  TList *hlist = new TList();

  TH2F *h_PID_T_Event = new TH2F("PID_T_Event", "PID T [ns] vs Event.", 5000, 0, 100000000, 5000, 900, 1200);
  hlist->Add(h_PID_T_Event);

  TH2F *h_PID_T_Event_C = new TH2F("PID_T_Event_C", "PID T (C) [ns] vs Event.", 5000, 0, 100000000, 5000, 900, 1200);
  hlist->Add(h_PID_T_Event_C);

  TH2F *h_Beam_T_Event = new TH2F("Beam_T_Event", "Beam T [ns] vs Event.", 1000, 0, 100000000, 100, 240, 290);
  hlist->Add(h_Beam_T_Event);

  TH2F *h_Beam_T_Event_C = new TH2F("Beam_T_Event_C", "Beam T (C) [ns] vs Event.", 1000, 0, 100000000, 100, 240, 290);
  hlist->Add(h_Beam_T_Event_C);

  TH2F *h_S1PID_T_Event = new TH2F("S1PID_T_Event", "S1PID T [ns] vs Event.", 5000, 0, 100000000, 5000, 580, 640);
  hlist->Add(h_S1PID_T_Event);

  TH2F *h_S1PID_T_Event_C = new TH2F("S1PID_T_Event_C", "S1PID T (C) [ns] vs Event.", 5000, 0, 100000000, 5000, 580, 640);
  hlist->Add(h_S1PID_T_Event_C);

  TH2F *h_AQ_Event = new TH2F("AQ_Event", "AQ vs Event.", 5000, 0, 100000000, 5000, 2.3, 3.);
  hlist->Add(h_AQ_Event);

  TH2F *h_AQ_Event_C = new TH2F("AQ_Event_C", "AQ (C) [ns] vs Event.", 5000, 0, 100000000, 5000, 2.3, 3.);
  hlist->Add(h_AQ_Event_C);

  // Check for beam energy.
  TH2F *h_TTT = new TH2F("TTT", "AQ (C) [ns] vs Event.", 30, -0.5, 29.5, 4, -0.5, 3.5);
  hlist->Add(h_TTT);

  TH2F *h_YY1 = new TH2F("YY1", "AQ (C) [ns] vs Event.", 30, -0.5, 29.5, 4, -0.5, 3.5);
  hlist->Add(h_YY1);

  TH1F *h_Ex_TTT[51][4];
  for (int i = 0; i < 51; i++)
  {
    h_Ex_TTT[i][0] = new TH1F(Form("Ex_TTT_beta_%d", i), Form("Ex_TTT_beta_%d", i), 140, -1, 6);
    hlist->Add(h_Ex_TTT[i][0]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_TTT[i][1] = new TH1F(Form("Ex_TTT_beta_%d_X", i), Form("Ex_TTT_beta_%d_X", i), 140, -1, 6);
    hlist->Add(h_Ex_TTT[i][1]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_TTT[i][2] = new TH1F(Form("Ex_TTT_beta_%d_Y", i), Form("Ex_TTT_beta_%d_Y", i), 140, -1, 6);
    hlist->Add(h_Ex_TTT[i][2]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_TTT[i][3] = new TH1F(Form("Ex_TTT_beta_%d_XY", i), Form("Ex_TTT_beta_%d_XY", i), 140, -1, 6);
    hlist->Add(h_Ex_TTT[i][3]);
  }

  TH1F *h_Ex_YY1[51][4];
  for (int i = 0; i < 51; i++)
  {
    h_Ex_YY1[i][0] = new TH1F(Form("Ex_YY1_beta_%d", i), Form("Ex_YY1_beta_%d", i), 140, -1, 6);
    hlist->Add(h_Ex_YY1[i][0]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_YY1[i][1] = new TH1F(Form("Ex_YY1_beta_%d_X", i), Form("Ex_YY1_beta_%d_X", i), 140, -1, 6);
    hlist->Add(h_Ex_YY1[i][1]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_YY1[i][2] = new TH1F(Form("Ex_YY1_beta_%d_Y", i), Form("Ex_YY1_beta_%d_Y", i), 140, -1, 6);
    hlist->Add(h_Ex_YY1[i][2]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_YY1[i][3] = new TH1F(Form("Ex_YY1_beta_%d_XY", i), Form("Ex_YY1_beta_%d_XY", i), 140, -1, 6);
    hlist->Add(h_Ex_YY1[i][3]);
  }

  TH1F *h_Ex_TINA[51][4];
  for (int i = 0; i < 51; i++)
  {
    h_Ex_TINA[i][0] = new TH1F(Form("Ex_TINA_beta_%d", i), Form("Ex_TINA_beta_%d", i), 140, -1, 6);
    hlist->Add(h_Ex_TINA[i][0]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_TINA[i][1] = new TH1F(Form("Ex_TINA_beta_%d_X", i), Form("Ex_TINA_beta_%d_X", i), 140, -1, 6);
    hlist->Add(h_Ex_TINA[i][1]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_TINA[i][2] = new TH1F(Form("Ex_TINA_beta_%d_Y", i), Form("Ex_TINA_beta_%d_Y", i), 140, -1, 6);
    hlist->Add(h_Ex_TINA[i][2]);
  }

  for (int i = 0; i < 51; i++)
  {
    h_Ex_TINA[i][3] = new TH1F(Form("Ex_TINA_beta_%d_XY", i), Form("Ex_TINA_beta_%d_XY", i), 140, -1, 6);
    hlist->Add(h_Ex_TINA[i][3]);
  }

  TH1F *h_Ex_TINA_Shift[31][31];
  TH1F *h_Ex_TTT_Shift[31][31];
  TH1F *h_Ex_YY1_Shift[31][31];
  TH1F *h_Ex_TINA_Shift_X[31][31];
  TH1F *h_Ex_TTT_Shift_X[31][31];
  TH1F *h_Ex_YY1_Shift_X[31][31];
  TH1F *h_Ex_TINA_Shift_Y[31][31];
  TH1F *h_Ex_TTT_Shift_Y[31][31];
  TH1F *h_Ex_YY1_Shift_Y[31][31];
  TH1F *h_Ex_TINA_Shift_XY[31][31];
  TH1F *h_Ex_TTT_Shift_XY[31][31];
  TH1F *h_Ex_YY1_Shift_XY[31][31];
  for (int i = 0; i < 31; i++)
  {
    for (int j = 0; j < 31; j++)
    {
      h_Ex_TINA_Shift[i][j] = new TH1F(Form("Ex_TINA_Shift_%d_%d", i, j), Form("Ex_TINA_Shift X = %d, Y = %d", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_TINA_Shift[i][j]);
      h_Ex_TTT_Shift[i][j] = new TH1F(Form("Ex_TTT_Shift_%d_%d", i, j), Form("Ex_TTT_Shift X = %d, Y = %d", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_TTT_Shift[i][j]);
      h_Ex_YY1_Shift[i][j] = new TH1F(Form("Ex_YY1_Shift_%d_%d", i, j), Form("Ex_YY1_Shift X = %d, Y = %d", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_YY1_Shift[i][j]);

      h_Ex_TINA_Shift_X[i][j] = new TH1F(Form("Ex_TINA_Shift_%d_%d_X", i, j), Form("Ex_TINA_Shift X = %d, Y = %d (X Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_TINA_Shift_X[i][j]);
      h_Ex_TTT_Shift_X[i][j] = new TH1F(Form("Ex_TTT_Shift_%d_%d_X", i, j), Form("Ex_TTT_Shift X = %d, Y = %d (X Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_TTT_Shift_X[i][j]);
      h_Ex_YY1_Shift_X[i][j] = new TH1F(Form("Ex_YY1_Shift_%d_%d_X", i, j), Form("Ex_YY1_Shift X = %d, Y = %d (X Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_YY1_Shift_X[i][j]);

      h_Ex_TINA_Shift_Y[i][j] = new TH1F(Form("Ex_TINA_Shift_%d_%d_Y", i, j), Form("Ex_TINA_Shift X = %d, Y = %d (Y Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_TINA_Shift_Y[i][j]);
      h_Ex_TTT_Shift_Y[i][j] = new TH1F(Form("Ex_TTT_Shift_%d_%d_Y", i, j), Form("Ex_TTT_Shift X = %d, Y = %d (Y Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_TTT_Shift_Y[i][j]);
      h_Ex_YY1_Shift_Y[i][j] = new TH1F(Form("Ex_YY1_Shift_%d_%d_Y", i, j), Form("Ex_YY1_Shift X = %d, Y = %d (Y Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_YY1_Shift_Y[i][j]);

      h_Ex_TINA_Shift_XY[i][j] = new TH1F(Form("Ex_TINA_Shift_%d_%d_XY", i, j), Form("Ex_TINA_Shift X = %d, Y = %d (XY Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_TINA_Shift_XY[i][j]);
      h_Ex_TTT_Shift_XY[i][j] = new TH1F(Form("Ex_TTT_Shift_%d_%d_XY", i, j), Form("Ex_TTT_Shift X = %d, Y = %d (XY Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_TTT_Shift_XY[i][j]);
      h_Ex_YY1_Shift_XY[i][j] = new TH1F(Form("Ex_YY1_Shift_%d_%d_XY", i, j), Form("Ex_YY1_Shift X = %d, Y = %d (XY Mirror)", i - 15, j - 15), 140, -1, 6);
      hlist->Add(h_Ex_YY1_Shift_XY[i][j]);
    }
  }

  TH2F *h_Ex = new TH2F("Ex", "Beam Shift [MeV/nucleon/10] vs Mirror.", 4, -0.5, 3.5, 51, -35.5, 15.5);
  hlist->Add(h_Ex);

  TH2F *h_TTT_Ex = new TH2F("TTT_Ex", "Beam Shift [MeV/nucleon/10] vs Mirror. (TTT)", 4, -0.5, 3.5, 51, -35.5, 15.5);
  hlist->Add(h_TTT_Ex);

  TH2F *h_YY1_Ex = new TH2F("YY1_Ex", "Beam Shift [MeV/nucleon/10] vs Mirror. (YY1)", 4, -0.5, 3.5, 51, -35.5, 15.5);
  hlist->Add(h_YY1_Ex);

  TH2F *h_Ex_Shift = new TH2F("Ex_Shift", "Delta X vs Delta Y.", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_Ex_Shift);

  TH2F *h_TTT_Ex_Shift = new TH2F("TTT_Ex_Shift", "Delta X vs Delta Y (TTT).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_TTT_Ex_Shift);

  TH2F *h_YY1_Ex_Shift = new TH2F("YY1_Ex_Shift", "Delta X vs Delta Y (YY1).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_YY1_Ex_Shift);

  TH2F *h_Ex_Shift_X = new TH2F("Ex_Shift_X", "Delta X vs Delta Y (Mirror X).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_Ex_Shift_X);

  TH2F *h_TTT_Ex_Shift_X = new TH2F("TTT_Ex_Shift_X", "Delta X vs Delta Y (TTT, Mirror X).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_TTT_Ex_Shift_X);

  TH2F *h_YY1_Ex_Shift_X = new TH2F("YY1_Ex_Shift_X", "Delta X vs Delta Y (YY1, Mirror X).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_YY1_Ex_Shift_X);

  TH2F *h_Ex_Shift_Y = new TH2F("Ex_Shift_Y", "Delta X vs Delta Y (Mirror Y).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_Ex_Shift_Y);

  TH2F *h_TTT_Ex_Shift_Y = new TH2F("TTT_Ex_Shift_Y", "Delta X vs Delta Y (TTT, Mirror Y).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_TTT_Ex_Shift_Y);

  TH2F *h_YY1_Ex_Shift_Y = new TH2F("YY1_Ex_Shift_Y", "Delta X vs Delta Y (YY1, Mirror Y).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_YY1_Ex_Shift_Y);

  TH2F *h_Ex_Shift_XY = new TH2F("Ex_Shift_XY", "Delta X vs Delta Y (Mirror XY).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_Ex_Shift_XY);

  TH2F *h_TTT_Ex_Shift_XY = new TH2F("TTT_Ex_Shift_XY", "Delta X vs Delta Y (TTT, Mirror XY).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_TTT_Ex_Shift_XY);

  TH2F *h_YY1_Ex_Shift_XY = new TH2F("YY1_Ex_Shift_XY", "Delta X vs Delta Y (YY1, Mirror XY).", 31, -15.5, 15.5, 31, -15.5, 15.5);
  hlist->Add(h_YY1_Ex_Shift_XY);

  TH2F *h_TTT_CsI[16];
  for (int i = 0; i < 16; i++)
  {
    h_TTT_CsI[i] = new TH2F(Form("TTT_CsI_%d", i), Form("Crystal_%d", i), 150, -5, 25, 300, -5, 25);
    hlist->Add(h_TTT_CsI[i]);
  }

  TH2F *h_TTT_CsI_all = new TH2F("TTT_CsI_all", "All Crystals.", 150, -5, 25, 300, -5, 25);
  hlist->Add(h_TTT_CsI_all);

  // TH2F* h_S0 = new TH2F("S0", "S0 Y [mm] vs S0 X [mm].", 101, -25.25, 25.25, 101, -25.25, 25.25);
  // hlist->Add(h_S0);

  // TH2F* h_S0_S1Gate = new TH2F("S0_S1Gate", "S0 Y [mm] vs S0 X [mm] (S1 Gate).", 101, -25.25, 25.25, 101, -25.25, 25.25);
  // hlist->Add(h_S0_S1Gate);

  // TH2F* h_Transmission = new TH2F("Transmission", "S0 Y [mm] vs S0 X [mm] (Transmission).", 101, -25.25, 25.25, 101, -25.25, 25.25);
  // hlist->Add(h_Transmission);

  TH2F *h_CsI_Gate = new TH2F("CsI_Gate", "TTT vs CsI Energy [MeV] (Gate on CsI Banana)", 2000, -3.5, 32.5, 2000, -0.5, 14.5);
  hlist->Add(h_CsI_Gate);

  // We declare the variables que want to fill in the new Tree.
  for (int k = 0; k < 7; k++)
  {
    string F3_T_name = "F3_T_" + to_string(k);
    tree_new->Branch(F3_T_name.c_str(), &F3_T[k], 320000);

    string PID_T_name = "PID_T_" + to_string(k);
    tree_new->Branch(PID_T_name.c_str(), &PID_T[k], 320000);

    string FE9_X_name = "FE9_X_" + to_string(k);
    tree_new->Branch(FE9_X_name.c_str(), &FE9_X[k], 320000);

    string FE9_Y_name = "FE9_Y_" + to_string(k);
    tree_new->Branch(FE9_Y_name.c_str(), &FE9_Y[k], 320000);

    string FE9_A_name = "FE9_A_" + to_string(k);
    tree_new->Branch(FE9_A_name.c_str(), &FE9_A[k], 320000);

    string FE9_B_name = "FE9_B_" + to_string(k);
    tree_new->Branch(FE9_B_name.c_str(), &FE9_B[k], 320000);

    string Beam_T_name = "Beam_T_" + to_string(k);
    tree_new->Branch(Beam_T_name.c_str(), &Beam_T[k], 320000);

    string Beam_E_name = "Beam_E_" + to_string(k);
    tree_new->Branch(Beam_E_name.c_str(), &E_Beam_C[k], 320000);

    string Beam_E_2_name = "Beam_E_2_" + to_string(k);
    tree_new->Branch(Beam_E_2_name.c_str(), &E_Beam_2[k], 320000);

    string S0_X_name = "S0_X_" + to_string(k);
    tree_new->Branch(S0_X_name.c_str(), &S0_X[k], 320000);

    string S0_Y_name = "S0_Y_" + to_string(k);
    tree_new->Branch(S0_Y_name.c_str(), &S0_Y[k], 320000);

    string S0_A_name = "S0_A_" + to_string(k);
    tree_new->Branch(S0_A_name.c_str(), &S0_A[k], 320000);

    string S0_B_name = "S0_B_" + to_string(k);
    tree_new->Branch(S0_B_name.c_str(), &S0_B[k], 320000);

    string S1PID_T_name = "S1PID_T_" + to_string(k);
    tree_new->Branch(S1PID_T_name.c_str(), &S1PID_T[k], 320000);

    string S1_X_name = "S1_X_" + to_string(k);
    tree_new->Branch(S1_X_name.c_str(), &S1_X[k], 320000);

    string S1_Y_name = "S1_Y_" + to_string(k);
    tree_new->Branch(S1_Y_name.c_str(), &S1_Y[k], 320000);

    string S1_A_name = "S1_A_" + to_string(k);
    tree_new->Branch(S1_A_name.c_str(), &S1_A[k], 320000);

    string S1_B_name = "S1_B_" + to_string(k);
    tree_new->Branch(S1_B_name.c_str(), &S1_B[k], 320000);

    string IC_X_name = "IC_X_" + to_string(k);
    tree_new->Branch(IC_X_name.c_str(), &IC_X[k], 320000);

    string IC_Y_name = "IC_Y_" + to_string(k);
    tree_new->Branch(IC_Y_name.c_str(), &IC_Y[k], 320000);

    string AQ_name = "AQ_" + to_string(k);
    tree_new->Branch(AQ_name.c_str(), &AQ[k], 320000);

    string AQ_2_name = "AQ_2_" + to_string(k);
    tree_new->Branch(AQ_2_name.c_str(), &AQ_2[k], 320000);

    // For My Analysis.
    string IC_Brho_name = "IC_Brho_" + to_string(k);
    tree_new->Branch(IC_Brho_name.c_str(), &IC_Brho_array[k], 320000);

    string beta_sharaq_name = "beta_sharaq_" + to_string(k);
    tree_new->Branch(beta_sharaq_name.c_str(), &beta_staraq_array[k], 320000);

    string beta_sharaq_f_oedo_name = "beta_sharaq_f_oedo_" + to_string(k); // beta SHARAQ from OEDO beam Energy
    tree_new->Branch(beta_sharaq_f_oedo_name.c_str(), &beta_sharaq_f_oedo_array[k], 320000);

    string FE9_E_name = "FE9_E_" + to_string(k);
    tree_new->Branch(FE9_E_name.c_str(), &FE9_E[k], 320000);

    string FE12_E_name = "FE12_E_" + to_string(k);
    tree_new->Branch(FE12_E_name.c_str(), &FE12_E[k], 320000);

    string S1_E_name = "S1_E_" + to_string(k);
    tree_new->Branch(S1_E_name.c_str(), &S1_E[k], 320000);

    string betaFE9_name = "betaFE9_" + to_string(k);
    tree_new->Branch(betaFE9_name.c_str(), &betaFE9[k], 320000);

    string betaFE12_name = "betaFE12_" + to_string(k);
    tree_new->Branch(betaFE12_name.c_str(), &betaFE12[k], 320000);

    string betaS1_name = "betaS1_" + to_string(k);
    tree_new->Branch(betaS1_name.c_str(), &betaS1[k], 320000);
  }

  // for (int k = 0; k < 30; k++)
  // {
  //   string IC_C_name = "IC_C_" + to_string(k);
  //   tree_new->Branch(IC_C_name.c_str(), &IC_C[k], 320000);
  // }

  // for (int k = 0; k < 30; k++)
  // {
  //   string IC_E_name = "IC_E_" + to_string(k);
  //   tree_new->Branch(IC_E_name.c_str(), &IC_E[k], 320000);
  // }

  for (int k = 0; k < 30; k++)
  {
    string IC_E_Cal_name = "IC_E_Cal_" + to_string(k);
    tree_new->Branch(IC_E_Cal_name.c_str(), &IC_E_Cal[k], 320000);
  }

  // for my analysis

  tree_new->Branch("splineRange", &splineRange, 320000);
  tree_new->Branch("splineTotalE", &splineTotalE, 320000);
  tree_new->Branch("splinePeakE", &splinePeakE, 320000);
  tree_new->Branch("x_scaling", &x_scaling, 320000);
  tree_new->Branch("y_scaling", &y_scaling, 320000);
  tree_new->Branch("chi2", &chi2, 320000);
  tree_new->Branch("avg3padE", &avg3padE, 320000);
  tree_new->Branch("a3EBetaS1", &a3EBetaS1, 320000);
  tree_new->Branch("sPeakES1X", &sPeakES1X, 320000);
  tree_new->Branch("sRangeBetaS1", &sRangeBetaS1, 320000);

  // tree_new->Branch("range_expt", &range_expt, 320000);
  // tree_new->Branch("range_simCoarse", &range_simCoarse, 320000);
  // tree_new->Branch("peak_sim", &peak_sim, 320000);
  // tree_new->Branch("rangeAll", &rangeAll, 320000);
  // tree_new->Branch("peakAll", &peakAll, 320000);
  // tree_new->Branch("rangeCa", &rangeCa, 320000);
  // tree_new->Branch("peakCa", &peakCa, 320000);
  // tree_new->Branch("deltaE", &deltaE, 320000);
  // tree_new->Branch("totalE", &totalE, 320000);
  // tree_new->Branch("peakS1YCa", &peakS1YCa, 320000);
  // tree_new->Branch("peakS1YXCa", &peakS1YXCa, 320000);
  // tree_new->Branch("deltaCa_st", &deltaCa_st, 320000);
  // tree_new->Branch("deltaCa_lise", &deltaCa_lise, 320000);
  // tree_new->Branch("peakS1YAll", &peakS1YAll, 320000);
  // tree_new->Branch("peakS1YXAll", &peakS1YXAll, 320000);
  // tree_new->Branch("deltaAll_st", &deltaAll_st, 320000);
  // tree_new->Branch("deltaAll_lise", &deltaAll_lise, 320000);
  // tree_new->Branch("deltaCa_stS1Y", &deltaCa_stS1Y, 320000);
  // tree_new->Branch("deltaAll_stS1Y", &deltaAll_stS1Y, 320000);

  // tree_new->Branch("TTT_F_ID", &TTT_Front_ID, 320000);
  // tree_new->Branch("TTT_F_E", &TTT_Front_E, 320000);
  // tree_new->Branch("TTT_B_ID", &TTT_Back_ID, 320000);
  // tree_new->Branch("TTT_B_E", &TTT_Back_E, 320000);
  // tree_new->Branch("YY1_ID", &YY1_Energy_ID, 320000);
  // tree_new->Branch("YY1_E", &YY1_Energy, 320000);
  // tree_new->Branch("CsI_E", &CsI_Energy, 320000);
  // // tree_new->Branch("TINA2_E", &TINA2_E, 320000);

  // tree_new->Branch("TTT_T", &TTT_T_RAW, 320000);
  // tree_new->Branch("YY1_T", &YY1_Time, 320000);
  // tree_new->Branch("CsI_T", &CsI_Time, 320000);

  // tree_new->Branch("TTT_X", &TTT_X, 320000);
  // tree_new->Branch("TTT_Y", &TTT_Y, 320000);
  // tree_new->Branch("TTT_Z", &TTT_Z, 320000);

  // tree_new->Branch("YY1_X", &YY1_X, 320000);
  // tree_new->Branch("YY1_Y", &YY1_Y, 320000);
  // tree_new->Branch("YY1_Z", &YY1_Z, 320000);

  // Tree entries.
  int tentries = tree->GetEntries();
  cout << "Total Number of Entries: " << tree->GetEntries() << endl;

  // int F3_Counter = 0;
  // int FE9_Counter = 0;
  // int FE12_Counter = 0;
  // int S1_Counter = 0;

  // int F3_50Ca_Counter = 0;
  // int FE9_50Ca_Counter = 0;
  // int FE12_50Ca_Counter = 0;
  // int S1_50Ca_Counter = 0;

  // int F3_50Ca_S0_Counter = 0;
  // int FE9_50Ca_S0_Counter = 0;
  // int FE12_50Ca_S0_Counter = 0;
  // int S1_50Ca_S0_Counter = 0;

  // int S0_Counter[101][101];
  // int S1_Counter[101][101];
  // double Transmission_Counter[101][101];

  // for (int i = 0; i < 101; i++)
  //   {
  //     for (int j = 0; j < 101; j++)
  // 	{
  // 	  S0_Counter[i][j] = 0;
  // 	  S1_Counter[i][j] = 0;
  // 	  Transmission_Counter[i][j] = 0;
  // 	}
  //   }

  // double Ex_Counter[4][51];
  // double Ex_TTT_Counter[4][51];
  // double Ex_YY1_Counter[4][51];
  // for (int i = 0; i < 4; i++)
  // {
  //   for (int j = 0; j < 51; j++)
  //   {
  //     Ex_Counter[i][j] = 0;
  //     Ex_TTT_Counter[i][j] = 0;
  //     Ex_YY1_Counter[i][j] = 0;
  //   }
  // }

  // double Ex_Shift_Counter[31][31];
  // double Ex_TTT_Shift_Counter[31][31];
  // double Ex_YY1_Shift_Counter[31][31];
  // double Ex_Shift_Counter_X[31][31];
  // double Ex_TTT_Shift_Counter_X[31][31];
  // double Ex_YY1_Shift_Counter_X[31][31];
  // double Ex_Shift_Counter_Y[31][31];
  // double Ex_TTT_Shift_Counter_Y[31][31];
  // double Ex_YY1_Shift_Counter_Y[31][31];
  // double Ex_Shift_Counter_XY[31][31];
  // double Ex_TTT_Shift_Counter_XY[31][31];
  // double Ex_YY1_Shift_Counter_XY[31][31];
  // for (int i = 0; i < 31; i++)
  // {
  //   for (int j = 0; j < 31; j++)
  //   {
  //     Ex_Shift_Counter[i][j] = 0;
  //     Ex_TTT_Shift_Counter[i][j] = 0;
  //     Ex_YY1_Shift_Counter[i][j] = 0;
  //     Ex_Shift_Counter_X[i][j] = 0;
  //     Ex_TTT_Shift_Counter_X[i][j] = 0;
  //     Ex_YY1_Shift_Counter_X[i][j] = 0;
  //     Ex_Shift_Counter_Y[i][j] = 0;
  //     Ex_TTT_Shift_Counter_Y[i][j] = 0;
  //     Ex_YY1_Shift_Counter_Y[i][j] = 0;
  //     Ex_Shift_Counter_XY[i][j] = 0;
  //     Ex_TTT_Shift_Counter_XY[i][j] = 0;
  //     Ex_YY1_Shift_Counter_XY[i][j] = 0;
  //   }
  // }

  // for (int i = fEvent; i < lEvent; i++)
  // for(int i = 48378060; i < 48378065; i++)
  for (int i = 0; i < tentries; i++)
  {
    tree->GetEntry(i);

    if (i % 100000 == 0)
    {
      cout << "Event Number = " << i << endl;
    }

    // F3DS Trigger.
    // if (F3_T < -2692.7 || F3_T > -2692.1)
    //  continue;

    // Physics Trigger and F3-S0 Timing correlation.
    // if (F3_T > -2692.7 && F3_T < -2692.1)
    //  continue;

    for (int j = 0; j < 7; j++)
    {
      h_PID_T_Event->Fill(i, PID_T[j]);
      // h_PID_T_Event_C->Fill(i,PID_T[j] + PID_Shift[i/29600]);

      h_S1PID_T_Event->Fill(i, S1PID_T[j]);
      // h_S1PID_T_Event_C->Fill(i,S1PID_T[j] + S1PID_Shift[i/59200]);

      h_Beam_T_Event->Fill(i, Beam_T[j]);
      // h_Beam_T_Event_C->Fill(i,Beam_T[j] + Beam_Shift[i/75000]);

      h_AQ_Event->Fill(i, AQ[j]);
      // h_AQ_Event_C->Fill(i,AQ[j] + AQ_Shift[i/370000]);

      // PID_T[j] = PID_T[j] + PID_Shift[i/29600];
      // Beam_T[j] = Beam_T[j] + Beam_Shift[i/75000];

      // S1PID_T[j] = S1PID_T[j] + S1PID_Shift[i/59200];

      double beta_PID = 68543.07 / (PID_T[j] - 644.49) / 0.000000001 / 299792458000;
      double gamma_PID = sqrt(1 / (1 - pow(beta_PID, 2)));
      double E_Beam_PID = M_50Ca * (gamma_PID - 1);

      double beta_SHARAQ = 14871.68 / (Beam_T[j] + 97.53 - 85 + 2 - 20 - 120 + 150 - 11 + 5.7) / 0.000000001 / 299792458000;
      double gamma_SHARAQ = sqrt(1 / (1 - pow(beta_SHARAQ, 2)));
      double E_Beam_SHARAQ = M_50Ca * (gamma_SHARAQ - 1);

      double E_Beam_temp = E_Beam_SHARAQ;

      E_Beam_temp = (Beam_param[0] * pow(E_Beam_temp / 50., 5) + Beam_param[1] * pow(E_Beam_temp / 50., 4) + Beam_param[2] * pow(E_Beam_temp / 50., 3) + Beam_param[3] * pow(E_Beam_temp / 50., 2) + Beam_param[4] * pow(E_Beam_temp / 50., 1) + Beam_param[5] * pow(E_Beam_temp / 50., 0)) * 50;
      E_Beam_temp = (Beam_param_2[0] * pow(E_Beam_temp / 50., 5) + Beam_param_2[1] * pow(E_Beam_temp / 50., 4) + Beam_param_2[2] * pow(E_Beam_temp / 50., 3) + Beam_param_2[3] * pow(E_Beam_temp / 50., 2) + Beam_param_2[4] * pow(E_Beam_temp / 50., 1) + Beam_param_2[5] * pow(E_Beam_temp / 50., 0)) * 50;

      double gamma_S0 = E_Beam_temp / M_50Ca + 1;
      double beta_S0 = sqrt(1 - 1 / pow(gamma_S0, 2));
      double E_Beam_S0 = M_50Ca * (gamma_S0 - 1);

      E_Beam_temp = (Beam_param_2[0] * pow(E_Beam_temp / 50, 5) + Beam_param_2[1] * pow(E_Beam_temp / 50, 4) + Beam_param_2[2] * pow(E_Beam_temp / 50, 3) + Beam_param_2[3] * pow(E_Beam_temp / 50, 2) + Beam_param_2[4] * pow(E_Beam_temp / 50, 1) + Beam_param_2[5] * pow(E_Beam_temp / 50, 0)) * 50;

      double gamma_OEDO_C = E_Beam_temp / M_50Ca + 1;
      double beta_OEDO_C = sqrt(1 - 1 / pow(gamma_OEDO_C, 2));
      double E_Beam_OEDO_C = M_50Ca * (gamma_OEDO_C - 1);

      double beta_OEDO = 9489.2 / (S1PID_T[j] - 579.06 + 110 - 2 + 30 + 7) / 0.000000001 / 299792458000;
      double gamma_OEDO = sqrt(1 / (1 - pow(beta_OEDO, 2)));
      double E_Beam_OEDO = M_50Ca * (gamma_OEDO - 1);

      // if (E_Beam_SHARAQ > 1.)
      //   {
      //     cout << endl;
      //     cout << "Event " << i << endl;
      //     cout << "SHARAQ Energy = " << E_Beam_SHARAQ/50. << endl;
      //     cout << "SHARAQ beta = " << beta_SHARAQ << endl;
      //     cout << "SHARAQ gamma = " << gamma_SHARAQ << endl;
      //     cout << "S0 Energy = " << E_Beam_S0/50. << endl;
      //     cout << "OEDO_C_Energy = " << E_Beam_OEDO_C/50. << endl;
      //     cout << "OEDO_Energy = " << E_Beam_OEDO/50. << endl;
      //   }

      IC_X[j] = S1_X[j];
      IC_Y[j] = S1_Y[j];

      // S1_X[j] = S1_X[j] + 1250*sin(S1_A[j]);
      // S1_Y[j] = S1_Y[j] + 1250*sin(S1_B[j]);

      // For My Analysis.
      IC_Brho_array[j] = -1000000;
      beta_staraq_array[j] = -1000000;
      beta_sharaq_f_oedo_array[j] = -1000000;

      AQ[j] = -1000000;

      AQ_2[j] = -1000000;

      // TTT_Theta[j] = -1000000;
      // TTT_Theta_X[j] = -1000000;
      // TTT_Theta_Y[j] = -1000000;
      // TTT_Theta_XY[j] = -1000000;
      // TTT_Phi[j] = -1000000;

      // YY1_Theta[j] = -1000000;
      // YY1_Theta_X[j] = -1000000;
      // YY1_Theta_Y[j] = -1000000;
      // YY1_Theta_XY[j] = -1000000;
      // YY1_Phi[j] = -1000000;

      // S0_X[j] = -S0_X[j];
      // S0_Y[j] = -S0_Y[j];
      // S0_A[j] = -S0_A[j];
      // S0_B[j] = -S0_B[j];

      // Beam Energy at S0 calculation.

      // E_Beam = E_Beam_SHARAQ;

      // double temp = Beam_param[0]*pow(E_Beam,8) + Beam_param[1]*pow(E_Beam,7) + Beam_param[2]*pow(E_Beam,6) + Beam_param[3]*pow(E_Beam,5) + Beam_param[4]*pow(E_Beam,4) + Beam_param[5]*pow(E_Beam,3) + Beam_param[6]*pow(E_Beam,2) + Beam_param[7]*pow(E_Beam,1) + Beam_param[8]*pow(E_Beam,0); //We make an ATIMA correction of the beam energy loss at FE12 SR-PPACs.
      // E_Loss = (E_Beam - temp)/sin(PI/2. + sqrt(pow(S0_A[j]/1000.,2) + pow(S0_B[j]/1000.,2)));
      // E_Beam_C[j] = E_Beam - E_Loss;

      E_Beam_C[j] = E_Beam_S0 / 50.;

      // temp = Beam_param_2[0]*pow(E_Beam_C[j],8) + Beam_param_2[1]*pow(E_Beam_C[j],7) + Beam_param_2[2]*pow(E_Beam_C[j],6) + Beam_param_2[3]*pow(E_Beam_C[j],5) + Beam_param_2[4]*pow(E_Beam_C[j],4) + Beam_param_2[5]*pow(E_Beam_C[j],3) + Beam_param_2[6]*pow(E_Beam_C[j],2) + Beam_param_2[7]*pow(E_Beam_C[j],1) + Beam_param_2[8]*pow(E_Beam_C[j],0); //We make an ATIMA correction of the beam energy at the CD2 target.
      // E_Loss = (E_Beam_C[j] - temp)/sin(PI/2. + sqrt(pow(S0_A[j]/1000.,2) + pow(S0_B[j]/1000.,2)));
      // E_Beam_C[j] = (E_Beam_C[j] - E_Loss)/50.;

      E_Beam_2[j] = E_Beam_OEDO / 50.;
      // E_Beam_C[j] = (E_Beam_C[j] - E_Loss)/49.;

      // if (PID_T[j] > 1)
      //   {
      //     cout << "beta = " << beta << endl;
      //     cout << "gamma = " << gamma << endl;
      //     cout << "E_Beam = " << E_Beam << endl;
      //     cout << "E_Beam_C = " << E_Beam_C[j] << endl;
      //     cout << endl;
      //   }

      // if (E_Beam > 1)
      //   {
      //     cout << endl;
      //     cout << E_Beam/50. << endl;
      //     cout << E_Beam_C[j] << endl;
      //   }

      // E_Beam_C[j] = E_Beam_C[j]/51;
      // E_Beam_C[j] = E_Beam_C[j]/49;

      // PID at S1.
      IC_Brho = IC_Brho_zero * (1 + (S1_X[j] - IC_X0 * S0_X[j]) / (IC_delta + IC_delta_a * S0_A[j])); // Equation provided by Thomas Chillery.

      // For My Analysis.
      IC_Brho_array[j] = IC_Brho;
      beta_staraq_array[j] = beta_OEDO;          // beta SHARAQ from FE12-S1 TOF
      beta_sharaq_f_oedo_array[j] = beta_OEDO_C; // beta SHARAQ from OEDO beam Energy

      // beta = sqrt(2*(E_Beam_S1)*50/M_50Ca); //Assumption of no energy loss.

      // double beta_1 = 14871.68/(Beam_T_1[j] + 97.53 - 85 + 2 - 20 - 120 + 150 - 11)/0.000000001/299792458000;
      // double gamma_1 = sqrt(1/(1-pow(beta_1,2)));
      // AQ_Uno[j] = IC_Brho*0.32184/(beta_1*gamma_1);

      AQ[j] = IC_Brho * 0.32184 / (beta_OEDO_C * gamma_OEDO_C);

      AQ[j] = AQ[j] - 0.001567 * S0_X[j];
      AQ[j] = AQ[j] + 0.000003749 * pow(S1_A[j] * 1000, 2) - 0.00006171 * S1_A[j] * 1000;
      AQ[j] = AQ[j] - 0.0003949 * S0_A[j] * 1000;
      AQ[j] = AQ[j] - 0.001567 * S0_X[j];
      AQ[j] = AQ[j] + 0.001719 * S0_X[j];
      AQ[j] = AQ[j] - 0.00001370 * pow(S0_A[j] * 1000, 2) - 0.0001998 * S0_A[j] * 1000;

      // AQ[j] = AQ[j] + AQ_Shift[i/370000];

      AQ[j] = AQ[j] * 1.1418 - 0.3587;

      AQ_2[j] = IC_Brho * 0.32184 / (beta_OEDO * gamma_OEDO);

      AQ_2[j] = AQ_2[j] - 0.000224 * S1_Y[j];
      AQ_2[j] = AQ_2[j] - 0.000955 * S0_X[j];
      AQ_2[j] = AQ_2[j] - 0.000180 * S1_X[j];
      AQ_2[j] = AQ_2[j] + 0.000779 * S0_A[j] * 1000;
      AQ_2[j] = AQ_2[j] - 0.000229 * S0_A[j] * 1000;
      AQ_2[j] = AQ_2[j] + 0.00000428 * pow(S1_A[j] * 1000, 2) + 0.0000262 * S1_A[j] * 1000;
      AQ_2[j] = AQ_2[j] + 0.000198 * S0_Y[j];
      AQ_2[j] = AQ_2[j] + 0.000114 * S1_B[j] * 1000;
      AQ_2[j] = AQ_2[j] - 0.000022609 * pow(S0_A[j] * 1000, 2) - 0.0001201 * S0_A[j] * 1000;
      AQ_2[j] = AQ_2[j] - 0.0002295 * S0_A[j] * 1000;
      AQ_2[j] = AQ_2[j] + 0.000009783 * pow(S0_A[j] * 1000, 2) + 0.00006245 * S0_A[j] * 1000;

      if (S1_Y[j] > -78.5 && S1_Y[j] < 78.5)
      {
        AQ_2[j] = AQ_2[j] + AQ_2_S1_Y[int(2 * S1_Y[j] + 160)];
      }

      if (S1_X[j] > -139.5 && S1_X[j] < 139.5)
      {
        AQ_2[j] = AQ_2[j] + AQ_2_S1_X[int(2 * S1_X[j] + 280)];
      }

      AQ_2[j] = AQ_2[j] * 1.1927 - 0.4819;
    }

    // // CsI Energy Calibration.
    // for (int j = 0; j < 16; j++)
    // {
    //   CsI_E[j] = (CsI_E[j] + CsI_N[j]) * CsI_M[j];
    // }

    // // TTT Energy Analysis.
    // TTT_Front_E = TTT_F_E[0];
    // TTT_Back_E = TTT_B_E[0];
    // TTT_Front_ID = TTT_F_ID[0];
    // TTT_Back_ID = TTT_B_ID[0];

    // if (std::abs(TTT_B_ID[0] - TTT_B_ID[1]) == 1 && TTT_B_E[1] > 1.)
    // {
    //   for (int j = 0; j < 512; j++)
    //   {
    //     if (TTT_B_ID[0] == j)
    //     {
    //       TTT_Back_E = 1.075 * ((TTT_B_E[0] + TTT_B_E[1]) * list[0][j]) - 0.247;
    //     }
    //   }
    // }

    // // TTT-CsI correlation.
    // for (int k = 0; k < 16; k++)
    // {
    //   if (TTT_Back_ID > Back_Min[k] && TTT_Back_ID < (Back_Min[k] + 65) && TTT_Front_ID > Front_Min[k] && TTT_Front_ID < (Front_Min[k] + 65))
    //   {
    //     CsI_Energy = CsI_E[ID[k]];
    //     h_TTT_CsI[ID[k]]->Fill(CsI_Energy, TTT_Front_E);
    //     h_TTT_CsI_all->Fill(CsI_Energy, TTT_Front_E);
    //     CsI_Time = CsI_T_RAW[Time_ID[k]];
    //   }
    // }

    // // YY1 Energy analysis
    // YY1_Energy = YY1_E[0];
    // YY1_Energy_ID = YY1_ID[0];

    // // YY1-Time correlation.
    // if (YY1_Energy_ID > -1 && YY1_Energy_ID < 16)
    // {
    //   YY1_Time = YY1_T_RAW[0];
    // }

    // if (YY1_Energy_ID > 15 && YY1_Energy_ID < 32)
    // {
    //   YY1_Time = YY1_T_RAW[1];
    // }

    // if (YY1_Energy_ID > 31 && YY1_Energy_ID < 48)
    // {
    //   YY1_Time = YY1_T_RAW[2];
    // }

    // if (YY1_Energy_ID > 47 && YY1_Energy_ID < 64)
    // {
    //   YY1_Time = YY1_T_RAW[3];
    // }

    // if (YY1_Energy_ID > 63 && YY1_Energy_ID < 80)
    // {
    //   YY1_Time = YY1_T_RAW[5];
    // }

    // if (YY1_Energy_ID > 79 && YY1_Energy_ID < 96)
    // {
    //   YY1_Time = YY1_T_RAW[4];
    // }

    // Scattering angle and incident angle.

    // TTT_X = -1000000;
    // TTT_Y = -1000000;
    // TTT_Z = -1000000;

    // if (TTT_Front_ID > -1 && TTT_Back_ID > -1 && TTT_Front_E > 0)
    // {
    //   TTT_X = f_Coordinates(TTT_Front_ID, TTT_Back_ID)[0];
    //   TTT_Y = f_Coordinates(TTT_Front_ID, TTT_Back_ID)[1];
    //   TTT_Z = f_Coordinates(TTT_Front_ID, TTT_Back_ID)[2];

    //   for (int j = 0; j < 7; j++)
    //   {
    //     TTT_Theta[j] = -1000000;
    //     TTT_Theta_X[j] = -1000000;
    //     TTT_Theta_Y[j] = -1000000;
    //     TTT_Theta_XY[j] = -1000000;
    //     TTT_Phi[j] = -1000000;
    //     TTT_Excitation[j] = -1000000;
    //     TTT_Excitation_X[j] = -1000000;
    //     TTT_Excitation_Y[j] = -1000000;
    //     TTT_Excitation_XY[j] = -1000000;
    //     TTT_Excitation_Mean[j] = -1000000;

    //     TINA2_E = -1000000;

    //     if (S0_X[j] > -100000 && S0_Y[j] > -100000 && S0_A[j] > -100000 && S0_B[j] > -100000)
    //     {
    //       TTT_Theta[j] = f_Theta(TTT_X - S0_X[j], TTT_Y - S0_Y[j], TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1);
    //       TTT_Phi[j] = f_Phi(TTT_X, TTT_Y);

    //       // cout << "hola" << endl;
    //       if ((CsI_Banana_Bottom[0] * pow(CsI_Energy * sin(TTT_Theta[j]), 6) + CsI_Banana_Bottom[1] * pow(CsI_Energy * sin(TTT_Theta[j]), 5) + CsI_Banana_Bottom[2] * pow(CsI_Energy * sin(TTT_Theta[j]), 4) + CsI_Banana_Bottom[3] * pow(CsI_Energy * sin(TTT_Theta[j]), 3) + CsI_Banana_Bottom[4] * pow(CsI_Energy * sin(TTT_Theta[j]), 2) + CsI_Banana_Bottom[5] * pow(CsI_Energy * sin(TTT_Theta[j]), 1) + CsI_Banana_Bottom[6] * pow(CsI_Energy * sin(TTT_Theta[j]), 0)) < TTT_Front_E && (CsI_Banana_Top[0] * pow(CsI_Energy * sin(TTT_Theta[j]), 6) + CsI_Banana_Top[1] * pow(CsI_Energy * sin(TTT_Theta[j]), 5) + CsI_Banana_Top[2] * pow(CsI_Energy * sin(TTT_Theta[j]), 4) + CsI_Banana_Top[3] * pow(CsI_Energy * sin(TTT_Theta[j]), 3) + CsI_Banana_Top[4] * pow(CsI_Energy * sin(TTT_Theta[j]), 2) + CsI_Banana_Top[5] * pow(CsI_Energy * sin(TTT_Theta[j]), 1) + CsI_Banana_Top[6] * pow(CsI_Energy * sin(TTT_Theta[j]), 0)) > TTT_Front_E && CsI_Energy > 0.2)
    //       {
    //         double TTT_Saver = TTT_Front_E;

    //         TTT_Front_E = f_Energy(TTT_Front_E, TTT_Theta[j]);

    //         h_CsI_Gate->Fill(CsI_Energy, TTT_Front_E);

    //         // cout << "hola" << endl;
    //         TTT_Theta_X[j] = f_Theta(TTT_X + S0_X[j], TTT_Y - S0_Y[j], TTT_Z, tan(-S0_A[j]), tan(S0_B[j]), 1);
    //         TTT_Excitation_X[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //         TTT_Theta_Y[j] = f_Theta(TTT_X - S0_X[j], TTT_Y + S0_Y[j], TTT_Z, tan(S0_A[j]), tan(-S0_B[j]), 1);
    //         TTT_Excitation_Y[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //         TTT_Theta_XY[j] = f_Theta(TTT_X + S0_X[j], TTT_Y + S0_Y[j], TTT_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1);
    //         TTT_Excitation_XY[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //         TTT_Theta[j] = f_Theta(TTT_X - S0_X[j], TTT_Y - S0_Y[j], TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1);
    //         TTT_Excitation[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //         TTT_Excitation_Mean[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], 15.2, sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));
    //         // cout << endl;
    //         // cout << "Event Number = " << i << endl;
    //         // cout << "TTT_E = " << TTT_Front_E << endl;
    //         // cout << "CsI_E = " << CsI_Energy << endl;

    //         if ((FE9_X[j] < (FS_50Ca_M_AB * PID_T[j] + FS_50Ca_N_AB)) && (FE9_X[j] > (FS_50Ca_M_BC * PID_T[j] + FS_50Ca_N_BC)) && (FE9_X[j] > (FS_50Ca_M_CD * PID_T[j] + FS_50Ca_N_CD)) && (FE9_X[j] < (FS_50Ca_M_DA * PID_T[j] + FS_50Ca_N_DA)))
    //         {
    //           if (F3_T[j] < (-8.913 * pow(TTT_Front_E, 2) + 102.08 * pow(TTT_Front_E, 1) - 3478) && F3_T[j] > (-4.259 * pow(TTT_Front_E, 2) + 66.11 * pow(TTT_Front_E, 1) - 3543.5))
    //           {
    //             if (TTT_Front_E > 1.7 && TTT_Front_E < 6.5 && abs(TTT_Front_E - TTT_Back_E) < 0.25)
    //             {
    //               if ((pow(AQ[j] - 2.515, 2) / pow(0.05, 2) + pow(AQ_2[j] - 2.547, 2) / pow(0.03, 2) < 1.) || (pow(AQ[j] - 2.65, 2) / pow(0.03, 2) + pow(AQ_2[j] - 2.675, 2) / pow(0.02, 2) < 1.))
    //               {
    //                 for (int k = 0; k < 51; k++)
    //                 {
    //                   h_Ex_TTT[k][0]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TTT[k][1]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TTT[k][2]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TTT[k][3]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                   h_Ex_TINA[k][0]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TINA[k][1]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TINA[k][2]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TINA[k][3]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                   if (f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                   {
    //                     Ex_Counter[0][k] = Ex_Counter[0][k] + 1;
    //                     Ex_TTT_Counter[0][k] = Ex_TTT_Counter[0][k] + 1;
    //                   }

    //                   if (f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                   {
    //                     Ex_Counter[1][k] = Ex_Counter[1][k] + 1;
    //                     Ex_TTT_Counter[1][k] = Ex_TTT_Counter[1][k] + 1;
    //                   }

    //                   if (f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                   {
    //                     Ex_Counter[2][k] = Ex_Counter[2][k] + 1;
    //                     Ex_TTT_Counter[2][k] = Ex_TTT_Counter[2][k] + 1;
    //                   }

    //                   if (f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                   {
    //                     Ex_Counter[3][k] = Ex_Counter[3][k] + 1;
    //                     Ex_TTT_Counter[3][k] = Ex_TTT_Counter[3][k] + 1;
    //                   }

    //                   if (k == 35)
    //                   {
    //                     for (int l = 0; l < 31; l++)
    //                     {
    //                       for (int m = 0; m < 31; m++)
    //                       {
    //                         h_Ex_TINA_Shift[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TTT_Shift[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TINA_Shift_X[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TTT_Shift_X[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TINA_Shift_Y[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TTT_Shift_Y[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TINA_Shift_XY[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TTT_Shift_XY[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                         if (f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                         {
    //                           Ex_Shift_Counter[l][m] = Ex_Shift_Counter[l][m] + 1;
    //                           Ex_TTT_Shift_Counter[l][m] = Ex_TTT_Shift_Counter[l][m] + 1;
    //                         }

    //                         if (f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(-S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                         {
    //                           Ex_Shift_Counter_X[l][m] = Ex_Shift_Counter_X[l][m] + 1;
    //                           Ex_TTT_Shift_Counter_X[l][m] = Ex_TTT_Shift_Counter_X[l][m] + 1;
    //                         }

    //                         if (f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                         {
    //                           Ex_Shift_Counter_Y[l][m] = Ex_Shift_Counter_Y[l][m] + 1;
    //                           Ex_TTT_Shift_Counter_Y[l][m] = Ex_TTT_Shift_Counter_Y[l][m] + 1;
    //                         }

    //                         if (f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                         {
    //                           Ex_Shift_Counter_XY[l][m] = Ex_Shift_Counter_XY[l][m] + 1;
    //                           Ex_TTT_Shift_Counter_XY[l][m] = Ex_TTT_Shift_Counter_XY[l][m] + 1;
    //                         }
    //                       }
    //                     }
    //                   }
    //                 }
    //               }
    //             }
    //           }
    //         }
    //         TTT_Front_E = TTT_Saver;
    //       }

    //       if (CsI_Energy < 0.2 && TTT_Front_E < 6.5)
    //       {
    //         TTT_Theta_X[j] = f_Theta(TTT_X + S0_X[j], TTT_Y - S0_Y[j], TTT_Z, tan(-S0_A[j]), tan(S0_B[j]), 1);
    //         TTT_Excitation_X[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //         TTT_Theta_Y[j] = f_Theta(TTT_X - S0_X[j], TTT_Y + S0_Y[j], TTT_Z, tan(S0_A[j]), tan(-S0_B[j]), 1);
    //         TTT_Excitation_Y[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //         TTT_Theta_XY[j] = f_Theta(TTT_X + S0_X[j], TTT_Y + S0_Y[j], TTT_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1);
    //         TTT_Excitation_XY[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //         TTT_Theta[j] = f_Theta(TTT_X - S0_X[j], TTT_Y - S0_Y[j], TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1);
    //         TTT_Excitation[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));
    //         TTT_Excitation_Mean[j] = f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], 15.2, sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));
    //         // cout << endl;
    //         // cout << "Event Number = " << i << endl;
    //         // cout << "TTT_E = " << TTT_Front_E << endl;
    //         // cout << "CsI_E = " << CsI_Energy << endl;

    //         if ((FE9_X[j] < (FS_50Ca_M_AB * PID_T[j] + FS_50Ca_N_AB)) && (FE9_X[j] > (FS_50Ca_M_BC * PID_T[j] + FS_50Ca_N_BC)) && (FE9_X[j] > (FS_50Ca_M_CD * PID_T[j] + FS_50Ca_N_CD)) && (FE9_X[j] < (FS_50Ca_M_DA * PID_T[j] + FS_50Ca_N_DA)))
    //         {
    //           if (F3_T[j] < (-8.913 * pow(TTT_Front_E, 2) + 102.08 * pow(TTT_Front_E, 1) - 3478) && F3_T[j] > (-4.259 * pow(TTT_Front_E, 2) + 66.11 * pow(TTT_Front_E, 1) - 3543.5))
    //           {
    //             if (TTT_Front_E > 1.7 && TTT_Front_E < 6.5 && abs(TTT_Front_E - TTT_Back_E) < 0.25)
    //             {
    //               if ((pow(AQ[j] - 2.515, 2) / pow(0.05, 2) + pow(AQ_2[j] - 2.547, 2) / pow(0.03, 2) < 1.) || (pow(AQ[j] - 2.65, 2) / pow(0.03, 2) + pow(AQ_2[j] - 2.675, 2) / pow(0.02, 2) < 1.))
    //               {
    //                 for (int k = 0; k < 51; k++)
    //                 {
    //                   h_Ex_TTT[k][0]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TTT[k][1]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TTT[k][2]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TTT[k][3]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                   h_Ex_TINA[k][0]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TINA[k][1]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TINA[k][2]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                   h_Ex_TINA[k][3]->Fill(f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                   if (f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, TTT_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                   {
    //                     Ex_Counter[0][k] = Ex_Counter[0][k] + 1;
    //                     Ex_TTT_Counter[0][k] = Ex_TTT_Counter[0][k] + 1;
    //                   }

    //                   if (f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, TTT_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                   {
    //                     Ex_Counter[1][k] = Ex_Counter[1][k] + 1;
    //                     Ex_TTT_Counter[1][k] = Ex_TTT_Counter[1][k] + 1;
    //                   }

    //                   if (f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, TTT_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                   {
    //                     Ex_Counter[2][k] = Ex_Counter[2][k] + 1;
    //                     Ex_TTT_Counter[2][k] = Ex_TTT_Counter[2][k] + 1;
    //                   }

    //                   if (f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, TTT_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                   {
    //                     Ex_Counter[3][k] = Ex_Counter[3][k] + 1;
    //                     Ex_TTT_Counter[3][k] = Ex_TTT_Counter[3][k] + 1;
    //                   }

    //                   if (k == 35)
    //                   {
    //                     for (int l = 0; l < 31; l++)
    //                     {
    //                       for (int m = 0; m < 31; m++)
    //                       {
    //                         h_Ex_TINA_Shift[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TTT_Shift[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TINA_Shift_X[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TTT_Shift_X[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TINA_Shift_Y[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TTT_Shift_Y[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TINA_Shift_XY[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                         h_Ex_TTT_Shift_XY[l][m]->Fill(f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                         if (f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                         {
    //                           Ex_Shift_Counter[l][m] = Ex_Shift_Counter[l][m] + 1;
    //                           Ex_TTT_Shift_Counter[l][m] = Ex_TTT_Shift_Counter[l][m] + 1;
    //                         }

    //                         if (f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y - S0_Y[j] + (m - 15), TTT_Z, tan(-S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                         {
    //                           Ex_Shift_Counter_X[l][m] = Ex_Shift_Counter_X[l][m] + 1;
    //                           Ex_TTT_Shift_Counter_X[l][m] = Ex_TTT_Shift_Counter_X[l][m] + 1;
    //                         }

    //                         if (f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X - S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                         {
    //                           Ex_Shift_Counter_Y[l][m] = Ex_Shift_Counter_Y[l][m] + 1;
    //                           Ex_TTT_Shift_Counter_Y[l][m] = Ex_TTT_Shift_Counter_Y[l][m] + 1;
    //                         }

    //                         if (f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(TTT_Front_E, f_Theta(TTT_X + S0_X[j] + (l - 15), TTT_Y + S0_Y[j] + (m - 15), TTT_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                         {
    //                           Ex_Shift_Counter_XY[l][m] = Ex_Shift_Counter_XY[l][m] + 1;
    //                           Ex_TTT_Shift_Counter_XY[l][m] = Ex_TTT_Shift_Counter_XY[l][m] + 1;
    //                         }
    //                       }
    //                     }
    //                   }
    //                 }
    //               }
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
    // }

    // YY1_X = -1000000;
    // YY1_Y = -1000000;
    // YY1_Z = -1000000;

    // if (YY1_Energy_ID > -1)
    // {
    //   YY1_X = YY1_r_list[YY1_ID_list[YY1_Energy_ID % 16]] * sin(YY1_theta_list[YY1_ID_list[YY1_Energy_ID % 16]]) * cos(YY1_phi_list[YY1_Energy_ID / 16]);
    //   YY1_Y = YY1_r_list[YY1_ID_list[YY1_Energy_ID % 16]] * sin(YY1_theta_list[YY1_ID_list[YY1_Energy_ID % 16]]) * sin(YY1_phi_list[YY1_Energy_ID / 16]);
    //   YY1_Z = YY1_r_list[YY1_ID_list[YY1_Energy_ID % 16]] * cos(YY1_theta_list[YY1_ID_list[YY1_Energy_ID % 16]]);

    //   // YY1_X = YY1_r_list[YY1_Energy_ID%16]*sin(YY1_theta_list[YY1_Energy_ID%16])*cos(YY1_phi_list[YY1_Energy_ID/16]);
    //   // YY1_Y = YY1_r_list[YY1_Energy_ID%16]*sin(YY1_theta_list[YY1_Energy_ID%16])*sin(YY1_phi_list[YY1_Energy_ID/16]);
    //   // YY1_Z = YY1_r_list[YY1_Energy_ID%16]*cos(YY1_theta_list[YY1_Energy_ID%16]);

    //   for (int j = 0; j < 7; j++)
    //   {
    //     YY1_Theta[j] = -1000000;
    //     YY1_Theta_X[j] = -1000000;
    //     YY1_Theta_Y[j] = -1000000;
    //     YY1_Theta_XY[j] = -1000000;
    //     YY1_Excitation[j] = -1000000;
    //     YY1_Excitation_X[j] = -1000000;
    //     YY1_Excitation_Y[j] = -1000000;
    //     YY1_Excitation_XY[j] = -1000000;
    //     YY1_Excitation_Mean[j] = -1000000;

    //     if (S0_X[j] > -100000 && S0_Y[j] > -100000 && S0_A[j] > -100000 && S0_B[j] > -100000)
    //     {
    //       YY1_Phi[j] = f_Phi(YY1_X, YY1_Y);

    //       YY1_Theta_X[j] = f_Theta(YY1_X + S0_X[j], YY1_Y - S0_Y[j], YY1_Z, tan(-S0_A[j]), tan(S0_B[j]), 1);
    //       YY1_Excitation_X[j] = f_Excitation_Energy(YY1_Energy, YY1_Theta_X[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //       YY1_Theta_Y[j] = f_Theta(YY1_X - S0_X[j], YY1_Y + S0_Y[j], YY1_Z, tan(S0_A[j]), tan(-S0_B[j]), 1);
    //       YY1_Excitation_Y[j] = f_Excitation_Energy(YY1_Energy, YY1_Theta_Y[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //       YY1_Theta_XY[j] = f_Theta(YY1_X + S0_X[j], YY1_Y + S0_Y[j], YY1_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1);
    //       YY1_Excitation_XY[j] = f_Excitation_Energy(YY1_Energy, YY1_Theta_XY[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //       YY1_Theta[j] = f_Theta(YY1_X - S0_X[j], YY1_Y - S0_Y[j], YY1_Z, tan(S0_A[j]), tan(S0_B[j]), 1);
    //       YY1_Excitation[j] = f_Excitation_Energy(YY1_Energy, YY1_Theta[j], E_Beam_C[j], sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));
    //       YY1_Excitation_Mean[j] = f_Excitation_Energy(YY1_Energy, YY1_Theta[j], 15.2, sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2)));

    //       if ((FE9_X[j] < (FS_50Ca_M_AB * PID_T[j] + FS_50Ca_N_AB)) && (FE9_X[j] > (FS_50Ca_M_BC * PID_T[j] + FS_50Ca_N_BC)) && (FE9_X[j] > (FS_50Ca_M_CD * PID_T[j] + FS_50Ca_N_CD)) && (FE9_X[j] < (FS_50Ca_M_DA * PID_T[j] + FS_50Ca_N_DA)))
    //       {
    //         if (F3_T[j] > -2500 && F3_T[j] < -2420)
    //         {
    //           if (YY1_Energy > 1.7 && YY1_Energy < 4.5)
    //           {
    //             if ((pow(AQ[j] - 2.515, 2) / pow(0.05, 2) + pow(AQ_2[j] - 2.547, 2) / pow(0.03, 2) < 1.) || (pow(AQ[j] - 2.65, 2) / pow(0.03, 2) + pow(AQ_2[j] - 2.675, 2) / pow(0.02, 2) < 1.))
    //             {
    //               for (int k = 0; k < 51; k++)
    //               {
    //                 h_Ex_YY1[k][0]->Fill(f_Excitation_Energy(YY1_Energy, YY1_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                 h_Ex_YY1[k][1]->Fill(f_Excitation_Energy(YY1_Energy, YY1_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                 h_Ex_YY1[k][2]->Fill(f_Excitation_Energy(YY1_Energy, YY1_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                 h_Ex_YY1[k][3]->Fill(f_Excitation_Energy(YY1_Energy, YY1_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                 h_Ex_TINA[k][0]->Fill(f_Excitation_Energy(YY1_Energy, YY1_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                 h_Ex_TINA[k][1]->Fill(f_Excitation_Energy(YY1_Energy, YY1_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                 h_Ex_TINA[k][2]->Fill(f_Excitation_Energy(YY1_Energy, YY1_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                 h_Ex_TINA[k][3]->Fill(f_Excitation_Energy(YY1_Energy, YY1_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                 if (f_Excitation_Energy(YY1_Energy, YY1_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(YY1_Energy, YY1_Theta[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                 {
    //                   Ex_Counter[0][k] = Ex_Counter[0][k] + 1;
    //                   Ex_YY1_Counter[0][k] = Ex_YY1_Counter[0][k] + 1;
    //                 }
    //                 if (f_Excitation_Energy(YY1_Energy, YY1_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(YY1_Energy, YY1_Theta_X[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                 {
    //                   Ex_Counter[1][k] = Ex_Counter[1][k] + 1;
    //                   Ex_YY1_Counter[1][k] = Ex_YY1_Counter[1][k] + 1;
    //                 }
    //                 if (f_Excitation_Energy(YY1_Energy, YY1_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(YY1_Energy, YY1_Theta_Y[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                 {
    //                   Ex_Counter[2][k] = Ex_Counter[2][k] + 1;
    //                   Ex_YY1_Counter[2][k] = Ex_YY1_Counter[2][k] + 1;
    //                 }
    //                 if (f_Excitation_Energy(YY1_Energy, YY1_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(YY1_Energy, YY1_Theta_XY[j], E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                 {
    //                   Ex_Counter[3][k] = Ex_Counter[3][k] + 1;
    //                   Ex_YY1_Counter[3][k] = Ex_YY1_Counter[3][k] + 1;
    //                 }

    //                 if (k == 35)
    //                 {
    //                   for (int l = 0; l < 31; l++)
    //                   {
    //                     for (int m = 0; m < 31; m++)
    //                     {
    //                       h_Ex_TINA_Shift[l][m]->Fill(f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X - S0_X[j] + (l - 15), YY1_Y - S0_Y[j] + (m - 15), YY1_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                       h_Ex_YY1_Shift[l][m]->Fill(f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X - S0_X[j] + (l - 15), YY1_Y - S0_Y[j] + (m - 15), YY1_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                       h_Ex_TINA_Shift_X[l][m]->Fill(f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X + S0_X[j] + (l - 15), YY1_Y - S0_Y[j] + (m - 15), YY1_Z, tan(-S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                       h_Ex_YY1_Shift_X[l][m]->Fill(f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X + S0_X[j] + (l - 15), YY1_Y - S0_Y[j] + (m - 15), YY1_Z, tan(-S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                       h_Ex_TINA_Shift_Y[l][m]->Fill(f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X - S0_X[j] + (l - 15), YY1_Y + S0_Y[j] + (m - 15), YY1_Z, tan(S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                       h_Ex_YY1_Shift_Y[l][m]->Fill(f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X - S0_X[j] + (l - 15), YY1_Y + S0_Y[j] + (m - 15), YY1_Z, tan(S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                       h_Ex_TINA_Shift_XY[l][m]->Fill(f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X + S0_X[j] + (l - 15), YY1_Y + S0_Y[j] + (m - 15), YY1_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));
    //                       h_Ex_YY1_Shift_XY[l][m]->Fill(f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X + S0_X[j] + (l - 15), YY1_Y + S0_Y[j] + (m - 15), YY1_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))));

    //                       if (f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X - S0_X[j] + (l - 15), YY1_Y - S0_Y[j] + (m - 15), YY1_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X - S0_X[j] + (l - 15), YY1_Y - S0_Y[j] + (m - 15), YY1_Z, tan(S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                       {
    //                         Ex_Shift_Counter[l][m] = Ex_Shift_Counter[l][m] + 1;
    //                         Ex_YY1_Shift_Counter[l][m] = Ex_YY1_Shift_Counter[l][m] + 1;
    //                       }

    //                       if (f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X + S0_X[j] + (l - 15), YY1_Y - S0_Y[j] + (m - 15), YY1_Z, tan(-S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X + S0_X[j] + (l - 15), YY1_Y - S0_Y[j] + (m - 15), YY1_Z, tan(-S0_A[j]), tan(S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                       {
    //                         Ex_Shift_Counter_X[l][m] = Ex_Shift_Counter_X[l][m] + 1;
    //                         Ex_YY1_Shift_Counter_X[l][m] = Ex_YY1_Shift_Counter_X[l][m] + 1;
    //                       }

    //                       if (f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X - S0_X[j] + (l - 15), YY1_Y + S0_Y[j] + (m - 15), YY1_Z, tan(S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X - S0_X[j] + (l - 15), YY1_Y + S0_Y[j] + (m - 15), YY1_Z, tan(S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                       {
    //                         Ex_Shift_Counter_Y[l][m] = Ex_Shift_Counter_Y[l][m] + 1;
    //                         Ex_YY1_Shift_Counter_Y[l][m] = Ex_YY1_Shift_Counter_Y[l][m] + 1;
    //                       }

    //                       if (f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X + S0_X[j] + (l - 15), YY1_Y + S0_Y[j] + (m - 15), YY1_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) > -0.5 && f_Excitation_Energy(YY1_Energy, f_Theta(YY1_X + S0_X[j] + (l - 15), YY1_Y + S0_Y[j] + (m - 15), YY1_Z, tan(-S0_A[j]), tan(-S0_B[j]), 1), E_Beam_C[j] + (0.1 * k - 3.5), sqrt(pow(S0_A[j], 2) + pow(S0_B[j], 2))) < 5.3)
    //                       {
    //                         Ex_Shift_Counter_XY[l][m] = Ex_Shift_Counter_XY[l][m] + 1;
    //                         Ex_YY1_Shift_Counter_XY[l][m] = Ex_YY1_Shift_Counter_XY[l][m] + 1;
    //                       }
    //                     }
    //                   }
    //                 }
    //               }
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
    // }
    tree_new->Fill();
  }

  // cout << F3_Counter << endl;
  // cout << FE9_Counter << endl;
  // cout << FE12_Counter << endl;
  // cout << S1_Counter << endl;

  // cout << F3_50Ca_Counter << endl;
  // cout << FE9_50Ca_Counter << endl;
  // cout << FE12_50Ca_Counter << endl;
  // cout << S1_50Ca_Counter << endl;

  // cout << F3_50Ca_S0_Counter << endl;
  // cout << FE9_50Ca_S0_Counter << endl;
  // cout << FE12_50Ca_S0_Counter << endl;
  // cout << S1_50Ca_S0_Counter << endl;

  // for (int k = 0; k < 101; k++)
  //   {
  //     for (int l = 0; l < 101; l++)
  // 	{
  // 	  if (S0_Counter[k][l] > 0)
  // 	    {
  // 	      Transmission_Counter[k][l] = (S1_Counter[k][l]*1.)/(S0_Counter[k][l]*1.)/20.;
  // 	      //cout << Transmission_X_Counter[k][l] << endl;
  // 	      //cout << Transmission_Y_Counter[k][l] << endl;
  // 	      h_S0->SetBinContent(1 + k,1 + l,S0_Counter[k][l]);
  // 	      h_S0_S1Gate->SetBinContent(1 + k,1 + l,S1_Counter[k][l]);
  // 	      h_Transmission->SetBinContent(1 + k,1 + l,Transmission_Counter[k][l]);
  // 	    }
  // 	}
  //   }

  // for (int i = 0; i < 4; i++)
  // {
  //   for (int j = 0; j < 51; j++)
  //   {
  //     h_Ex->SetBinContent(i + 1, j + 1, Ex_Counter[i][j]);
  //     h_TTT_Ex->SetBinContent(i + 1, j + 1, Ex_TTT_Counter[i][j]);
  //     h_YY1_Ex->SetBinContent(i + 1, j + 1, Ex_YY1_Counter[i][j]);
  //   }
  // }

  // for (int i = 0; i < 31; i++)
  // {
  //   for (int j = 0; j < 31; j++)
  //   {
  //     h_Ex_Shift->SetBinContent(i + 1, j + 1, Ex_Shift_Counter[i][j]);
  //     h_TTT_Ex_Shift->SetBinContent(i + 1, j + 1, Ex_TTT_Shift_Counter[i][j]);
  //     h_YY1_Ex_Shift->SetBinContent(i + 1, j + 1, Ex_YY1_Shift_Counter[i][j]);
  //     h_Ex_Shift_X->SetBinContent(i + 1, j + 1, Ex_Shift_Counter_X[i][j]);
  //     h_TTT_Ex_Shift_X->SetBinContent(i + 1, j + 1, Ex_TTT_Shift_Counter_X[i][j]);
  //     h_YY1_Ex_Shift_X->SetBinContent(i + 1, j + 1, Ex_YY1_Shift_Counter_X[i][j]);
  //     h_Ex_Shift_Y->SetBinContent(i + 1, j + 1, Ex_Shift_Counter_Y[i][j]);
  //     h_TTT_Ex_Shift_Y->SetBinContent(i + 1, j + 1, Ex_TTT_Shift_Counter_Y[i][j]);
  //     h_YY1_Ex_Shift_Y->SetBinContent(i + 1, j + 1, Ex_YY1_Shift_Counter_Y[i][j]);
  //     h_Ex_Shift_XY->SetBinContent(i + 1, j + 1, Ex_Shift_Counter_XY[i][j]);
  //     h_TTT_Ex_Shift_XY->SetBinContent(i + 1, j + 1, Ex_TTT_Shift_Counter_XY[i][j]);
  //     h_YY1_Ex_Shift_XY->SetBinContent(i + 1, j + 1, Ex_YY1_Shift_Counter_XY[i][j]);
  //   }
  // }

  file->Close();
  ofile_new->cd();
  // hlist->Write();
  tree_new->Write("", TObject::kOverwrite);
  ofile_new->Close();
}

vector<double> f_Coordinates(int ID_Front, int ID_Back)
{

  double TTT_X;
  double TTT_Y;
  double TTT_Z;

  // Some constants of TINA extracted from NPTool simulations.
  double TTT_length = 100.42;                                // Length of a TTT side.
  double TTT_active_length = 97.22;                          // Length of the active region of a TTT side.
  double TTT_strip_width = TTT_active_length / 128.;         // Width of each strip.
  double TTT_target_to_strip = -14.14 - TTT_strip_width / 2; // Distance from target to center of first back strip.

  // Output.
  vector<double> Position;

  // From the front strip we get the X and Y coordinates of the proton.
  if (ID_Front < 128)
  { // If the pixel is in the first TTT.
    TTT_Y = -55;
    TTT_X = (TTT_active_length - TTT_strip_width) / 2 - (TTT_strip_width * (ID_Front));
  }
  if (ID_Front > 127 && ID_Front < 256)
  { // If the pixel is in the second TTT.
    TTT_X = -55;
    TTT_Y = -(TTT_active_length - TTT_strip_width) / 2 + (TTT_strip_width * (ID_Front - 128));
  }
  if (ID_Front > 255 && ID_Front < 384)
  { // If the pixel is in the third TTT.
    TTT_Y = 55;
    TTT_X = -(TTT_active_length - TTT_strip_width) / 2 + (TTT_strip_width * (ID_Front - 256));
  }
  if (ID_Front > 383)
  { // If the pixel is in the fourth TTT.
    TTT_X = 55;
    TTT_Y = (TTT_active_length - TTT_strip_width) / 2 - (TTT_strip_width * (ID_Front - 384));
  }
  // From the back strip we get the Z coordinate.
  if (ID_Back < 128)
  { // If the pixel is in the first TTT.
    TTT_Z = TTT_target_to_strip - (TTT_strip_width * ID_Back);
  }
  if (ID_Back > 127 && ID_Back < 256)
  { // If the pixel is in the second TTT.
    TTT_Z = TTT_target_to_strip - (TTT_strip_width * (ID_Back - 128));
  }
  if (ID_Back > 255 && ID_Back < 384)
  { // If the pixel is in the third TTT.
    TTT_Z = TTT_target_to_strip - (TTT_strip_width * (ID_Back - 256));
  }
  if (ID_Back > 383)
  { // If the pixel is in the fourth TTT.
    TTT_Z = TTT_target_to_strip - (TTT_strip_width * (ID_Back - 384));
  }

  Position.push_back(TTT_X);
  Position.push_back(TTT_Y);
  Position.push_back(TTT_Z);

  return Position;
}
double f_Theta(double X, double Y, double Z, double X_0, double Y_0, double Z_0)
{
  double Theta;
  Theta = acos((X_0 * X + Y_0 * Y + Z_0 * Z) / ((sqrt(pow(X_0, 2) + pow(Y_0, 2) + pow(Z_0, 2)) * sqrt(pow(X, 2) + pow(Y, 2) + pow(Z, 2)))));
  return Theta;
}

double f_Phi(double X, double Y)
{
  double Phi;
  double PI = 3.14159265358979323846;

  if (X < 0)
  {
    Phi = (PI + atan(Y / X));
  }
  if (X > 0)
  {
    if (Y > 0)
    {
      Phi = (atan(Y / X));
    }
    if (Y < 0)
    {
      Phi = (2 * PI + atan(Y / X));
    }
  }
  if (X == 0)
  {
    if (Y > 0)
    {
      Phi = PI / 2.;
    }
    if (Y < 0)
    {
      Phi = -3 * PI / 2.;
    }
  }
  return Phi;
}

double f_Excitation_Energy(double Proton_Energy, double Scattering_Angle, double Beam_Energy, double Incident_Angle)
{

  // Masses of the nuclei.

  double M_49Ca = 45601.906;
  double M_50Ca = 46535.0915;
  double M_51Ca = 47469.8421;
  double M_2H = 1875.61689;
  double M_1H = 938.273993;

  double Qvalue = M_50Ca + M_2H - M_51Ca - M_1H; // Q-value of the 50Ca(d,p)51Ca neutron transfer reaction.
  // double Qvalue = M_49Ca + M_2H - M_50Ca - M_1H; //Q-value of the 49Ca(d,p)50Ca neutron transfer reaction.

  // Total energies of the nuclei.
  double E_50Ca = 0;
  double E_51Ca = 0;
  double E_1H = 0;
  double E_2H = 0;

  // Kinetic energies of the nuclei.
  double T_50Ca = 0;
  double T_51Ca = 0;
  double T_2H = 0;
  double T_1H = 0;

  // Momenta of the nuclei.
  double P_50Ca = 0;
  double Px_50Ca = 0;
  double Py_50Ca = 0;

  double P_51Ca = 0;
  double Px_51Ca = 0;
  double Py_51Ca = 0;

  double P_1H = 0;
  double Px_1H = 0;
  double Py_1H = 0;

  double P_2H = 0;
  double Px_2H = 0;
  double Py_2H = 0;

  double Invariant = 0;
  double Excitation = 0;

  double PI = 3.14159265358979323846;

  T_50Ca = Beam_Energy * 50;
  // T_50Ca = Beam_Energy*49;

  E_50Ca = T_50Ca + M_50Ca;
  P_50Ca = sqrt(pow(E_50Ca, 2) - pow(M_50Ca, 2));

  // E_50Ca = T_50Ca + M_49Ca;
  // P_50Ca = sqrt(pow(E_50Ca,2) - pow(M_49Ca,2));

  // Px_50Ca = P_50Ca*cos(Incident_Angle);
  // Py_50Ca = P_50Ca*sin(Incident_Angle);

  Px_50Ca = P_50Ca;

  // Proton.
  T_1H = Proton_Energy;
  E_1H = T_1H + M_1H;
  P_1H = sqrt(pow(E_1H, 2) - pow(M_1H, 2));
  Px_1H = P_1H * cos(Scattering_Angle);
  Py_1H = P_1H * sin(Scattering_Angle);

  // Excitation energy.
  E_51Ca = E_50Ca + M_2H - E_1H;
  T_51Ca = E_51Ca - M_51Ca;
  Px_51Ca = Px_50Ca - Px_1H;
  Py_51Ca = Py_50Ca - Py_1H;

  Invariant = sqrt(pow(E_51Ca, 2) - pow(Px_51Ca, 2) - pow(Py_51Ca, 2));

  Excitation = Invariant - M_51Ca;
  // cout << T_51Ca << endl;
  return Excitation;
}

double f_list(int a, int b)
{
  FILE *myfile;
  double myvariable;
  string SRPPAC;
  string sector;
  double output;
  int j;
  int k;

  string input = "/u/ddas/software/work/artemis-oedo/output/Analysis/TTT_B_param.txt";
  myfile = fopen(input.c_str(), "r");
  for (j = 0; j < HEIGHT; j++)
  {
    for (k = 0; k < WIDTH; k++)
    {
      fscanf(myfile, "%lf", &myvariable);
      if (j == b && k == a)
      {
        output = myvariable;
      }
    }
  }
  fclose(myfile);
  return output;
}

// double f_PID_Shift(int a)
// {
//   FILE *myfile;
//   double myvariable;
//   string SRPPAC;
//   string sector;
//   double output;
//   int j;
//   int k;
//    str ing input = "/home/jupiter/test/work/artemis-oedo/macro/Analysis/F3/PID_T_Event.txt";
//   myfile=fopen(input.c_str(), "r");
//   for(j = 0; j < 2500; j++)
//     {
//       for (k = 0 ; k < 4; k++)
// 	{
//           fscanf(myfile,"%lf",&myvariable);
//           if ( j == a && k == 2)
//             {
//               output = myvariable;
//             }
// 	}
//     }
// fclose(myfile);
// return output;
// }

// double f_Beam_Shift(int a)
// {
//   FILE *myfile;
//   double myvariable;
//   string SRPPAC;
//   string sector;
//   double output;
//   int j;
//   int k;

//   string input = "/home/jupiter/test/work/artemis-oedo/macro/Analysis/F3/Beam_T_Event.txt";
//   myfile=fopen(input.c_str(), "r");
//   for(j = 0; j < 1000; j++)
//     {
//       for (k = 0 ; k < 4; k++)
// 	{
//           fscanf(myfile,"%lf",&myvariable);
//           if ( j == a && k == 2)
//             {
//               output = myvariable;
//             }
// 	}
//     }
// fclose(myfile);
// return output;
// }

// double f_S1PID_Shift(int a)
// {
//   FILE *myfile;
//   double myvariable;
//   string SRPPAC;
//   string sector;
//   double output;
//   int j;
//   int k;

//   string input = "/home/jupiter/test/work/artemis-oedo/macro/Analysis/F3/S1PID_T_Event.txt";
//   myfile=fopen(input.c_str(), "r");
//   for(j = 0; j < 1250; j++)
//     {
//       for (k = 0 ; k < 4; k++)
// 	{
//           fscanf(myfile,"%lf",&myvariable);
//           if ( j == a && k == 2)
//             {
//               output = myvariable;
//             }
// 	}
//     }
// fclose(myfile);
// return output;
// }

// double f_AQ_Shift(int a)
// {
//   FILE *myfile;
//   double myvariable;
//   string SRPPAC;
//   string sector;
//   double output;
//   int j;
//   int k;

//   string input = "/home/jupiter/test/work/artemis-oedo/macro/Analysis/F3/AQ_Event_2024.txt";
//   myfile=fopen(input.c_str(), "r");
//   for(j = 0; j < 200; j++)
//     {
//       for (k = 0 ; k < 4; k++)
// 	{
//           fscanf(myfile,"%lf",&myvariable);
//           if ( j == a && k == 2)
//             {
//               output = myvariable;
//             }
// 	}
//     }
// fclose(myfile);
// return output;
// }

double f_AQ_2_S1_Y(int a)
{
  FILE *myfile;
  double myvariable;
  string SRPPAC;
  string sector;
  double output;
  int j;
  int k;

  string input = "/u/ddas/software/work/artemis-oedo/output/Analysis/AQ_2_S1_Y.txt";
  myfile = fopen(input.c_str(), "r");
  for (j = 0; j < 319; j++)
  {
    for (k = 0; k < 4; k++)
    {
      fscanf(myfile, "%lf", &myvariable);
      if (j == a && k == 2)
      {
        output = myvariable;
      }
    }
  }
  fclose(myfile);
  return output;
}

double f_AQ_2_S1_X(int a)
{
  FILE *myfile;
  double myvariable;
  string SRPPAC;
  string sector;
  double output;
  int j;
  int k;

  string input = "/u/ddas/software/work/artemis-oedo/output/Analysis/AQ_2_S1_X.txt";
  myfile = fopen(input.c_str(), "r");
  for (j = 0; j < 559; j++)
  {
    for (k = 0; k < 4; k++)
    {
      fscanf(myfile, "%lf", &myvariable);
      if (j == a && k == 2)
      {
        output = myvariable;
      }
    }
  }
  fclose(myfile);
  return output;
}

double f_Energy(double TTT_E, double Theta)
{
  double Energy;
  Energy = 0.124546059382704 * pow(TTT_E * sin(Theta), 4) - 2.11115981536274 * pow(TTT_E * sin(Theta), 3) + 13.6871829996149 * pow(TTT_E * sin(Theta), 2) - 41.3737617194654 * pow(TTT_E * sin(Theta), 1) + 57.0318203606131;
  return Energy;
}
