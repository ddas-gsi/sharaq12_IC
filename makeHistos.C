#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <stdio.h>

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

void makeHistos()
{
    // Physics Run for 2024
    // std::vector<std::string> RUN_list = {"1052", "1053", "1055"}; // After Carlos correction
    std::vector<std::string> RUN_list = {"1056", "1057", "1058", "1059"}; // for Carlos Presentation in OEDO workshop on 24.07.2025. this is with prev way of reduced_chi2
    // std::vector<std::string> RUN_list = {"1052", "1053", "1054", "1055"};

    // Some constants
    const Double_t chi2_mean = 2.344; // Mean of reduced chi2 distribution
    const Double_t chi2_sigma = 1.25; // Standard deviation of reduced chi2 distribution
    // const Double_t chi2_mean = 0.83; // Mean of reduced chi2 distribution
    // const Double_t chi2_sigma = 0.4; // Standard deviation of reduced chi2 distribution

    // We TChain the trees
    TChain *tree_chained = new TChain("tree_new");

    // Create TChain
    for (int fn = 0; fn < RUN_list.size(); fn++)
    {
        cout << "Reading FileNumber: " << RUN_list[fn] << endl;

        string inputfile = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/AoQ/zetCorr/" + RUN_list[fn] + "_Spline_2024.root";
        tree_chained->Add(inputfile.c_str());
    }

    TFile *histFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/AoQ/zetCorr/histSpline_2024.root", "recreate");

    // Read the Graph Cut-Files
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_50Ca20.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_51Sc21.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_49K19.cxx");

    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_50Ca20.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_50Ca19.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_50Ca18.cxx");

    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_51Ca20.cxx");

    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_51Sc21.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_51Sc20.cxx");

    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/50Ca20_20.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/50Ca20_19.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/50Ca20_18.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/51Ca20_20.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/51Sc21_21.cxx");

    TCutG *FE9cut50Ca20 = (TCutG *)gROOT->FindObject("FE9pidCut_50Ca_1053_2024");
    TCutG *FE9cut51Sc21 = (TCutG *)gROOT->FindObject("FE9pidCut_51Sc_1053_2024");
    TCutG *FE9cut49K19 = (TCutG *)gROOT->FindObject("FE9pidCut_49K_1053_2024");

    TCutG *S1cut50Ca20 = (TCutG *)gROOT->FindObject("s1pidCut_50Ca20_1053_2024");
    TCutG *S1cut50Ca19 = (TCutG *)gROOT->FindObject("s1pidCut_50Ca19_1053_2024");
    TCutG *S1cut50Ca18 = (TCutG *)gROOT->FindObject("s1pidCut_50Ca18_1053_2024");

    TCutG *S1cut51Ca20 = (TCutG *)gROOT->FindObject("s1pidCut_51Ca20_1053_2024"); // for 51Ca20

    TCutG *S1cut51Sc21 = (TCutG *)gROOT->FindObject("s1pidCut_51Sc21_1053_2024");
    TCutG *S1cut51Sc20 = (TCutG *)gROOT->FindObject("s1pidCut_51Sc20_1053_2024");

    TCutG *AoQ_51Ca20_20 = (TCutG *)gROOT->FindObject("51Ca20_20");
    TCutG *AoQ_50Ca20_20 = (TCutG *)gROOT->FindObject("50Ca20_20");
    TCutG *AoQ_50Ca20_19 = (TCutG *)gROOT->FindObject("50Ca20_19");
    TCutG *AoQ_50Ca20_18 = (TCutG *)gROOT->FindObject("50Ca20_18");
    TCutG *AoQ_51Sc21_21 = (TCutG *)gROOT->FindObject("51Sc21_21");

    // Variable Required for Tree Reading
    Double_t PID_T_0;        // Time of Flight between F3 and FE9.
    Double_t FE9_X_0;        // Horizontal position at FE9 focal plane.
    Double_t S1PID_T_0;      // Time at S1 (average between anodes).
    Double_t S1_X_0;         // Horizontal position at S1 focal plane.
    Double_t AQ_0;           // AoQ at IC
    Double_t AQ_2_0;         // AoQ at S1 focal plane.
    Double_t Z_a3EBetaS1;    // Z calculated from avgOfFirst3Pads
    Double_t Z_sPeakES1X;    // Z calculated from Spline Peak with S1X correction
    Double_t A_sRangeBetaS1; // A calculated from sRangeBetaS1
    Double_t chi2;           // Chi2 value for the fit

    // Activate the Branches
    tree_chained->SetBranchAddress("PID_T_0", &PID_T_0);
    tree_chained->SetBranchAddress("FE9_X_0", &FE9_X_0);
    tree_chained->SetBranchAddress("S1PID_T_0", &S1PID_T_0);
    tree_chained->SetBranchAddress("S1_X_0", &S1_X_0);
    tree_chained->SetBranchAddress("AQ_0", &AQ_0);
    tree_chained->SetBranchAddress("AQ_2_0", &AQ_2_0);
    tree_chained->SetBranchAddress("Z_a3EBetaS1", &Z_a3EBetaS1);
    tree_chained->SetBranchAddress("Z_sPeakES1X", &Z_sPeakES1X);
    tree_chained->SetBranchAddress("A_sRangeBetaS1", &A_sRangeBetaS1);
    tree_chained->SetBranchAddress("chi2", &chi2);

    // // Crate a TCanvas
    // TCanvas *c1 = new TCanvas("c1", "PID Canvas", 2500, 1800);
    // c1->Divide(5, 5);

    // Create a TList for Histograms
    TList *hlist = new TList();

    // Create Histograms
    // PID for All Isotopes
    // Use the Branch AQ_2_0
    TH2F *Z_vs_AoQ_all_from_avg3Pads_wo_Chi2 = new TH2F("Z_vs_AoQ_all_from_avg3Pads_wo_Chi2", "Z vs AoQ All; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_all_from_avg3Pads_wo_Chi2);

    TH2F *Z_vs_AoQ_all_from_avg3Pads = new TH2F("Z_vs_AoQ_all_from_avg3Pads", "Z vs AoQ All; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_all_from_avg3Pads);

    // Use the Branch AQ_0
    TH2F *Z_vs_AoQ_all_from_avg3Pads_AQ_0 = new TH2F("Z_vs_AoQ_all_from_avg3Pads_AQ_0", "Z vs AoQ All; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_all_from_avg3Pads_AQ_0);

    TH2F *Z_vs_AoQ_all_from_SplinePeak = new TH2F("Z_vs_AoQ_all_from_SplinePeak", "Z vs AoQ All; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_all_from_SplinePeak);

    TH2F *A_vs_AoQ = new TH2F("A_vs_AoQ", "A vs AoQ All; AoQ; A", 1000, 2, 3, 1000, 40, 60);
    hlist->Add(A_vs_AoQ);

    // PID for 50Ca
    TH2F *Z_vs_AoQ_50Ca_from_avg3Pads = new TH2F("Z_vs_AoQ_50Ca_from_avg3Pads", "Z vs AoQ 50Ca; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_50Ca_from_avg3Pads);

    TH2F *Z_vs_AoQ_50Ca_from_SplinePeak = new TH2F("Z_vs_AoQ_50Ca_from_SplinePeak", "Z vs AoQ 50Ca; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_50Ca_from_SplinePeak);

    // PID for 51Ca
    TH2F *Z_vs_AoQ_51Ca_from_avg3Pads = new TH2F("Z_vs_AoQ_51Ca_from_avg3Pads", "Z vs AoQ 51Ca; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_51Ca_from_avg3Pads);

    TH2F *Z_vs_AoQ_51Ca_from_SplinePeak = new TH2F("Z_vs_AoQ_51Ca_from_SplinePeak", "Z vs AoQ 51Ca; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_51Ca_from_SplinePeak);

    // PID for 51Sc
    TH2F *Z_vs_AoQ_51Sc_from_avg3Pads = new TH2F("Z_vs_AoQ_51Sc_from_avg3Pads", "Z vs AoQ 51Sc; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_51Sc_from_avg3Pads);

    TH2F *Z_vs_AoQ_51Sc_from_SplinePeak = new TH2F("Z_vs_AoQ_51Sc_from_SplinePeak", "Z vs AoQ 51Sc; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_51Sc_from_SplinePeak);

    // PID for 49K
    TH2F *Z_vs_AoQ_49K_from_avg3Pads = new TH2F("Z_vs_AoQ_49K_from_avg3Pads", "Z vs AoQ 49K; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_49K_from_avg3Pads);

    TH2F *Z_vs_AoQ_49K_from_SplinePeak = new TH2F("Z_vs_AoQ_49K_from_SplinePeak", "Z vs AoQ 49K; AoQ; Z", 1000, 2, 3, 1000, 10, 30);
    hlist->Add(Z_vs_AoQ_49K_from_SplinePeak);

    // Run the Histogram Creation
    // Get the total Number of Events
    int nEvents = tree_chained->GetEntries();
    cout << "Total Number of Entries: " << nEvents << endl;

    // Loop through Events
    for (int i = 0; i < nEvents; i++)
    {
        if (i % 100000 == 0)
        {
            cout << "Event Number = " << i << endl;
        }

        tree_chained->GetEntry(i);

        // Fill Histogram without Chi2 Cut for All Isotopes
        Z_vs_AoQ_all_from_avg3Pads_wo_Chi2->Fill(AQ_2_0, Z_a3EBetaS1);

        if (abs(chi2 - chi2_mean) < chi2_sigma)
        {
            // Fill Histograms for All Isotopes
            Z_vs_AoQ_all_from_avg3Pads->Fill(AQ_2_0, Z_a3EBetaS1);
            Z_vs_AoQ_all_from_avg3Pads_AQ_0->Fill(AQ_0, Z_a3EBetaS1);
            Z_vs_AoQ_all_from_SplinePeak->Fill(AQ_2_0, Z_sPeakES1X);
            A_vs_AoQ->Fill(AQ_2_0, A_sRangeBetaS1);

            // Fill Histograms for 50Ca
            if (FE9cut50Ca20->IsInside(PID_T_0, FE9_X_0) && abs(chi2 - chi2_mean) < chi2_sigma)
            {
                Z_vs_AoQ_50Ca_from_avg3Pads->Fill(AQ_2_0, Z_a3EBetaS1);
                Z_vs_AoQ_50Ca_from_SplinePeak->Fill(AQ_2_0, Z_sPeakES1X);

                // Fill Histograms for 51Ca
                if (S1cut50Ca20->IsInside(S1PID_T_0, S1_X_0) && (2.52 < AQ_2_0 && AQ_2_0 < 2.6))
                {
                    Z_vs_AoQ_51Ca_from_avg3Pads->Fill(AQ_2_0, Z_a3EBetaS1);
                    Z_vs_AoQ_51Ca_from_SplinePeak->Fill(AQ_2_0, Z_sPeakES1X);
                }
            }

            // Fill Histograms for 51Sc
            else if (FE9cut51Sc21->IsInside(PID_T_0, FE9_X_0) && abs(chi2 - chi2_mean) < chi2_sigma)
            {
                Z_vs_AoQ_51Sc_from_avg3Pads->Fill(AQ_2_0, Z_a3EBetaS1);
                Z_vs_AoQ_51Sc_from_SplinePeak->Fill(AQ_2_0, Z_sPeakES1X);
            }

            // Fill Histograms for 49K
            else if (FE9cut49K19->IsInside(PID_T_0, FE9_X_0) && abs(chi2 - chi2_mean) < chi2_sigma)
            {
                Z_vs_AoQ_49K_from_avg3Pads->Fill(AQ_2_0, Z_a3EBetaS1);
                Z_vs_AoQ_49K_from_SplinePeak->Fill(AQ_2_0, Z_sPeakES1X);
            }
        }
    }

    // Write the Histograms
    histFile->cd();
    // c1->Write();
    hlist->Write("HistogramList", TObject::kSingleKey);
    histFile->Close();
}
