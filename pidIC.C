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

void pidIC()
{
	int frun = 0;
	int lrun = 2;

	// Physics Run for 2024
	// std::vector<std::string> RUN_list = {"1052", "1053", "1054", "1055", "1056", "1057", "1058"};
	std::vector<std::string> RUN_list = {"1052", "1053", "1054", "1055"};

	// We TChain the trees
	TChain *tree_chained = new TChain("tree_new");

	// for (int fn = frun; fn < lrun; fn++)
	for (int fn = 0; fn < RUN_list.size(); fn++)
	{
		cout << "Reading FileNumber: " << RUN_list[fn] << endl;

		// string inputfile = "/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/" + RUN_list[fn] + "_Analysis_2024.root";
		string inputfile = "/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/" + RUN_list[fn] + "_Analysis_2024_new1404.root";
		tree_chained->Add(inputfile.c_str());
	}

	TFile *histFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/histograms_2024_new1404.root", "recreate");

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

	// Check the cut files:
	if (!FE9cut50Ca20)
	{
		std::cerr << "Error: Could not load the FE9cut50Ca20 graphical cut!" << std::endl;
		return; // Exit with error code
	}
	if (!S1cut50Ca20)
	{
		std::cerr << "Error: Could not load the S1cut50Ca20 graphical cut!" << std::endl;
		return; // Exit with error code
	}

	// Variable Required for Tree Reading
	Double_t PID_T_0;	// Time of Flight between F3 and FE9.
	Double_t FE9_X_0;	// Horizontal position at FE9 focal plane.
	Double_t S1PID_T_0; // Time at S1 (average between anodes).
	Double_t S1_X_0;	// Horizontal position at S1 focal plane.
	Double_t AQ_0;		// AoQ at IC
	Double_t AQ_2_0;	// AoQ at S1 focal plane.

	Double_t peakS1YAll;
	Double_t peakS1YXAll;
	Double_t deltaAll_st;
	Double_t deltaAll_stS1Y;

	// Activate the Branches
	tree_chained->SetBranchAddress("PID_T_0", &PID_T_0);
	tree_chained->SetBranchAddress("FE9_X_0", &FE9_X_0);
	tree_chained->SetBranchAddress("S1PID_T_0", &S1PID_T_0);
	tree_chained->SetBranchAddress("S1_X_0", &S1_X_0);
	tree_chained->SetBranchAddress("AQ_0", &AQ_0);
	tree_chained->SetBranchAddress("AQ_2_0", &AQ_2_0);

	// tree_chained->SetBranchStatus("peakS1YAll", 1);
	tree_chained->SetBranchAddress("peakS1YAll", &peakS1YAll);
	// tree_chained->SetBranchStatus("peakS1YXAll", 1);
	tree_chained->SetBranchAddress("peakS1YXAll", &peakS1YXAll);
	// tree_chained->SetBranchStatus("deltaAll_st", 1);
	tree_chained->SetBranchAddress("deltaAll_st", &deltaAll_st);
	// tree_chained->SetBranchStatus("deltaAll_stS1Y", 1);
	tree_chained->SetBranchAddress("deltaAll_stS1Y", &deltaAll_stS1Y);

	// Crate a TCanvas
	TCanvas *c1 = new TCanvas("c1", "PID Canvas", 2500, 1800);
	c1->Divide(5, 5);

	// Create a TList for Histograms
	TList *hlist = new TList();

	// Add Graphical Cuts to the TList
	hlist->Add(FE9cut50Ca20);
	hlist->Add(S1cut50Ca20);
	hlist->Add(S1cut50Ca19);
	hlist->Add(S1cut50Ca18);
	hlist->Add(S1cut51Ca20);
	hlist->Add(FE9cut51Sc21);
	hlist->Add(S1cut51Sc21);
	hlist->Add(S1cut51Sc20);
	hlist->Add(FE9cut49K19);
	hlist->Add(AoQ_51Ca20_20);
	hlist->Add(AoQ_50Ca20_20);
	hlist->Add(AoQ_50Ca20_19);
	hlist->Add(AoQ_50Ca20_18);
	hlist->Add(AoQ_51Sc21_21);

	// Create a Histogram
	TH2F *h2D_FE9_PID = new TH2F("FE9_PID", "FE9_X_0:PID_T_0; TOF (ns); X Position", 1000, 1050, 1150, 1000, -50, 50);
	hlist->Add(h2D_FE9_PID);
	TH2F *h2D_S1_PID = new TH2F("S1_PID", "S1_X_0:S1PID_T_0; TOF (ns); X Position", 1000, 600, 700, 1000, -200, 200);
	hlist->Add(h2D_S1_PID);
	TH2F *h2D_S1_PID_50Ca = new TH2F("h2D_S1_PID_50Ca", "S1_X_0:S1PID_T_0 for 50Ca with FE9 filter; TOF (ns); X Position", 1000, 600, 700, 1000, -200, 200);
	hlist->Add(h2D_S1_PID_50Ca);
	TH2F *h2D_S1_PID_51Ca20 = new TH2F("h2D_S1_PID_51Ca20", "S1_X_0:S1PID_T_0 for 51Ca20; TOF (ns); X Position", 1000, 600, 700, 1000, -200, 200);
	hlist->Add(h2D_S1_PID_51Ca20);
	TH2F *h2D_S1_PID_51Sc = new TH2F("h2D_S1_PID_51Sc", "S1_X_0:S1PID_T_0 for 51Sc with FE9 filter; TOF (ns); X Position", 1000, 600, 700, 1000, -200, 200);
	hlist->Add(h2D_S1_PID_51Sc);
	TH2F *h2D_S1_PID_49K = new TH2F("h2D_S1_PID_49K", "S1_X_0:S1PID_T_0 for 49K with FE9 filter; TOF (ns); X Position", 1000, 600, 700, 1000, -200, 200);
	hlist->Add(h2D_S1_PID_49K);
	TH2F *h2D_AoQ = new TH2F("AoQ", "AQ_0 vs AQ_2_0; AQ_2_0; AQ_0", 1000, 2, 3, 1000, 2, 3);
	hlist->Add(h2D_AoQ);

	// pid plots for all isotopes
	TH2F *h2D_delta_vs_AoQ_2 = new TH2F("h2D_delta_vs_AoQ_2", "delta vs AoQ_2; AoQ_2; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_delta_vs_AoQ_2);
	TH2F *h2D_delta_vs_AoQ = new TH2F("h2D_delta_vs_AoQ", "delta vs AoQ; AoQ; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_delta_vs_AoQ);

	// For 50Ca, 51Ca, 51Sc and 49K pid plot for different charge states
	TH2F *h2D_50Ca_delta_vs_AoQ_2 = new TH2F("h2D_50Ca_delta_vs_AoQ_2", "50Ca delta vs AoQ_2; AoQ_2; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_50Ca_delta_vs_AoQ_2);
	TH2F *h2D_50Ca_delta_vs_AoQ = new TH2F("h2D_50Ca_delta_vs_AoQ", "50Ca delta vs AoQ; AoQ; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_50Ca_delta_vs_AoQ);
	TH2F *h2D_51Sc_delta_vs_AoQ_2 = new TH2F("h2D_51Sc_delta_vs_AoQ_2", "51Sc delta vs AoQ_2; AoQ_2; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_51Sc_delta_vs_AoQ_2);
	TH2F *h2D_51Sc_delta_vs_AoQ = new TH2F("h2D_51Sc_delta_vs_AoQ", "51Sc delta vs AoQ; AoQ; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_51Sc_delta_vs_AoQ);
	TH2F *h2D_49K_delta_vs_AoQ_2 = new TH2F("h2D_49K_delta_vs_AoQ_2", "49K delta vs AoQ_2; AoQ_2; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_49K_delta_vs_AoQ_2);
	TH2F *h2D_49K_delta_vs_AoQ = new TH2F("h2D_49K_delta_vs_AoQ", "49K delta vs AoQ; AoQ; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_49K_delta_vs_AoQ);
	TH2F *h2D_51Ca_delta_vs_AoQ_2 = new TH2F("h2D_51Ca_delta_vs_AoQ_2", "51Ca delta vs AoQ_2 with AoQ-FE9 filter; AoQ_2; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_51Ca_delta_vs_AoQ_2);
	TH2F *h2D_51Ca_delta_vs_AoQ = new TH2F("h2D_51Ca_delta_vs_AoQ", "51Ca delta vs AoQ with AoQ-FE9 filter; AoQ; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_51Ca_delta_vs_AoQ);
	TH2F *h2D_51Ca_delta_vs_AoQ_2_wS1 = new TH2F("h2D_51Ca_delta_vs_AoQ_2_wS1", "51Ca delta vs AoQ_2 with AoQ-FE9-S1 filter; AoQ_2; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_51Ca_delta_vs_AoQ_2_wS1);
	TH2F *h2D_51Ca_delta_vs_AoQ_wS1 = new TH2F("h2D_51Ca_delta_vs_AoQ_wS1", "51Ca delta vs AoQ with AoQ-FE9-S1 filter; AoQ; delta", 1000, 2, 3, 1000, 0, 300);
	hlist->Add(h2D_51Ca_delta_vs_AoQ_wS1);

	TH2F *h2D_pidAll = new TH2F("h2D_pidAll", "pidAll; delta; peak", 1000, 0, 300, 1000, 0, 70);
	hlist->Add(h2D_pidAll);

	TH2F *h2D_pidCa20 = new TH2F("h2D_pidCa20", "pidCa20+ with FE9-S1 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with FE9-S1 filter
	hlist->Add(h2D_pidCa20);
	TH2F *h2D_pidCa20_AoQ = new TH2F("h2D_pidCa20_AoQ", "pidCa20+ with AoQ-FE9-S1 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with AoQ-FE9-S1 filter
	hlist->Add(h2D_pidCa20_AoQ);

	TH2F *h2D_pidCa19 = new TH2F("h2D_pidCa19", "pidCa19+ with FE9-S1 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with FE9-S1 filter
	hlist->Add(h2D_pidCa19);
	TH2F *h2D_pidCa19_AoQ = new TH2F("h2D_pidCa19_AoQ", "pidCa19+ with AoQ-FE9-S1 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with AoQ-FE9-S1 filter
	hlist->Add(h2D_pidCa19_AoQ);

	TH2F *h2D_pidCa18 = new TH2F("h2D_pidCa18", "pidCa18+ with FE9-S1 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with FE9-S1 filter
	hlist->Add(h2D_pidCa18);
	TH2F *h2D_pidCa18_AoQ = new TH2F("h2D_pidCa18_AoQ", "pidCa18+ with AoQ-FE9-S1 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with AoQ-FE9-S1 filter
	hlist->Add(h2D_pidCa18_AoQ);

	TH2F *h2D_pid51Ca20_AoQ = new TH2F("h2D_pid51Ca20_AoQ", "pid 51Ca20+ with AoQ-FE9 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with AoQ-FE9 filter
	hlist->Add(h2D_pid51Ca20_AoQ);
	TH2F *h2D_pid51Ca20_AoQ_wS1 = new TH2F("h2D_pid51Ca20_AoQ_wS1", "pid 51Ca20+ with AoQ-FE9-S1 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with AoQ-FE9 filter
	hlist->Add(h2D_pid51Ca20_AoQ_wS1);

	TH2F *h2D_pidScAll = new TH2F("h2D_pidScAll", "pidScAll; delta; peak", 1000, 0, 300, 1000, 0, 70); // with FE9 filter
	hlist->Add(h2D_pidScAll);
	TH2F *h2D_pidScAll_AoQ = new TH2F("h2D_pidScAll_AoQ", "pidScAll with AoQ-FE9 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with AoQ-FE9 filter
	hlist->Add(h2D_pidScAll_AoQ);

	TH2F *h2D_pidKAll = new TH2F("h2D_pidKAll", "pidKAll with FE9 filter; delta; peak", 1000, 0, 300, 1000, 0, 70); // with FE9 filter
	hlist->Add(h2D_pidKAll);

	TH1F *h1D_delta50Ca20 = new TH1F("h1D_delta50Ca20", "delta50Ca20; delta", 1000, 0, 300);
	hlist->Add(h1D_delta50Ca20);
	TH1F *h1D_peak50Ca20 = new TH1F("h1D_peak50Ca20", "peak50Ca20; peak", 1000, 0, 70);
	hlist->Add(h1D_peak50Ca20);
	TH1F *h1D_delta50Ca19 = new TH1F("h1D_delta50Ca19", "delta50Ca19; delta", 1000, 0, 300);
	hlist->Add(h1D_delta50Ca19);
	TH1F *h1D_peak50Ca19 = new TH1F("h1D_peak50Ca19", "peak50Ca19; peak", 1000, 0, 70);
	hlist->Add(h1D_peak50Ca19);
	TH1F *h1D_delta50Ca18 = new TH1F("h1D_delta50Ca18", "delta50Ca18; delta", 1000, 0, 300);
	hlist->Add(h1D_delta50Ca18);
	TH1F *h1D_peak50Ca18 = new TH1F("h1D_peak50Ca18", "peak50Ca18; peak", 1000, 0, 70);
	hlist->Add(h1D_peak50Ca18);

	TH1F *h1D_delta51Ca20 = new TH1F("h1D_delta51Ca20", "delta51Ca20; delta", 1000, 0, 300);
	hlist->Add(h1D_delta51Ca20);
	TH1F *h1D_peak51Ca20 = new TH1F("h1D_peak51Ca20", "peak51Ca20; peak", 1000, 0, 70);
	hlist->Add(h1D_peak51Ca20);

	TH1F *h1D_delta51Sc21 = new TH1F("h1D_delta51Sc21", "delta51Sc21; delta", 1000, 0, 300);
	hlist->Add(h1D_delta51Sc21);
	TH1F *h1D_peak51Sc21 = new TH1F("h1D_peak51Sc21", "peak51Sc21; peak", 1000, 0, 70);
	hlist->Add(h1D_peak51Sc21);
	TH1F *h1D_delta51Sc20 = new TH1F("h1D_delta51Sc20", "delta51Sc20; delta", 1000, 0, 300);
	hlist->Add(h1D_delta51Sc20);
	TH1F *h1D_peak51Sc20 = new TH1F("h1D_peak51Sc20", "peak51Sc20; peak", 1000, 0, 70);
	hlist->Add(h1D_peak51Sc20);

	TH1F *h1D_delta49K19 = new TH1F("h1D_delta49K19", "delta49K19; delta", 1000, 0, 300);
	hlist->Add(h1D_delta49K19);
	TH1F *h1D_peak49K19 = new TH1F("h1D_peak49K19", "peak49K19; peak", 1000, 0, 70);
	hlist->Add(h1D_peak49K19);

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

		// cout << "peakS1YAll: " << peakS1YAll << "peakS1YXAll: " << peakS1YXAll << "deltaAll_st: " << deltaAll_st << "deltaAll_stS1Y: " << deltaAll_stS1Y << endl;

		h2D_FE9_PID->Fill(PID_T_0, FE9_X_0);
		h2D_S1_PID->Fill(S1PID_T_0, S1_X_0);
		h2D_AoQ->Fill(AQ_2_0, AQ_0);

		h2D_pidAll->Fill(deltaAll_stS1Y, peakS1YXAll);

		h2D_delta_vs_AoQ_2->Fill(AQ_2_0, deltaAll_stS1Y);
		h2D_delta_vs_AoQ->Fill(AQ_0, deltaAll_stS1Y);

		if (FE9cut50Ca20->IsInside(PID_T_0, FE9_X_0))
		{
			h2D_S1_PID_50Ca->Fill(S1PID_T_0, S1_X_0);
			h2D_50Ca_delta_vs_AoQ_2->Fill(AQ_2_0, deltaAll_stS1Y);
			h2D_50Ca_delta_vs_AoQ->Fill(AQ_0, deltaAll_stS1Y);

			if (S1cut50Ca20->IsInside(S1PID_T_0, S1_X_0))
			{
				h2D_pidCa20->Fill(deltaAll_stS1Y, peakS1YXAll);
				h1D_delta50Ca20->Fill(deltaAll_stS1Y);
				h1D_peak50Ca20->Fill(peakS1YXAll);

				if (AoQ_50Ca20_20->IsInside(AQ_2_0, AQ_0))
				{
					h2D_pidCa20_AoQ->Fill(deltaAll_stS1Y, peakS1YXAll);
				}
			}

			else if (S1cut50Ca19->IsInside(S1PID_T_0, S1_X_0))
			{
				h2D_pidCa19->Fill(deltaAll_stS1Y, peakS1YXAll);
				h1D_delta50Ca19->Fill(deltaAll_stS1Y);
				h1D_peak50Ca19->Fill(peakS1YXAll);

				if (AoQ_50Ca20_19->IsInside(AQ_2_0, AQ_0))
				{
					h2D_pidCa19_AoQ->Fill(deltaAll_stS1Y, peakS1YXAll);
				}
			}

			else if (S1cut50Ca18->IsInside(S1PID_T_0, S1_X_0))
			{
				h2D_pidCa18->Fill(deltaAll_stS1Y, peakS1YXAll);
				h1D_delta50Ca18->Fill(deltaAll_stS1Y);
				h1D_peak50Ca18->Fill(peakS1YXAll);

				if (AoQ_50Ca20_18->IsInside(AQ_2_0, AQ_0))
				{
					h2D_pidCa18_AoQ->Fill(deltaAll_stS1Y, peakS1YXAll);
				}
			}

			if (AoQ_51Ca20_20->IsInside(AQ_2_0, AQ_0))
			{
				h2D_pid51Ca20_AoQ->Fill(deltaAll_stS1Y, peakS1YXAll);
				h2D_51Ca_delta_vs_AoQ_2->Fill(AQ_2_0, deltaAll_stS1Y);
				h2D_51Ca_delta_vs_AoQ->Fill(AQ_0, deltaAll_stS1Y);
				h2D_S1_PID_51Ca20->Fill(S1PID_T_0, S1_X_0);
				h1D_delta51Ca20->Fill(deltaAll_stS1Y);
				h1D_peak51Ca20->Fill(peakS1YXAll);

				if (S1cut51Ca20->IsInside(S1PID_T_0, S1_X_0))
				{
					h2D_51Ca_delta_vs_AoQ_2_wS1->Fill(AQ_2_0, deltaAll_stS1Y);
					h2D_51Ca_delta_vs_AoQ_wS1->Fill(AQ_0, deltaAll_stS1Y);
					h2D_pid51Ca20_AoQ_wS1->Fill(deltaAll_stS1Y, peakS1YXAll);
					// h2D_S1_PID_51Ca20->Fill(S1PID_T_0, S1_X_0);
				}
			}
		}

		else if (FE9cut51Sc21->IsInside(PID_T_0, FE9_X_0))
		{
			h2D_S1_PID_51Sc->Fill(S1PID_T_0, S1_X_0);
			h2D_51Sc_delta_vs_AoQ_2->Fill(AQ_2_0, deltaAll_stS1Y);
			h2D_51Sc_delta_vs_AoQ->Fill(AQ_0, deltaAll_stS1Y);

			if (S1cut51Sc21->IsInside(S1PID_T_0, S1_X_0))
			{
				h2D_pidScAll->Fill(deltaAll_stS1Y, peakS1YXAll);
				h1D_delta51Sc21->Fill(deltaAll_stS1Y);
				h1D_peak51Sc21->Fill(peakS1YXAll);

				if (AoQ_51Sc21_21->IsInside(AQ_2_0, AQ_0))
				{
					h2D_pidScAll_AoQ->Fill(deltaAll_stS1Y, peakS1YXAll);
				}
			}

			else if (S1cut51Sc20->IsInside(S1PID_T_0, S1_X_0))
			{
				h1D_delta51Sc20->Fill(deltaAll_stS1Y);
				h1D_peak51Sc20->Fill(peakS1YXAll);
			}
		}
		else if (FE9cut49K19->IsInside(PID_T_0, FE9_X_0))
		{
			h2D_S1_PID_49K->Fill(S1PID_T_0, S1_X_0);
			h2D_pidKAll->Fill(deltaAll_stS1Y, peakS1YXAll);
			h1D_delta49K19->Fill(deltaAll_stS1Y);
			h1D_peak49K19->Fill(peakS1YXAll);
			h2D_49K_delta_vs_AoQ_2->Fill(AQ_2_0, deltaAll_stS1Y);
			h2D_49K_delta_vs_AoQ->Fill(AQ_0, deltaAll_stS1Y);
		}
	}

	// Draw the Histogram
	c1->cd(1);
	h2D_FE9_PID->Draw("colz");
	FE9cut50Ca20->Draw("same");
	FE9cut50Ca20->SetLineColor(kRed);
	FE9cut50Ca20->SetLineWidth(2);
	FE9cut50Ca20->SetLineStyle(2);
	FE9cut51Sc21->Draw("same");
	FE9cut51Sc21->SetLineColor(kRed);
	FE9cut51Sc21->SetLineWidth(2);
	FE9cut51Sc21->SetLineStyle(2);
	FE9cut49K19->Draw("same");
	FE9cut49K19->SetLineColor(kRed);
	FE9cut49K19->SetLineWidth(2);
	FE9cut49K19->SetLineStyle(2);

	c1->cd(2);
	h2D_S1_PID->Draw("colz");
	S1cut50Ca20->Draw("same");
	S1cut50Ca20->SetLineColor(kRed);
	S1cut50Ca20->SetLineWidth(2);
	S1cut50Ca20->SetLineStyle(2);
	S1cut50Ca19->Draw("same");
	S1cut50Ca19->SetLineColor(kRed);
	S1cut50Ca19->SetLineWidth(2);
	S1cut50Ca19->SetLineStyle(2);
	S1cut50Ca18->Draw("same");
	S1cut50Ca18->SetLineColor(kRed);
	S1cut50Ca18->SetLineWidth(2);
	S1cut50Ca18->SetLineStyle(2);

	c1->cd(3);
	h2D_AoQ->Draw("colz");
	AoQ_50Ca20_20->Draw("same");
	AoQ_50Ca20_20->SetLineColor(kRed);
	AoQ_50Ca20_20->SetLineWidth(2);
	AoQ_50Ca20_19->Draw("same");
	AoQ_50Ca20_19->SetLineColor(kRed);
	AoQ_50Ca20_19->SetLineWidth(2);
	AoQ_50Ca20_18->Draw("same");
	AoQ_50Ca20_18->SetLineColor(kRed);
	AoQ_50Ca20_18->SetLineWidth(2);
	AoQ_51Ca20_20->Draw("same");
	AoQ_51Ca20_20->SetLineColor(kRed);
	AoQ_51Ca20_20->SetLineWidth(2);
	AoQ_51Sc21_21->Draw("same");
	AoQ_51Sc21_21->SetLineColor(kRed);
	AoQ_51Sc21_21->SetLineWidth(2);

	c1->cd(4);
	h2D_delta_vs_AoQ_2->Draw("colz");

	c1->cd(5);
	h2D_delta_vs_AoQ->Draw("colz");

	c1->cd(6);
	h2D_pidCa20->Draw("colz");

	c1->cd(7);
	h2D_pidCa20_AoQ->Draw("colz");

	c1->cd(8);
	h2D_pid51Ca20_AoQ->Draw("colz");

	c1->cd(9);
	h2D_pidCa19->Draw("colz");

	c1->cd(10);
	h2D_pidCa19_AoQ->Draw("colz");

	c1->cd(11);
	h2D_pidCa18->Draw("colz");

	c1->cd(12);
	h2D_pidCa18_AoQ->Draw("colz");

	c1->cd(13);
	h2D_pidScAll->Draw("colz");

	c1->cd(14);
	h2D_pidScAll_AoQ->Draw("colz");

	c1->cd(15);
	h2D_pidKAll->Draw("colz");

	c1->cd(16);
	h2D_pidAll->Draw("colz");

	c1->cd(17);
	h1D_delta50Ca20->Draw();
	h1D_delta51Sc21->Draw("same");
	h1D_delta51Sc21->SetLineColor(kRed);
	h1D_delta49K19->Draw("same");
	h1D_delta49K19->SetLineColor(kGreen);

	c1->cd(18);
	h1D_peak50Ca20->Draw();
	h1D_peak51Sc21->Draw("same");
	h1D_peak51Sc21->SetLineColor(kRed);
	h1D_peak49K19->Draw("same");
	h1D_peak49K19->SetLineColor(kGreen);

	c1->cd(19);
	h2D_50Ca_delta_vs_AoQ_2->Draw("colz");

	c1->cd(20);
	h2D_51Sc_delta_vs_AoQ_2->Draw("colz");

	c1->cd(21);
	h2D_49K_delta_vs_AoQ_2->Draw("colz");

	c1->cd(22);
	h2D_pid51Ca20_AoQ_wS1->Draw("colz");

	c1->cd(23);
	h2D_51Ca_delta_vs_AoQ_2->Draw("colz");

	c1->cd(24);
	h2D_51Ca_delta_vs_AoQ_2_wS1->Draw("colz");

	c1->cd(25);
	h2D_S1_PID_51Ca20->Draw("colz");
	S1cut51Ca20->Draw("same");
	S1cut51Ca20->SetLineColor(kRed);
	S1cut51Ca20->SetLineWidth(2);
	S1cut51Ca20->SetLineStyle(2);

	c1->Update();
	c1->Draw();

	// c1->WaitPrimitive();

	// Write the Histograms
	histFile->cd();
	c1->Write();
	hlist->Write("HistogramList", TObject::kSingleKey);
	histFile->Close();
}
