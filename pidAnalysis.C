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

void pidAnalysis()
{
    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/1052_Analysis_2024_new1404.root");
    TTree *tree = (TTree *)file->Get("tree_new");

    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_50Ca20.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_50Ca20.cxx");
    TCutG *FE9cut50Ca20 = (TCutG *)gROOT->FindObject("FE9pidCut_50Ca_1053_2024");
    TCutG *S1cut50Ca20 = (TCutG *)gROOT->FindObject("s1pidCut_50Ca20_1053_2024");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/delta_AQ2_cut_center.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/delta_AQ2_cut_upper.cxx");
    gROOT->ProcessLine(".L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/delta_AQ2_cut_lower.cxx");
    TCutG *delta_AQ2_cut_center = (TCutG *)gROOT->FindObject("delta_AQ2_cut_center");
    TCutG *delta_AQ2_cut_upper = (TCutG *)gROOT->FindObject("delta_AQ2_cut_upper");
    TCutG *delta_AQ2_cut_lower = (TCutG *)gROOT->FindObject("delta_AQ2_cut_lower");

    // Initialize the branches to read
    Double_t PID_T_0;
    Double_t FE9_X_0;
    Double_t S1PID_T_0;
    Double_t S1_X_0;
    Double_t AQ_2_0;
    Double_t peakS1YXAll;
    Double_t deltaAll_stS1Y;
    Double_t rangeAll;
    Double_t peakAll;

    tree->SetBranchAddress("PID_T_0", &PID_T_0);
    tree->SetBranchAddress("FE9_X_0", &FE9_X_0);
    tree->SetBranchAddress("S1PID_T_0", &S1PID_T_0);
    tree->SetBranchAddress("S1_X_0", &S1_X_0);
    tree->SetBranchAddress("AQ_2_0", &AQ_2_0);
    tree->SetBranchAddress("peakS1YXAll", &peakS1YXAll);
    tree->SetBranchAddress("deltaAll_stS1Y", &deltaAll_stS1Y);
    tree->SetBranchAddress("rangeAll", &rangeAll);
    tree->SetBranchAddress("peakAll", &peakAll);

    Double_t IC_E_Cal[30];
    for (int k = 0; k < 30; k++)
    {
        string IC_E_Cal_SetBranch = "IC_E_Cal_" + to_string(k);
        tree->SetBranchAddress(IC_E_Cal_SetBranch.c_str(), &IC_E_Cal[k]); // .c_str() converts any string to CONSTANT
    }

    TCanvas *c1 = new TCanvas("c1", "IC Canvas", 1800, 1350);
    c1->Divide(3, 3);
    TList *hlist = new TList();

    // Define the histograms
    TH2F *h2D_delta_vs_AoQ2 = new TH2F("h2D_50Ca_delta_vs_AoQ2", "50Ca deltaAll_stS1Y vs AoQ2; AoQ2; deltaAll_stS1Y", 1000, 2.3, 3, 1000, 0, 300);
    hlist->Add(h2D_delta_vs_AoQ2);
    TH2F *h2D_IC = new TH2F("h2D_IC", "IC Calibrated Energy; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60);
    hlist->Add(h2D_IC);
    TH2F *h2D_IC_50Ca_Centered = new TH2F("h2D_IC_50Ca_Centered", "IC Calibrated Energy 50Ca Centered; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60);
    hlist->Add(h2D_IC_50Ca_Centered);
    TH2F *h2D_IC_50Ca_Upper = new TH2F("h2D_IC_50Ca_Upper", "IC Calibrated Energy 50Ca Upper; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60);
    hlist->Add(h2D_IC_50Ca_Upper);
    TH2F *h2D_IC_50Ca_Lower = new TH2F("h2D_IC_50Ca_Lower", "IC Calibrated Energy 50Ca Lower; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60);
    hlist->Add(h2D_IC_50Ca_Lower);

    TH1F *h1D_range_50Ca_Centered = new TH1F("h1D_range_50Ca_Centered", "Range 50Ca Centered", 1000, 0, 900);
    hlist->Add(h1D_range_50Ca_Centered);
    TH1F *h1D_range_50Ca_Upper = new TH1F("h1D_range_50Ca_Upper", "Range 50Ca Upper", 1000, 0, 900);
    hlist->Add(h1D_range_50Ca_Upper);
    TH1F *h1D_range_50Ca_Lower = new TH1F("h1D_range_50Ca_Lower", "Range 50Ca Lower", 1000, 0, 900);
    hlist->Add(h1D_range_50Ca_Lower);
    TH1F *h1D_deltaAll_stS1Y_50Ca_Centered = new TH1F("h1D_deltaAll_stS1Y_50Ca_Centered", "deltaAll_stS1Y 50Ca Centered", 1000, 0, 300);
    hlist->Add(h1D_deltaAll_stS1Y_50Ca_Centered);
    TH1F *h1D_deltaAll_stS1Y_50Ca_Upper = new TH1F("h1D_deltaAll_stS1Y_50Ca_Upper", "deltaAll_stS1Y 50Ca Upper", 1000, 0, 300);
    hlist->Add(h1D_deltaAll_stS1Y_50Ca_Upper);
    TH1F *h1D_deltaAll_stS1Y_50Ca_Lower = new TH1F("h1D_deltaAll_stS1Y_50Ca_Lower", "deltaAll_stS1Y 50Ca Lower", 1000, 0, 300);
    hlist->Add(h1D_deltaAll_stS1Y_50Ca_Lower);
    TH1F *h1D_peakAll_50Ca_Centered = new TH1F("h1D_peakAll_50Ca_Centered", "peakAll 50Ca Centered", 1000, 0, 70);
    hlist->Add(h1D_peakAll_50Ca_Centered);
    TH1F *h1D_peakAll_50Ca_Upper = new TH1F("h1D_peakAll_50Ca_Upper", "peakAll 50Ca Upper", 1000, 0, 70);
    hlist->Add(h1D_peakAll_50Ca_Upper);
    TH1F *h1D_peakAll_50Ca_Lower = new TH1F("h1D_peakAll_50Ca_Lower", "peakAll 50Ca Lower", 1000, 0, 70);
    hlist->Add(h1D_peakAll_50Ca_Lower);
    TH1F *h1D_peakAllS1YX_50Ca_Centered = new TH1F("h1D_peakAllS1YX_50Ca_Centered", "peakAllS1YX 50Ca Centered", 1000, 0, 70);
    hlist->Add(h1D_peakAllS1YX_50Ca_Centered);
    TH1F *h1D_peakAllS1YX_50Ca_Upper = new TH1F("h1D_peakAllS1YX_50Ca_Upper", "peakAllS1YX 50Ca Upper", 1000, 0, 70);
    hlist->Add(h1D_peakAllS1YX_50Ca_Upper);
    TH1F *h1D_peakAllS1YX_50Ca_Lower = new TH1F("h1D_peakAllS1YX_50Ca_Lower", "peakAllS1YX 50Ca Lower", 1000, 0, 70);
    hlist->Add(h1D_peakAllS1YX_50Ca_Lower);

    int entries = tree->GetEntries();
    cout << "Number of Entries: " << entries << endl;
    int count = 0;
    int countBad = 0;
    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);

        if (i % 100000 == 0)
        {
            cout << "Event Number = " << i << endl;
        }

        if (FE9cut50Ca20->IsInside(PID_T_0, FE9_X_0))
        {
            h2D_delta_vs_AoQ2->Fill(AQ_2_0, deltaAll_stS1Y);

            for (int j = 0; j < 30; j++)
            {
                h2D_IC->Fill(j, IC_E_Cal[j]);
            }

            if (delta_AQ2_cut_center->IsInside(AQ_2_0, deltaAll_stS1Y))
            {
                h1D_range_50Ca_Centered->Fill(rangeAll);
                h1D_deltaAll_stS1Y_50Ca_Centered->Fill(deltaAll_stS1Y);
                // h1D_peakAll_50Ca_Centered->Fill(peakAll);
                // h1D_peakAllS1YX_50Ca_Centered->Fill(peakS1YXAll);
                for (int j = 0; j < 30; j++)
                {
                    h2D_IC_50Ca_Centered->Fill(j, IC_E_Cal[j]);
                }
            }
            else if (delta_AQ2_cut_upper->IsInside(AQ_2_0, deltaAll_stS1Y))
            {
                h1D_range_50Ca_Upper->Fill(rangeAll);
                h1D_deltaAll_stS1Y_50Ca_Upper->Fill(deltaAll_stS1Y);
                // h1D_peakAll_50Ca_Upper->Fill(peakAll);
                // h1D_peakAllS1YX_50Ca_Upper->Fill(peakS1YXAll);
                for (int j = 0; j < 30; j++)
                {
                    h2D_IC_50Ca_Upper->Fill(j, IC_E_Cal[j]);
                }
            }
            else if (delta_AQ2_cut_lower->IsInside(AQ_2_0, deltaAll_stS1Y))
            {
                h1D_range_50Ca_Lower->Fill(rangeAll);
                h1D_deltaAll_stS1Y_50Ca_Lower->Fill(deltaAll_stS1Y);
                // h1D_peakAll_50Ca_Lower->Fill(peakAll);
                // h1D_peakAllS1YX_50Ca_Lower->Fill(peakS1YXAll);
                for (int j = 0; j < 30; j++)
                {
                    h2D_IC_50Ca_Lower->Fill(j, IC_E_Cal[j]);
                }
            }
        }
    }

    c1->cd(1);
    h2D_delta_vs_AoQ2->Draw("colz");
    delta_AQ2_cut_center->Draw("same");
    delta_AQ2_cut_center->SetLineColor(kRed);
    delta_AQ2_cut_center->SetLineWidth(2);
    delta_AQ2_cut_center->SetLineStyle(2);
    delta_AQ2_cut_upper->Draw("same");
    delta_AQ2_cut_upper->SetLineColor(kRed);
    delta_AQ2_cut_upper->SetLineWidth(2);
    delta_AQ2_cut_upper->SetLineStyle(2);
    delta_AQ2_cut_lower->Draw("same");
    delta_AQ2_cut_lower->SetLineColor(kRed);
    delta_AQ2_cut_lower->SetLineWidth(2);
    delta_AQ2_cut_lower->SetLineStyle(2);

    c1->cd(2);
    h2D_IC->Draw("colz");
    c1->cd(3);
    h2D_IC_50Ca_Centered->Draw("colz");
    c1->cd(4);
    h2D_IC_50Ca_Upper->Draw("colz");
    c1->cd(5);
    h2D_IC_50Ca_Lower->Draw("colz");
    c1->cd(6);
    h1D_range_50Ca_Centered->Draw();
    // h1D_range_50Ca_Centered->Fit("gaus");
    c1->cd(7);
    h1D_range_50Ca_Upper->Draw();
    // h1D_range_50Ca_Upper->Fit("gaus");
    c1->cd(8);
    h1D_range_50Ca_Lower->Draw();
    // h1D_range_50Ca_Lower->Fit("gaus");

    c1->Update();

    // Save the histograms to a file
    TFile *histFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/pidAnaHist.root", "RECREATE");
    histFile->cd();
    // c1->Write();
    hlist->Write();
    histFile->Close();
    file->Close();
}