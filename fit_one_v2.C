#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"

// // Energy Loss for BeamEnergy 15.5 MeV
// vector<double> EnergyLoss15_5MeV = {19.307, 19.777, 20.276, 20.82, 21.413, 22.06, 22.77,
//                                     23.554, 24.425, 25.4, 26.502, 27.759, 29.211, 30.905,
//                                     32.892, 35.247, 37.76, 35.479, 12.255, 0, 0, 0, 0, 0,
//                                     0, 0, 0, 0, 0, 0};

vector<double> Segment = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                          17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

const int n = 30;

void fit_one_v2()
{

    // ---------------------------------------------------------------------
    // Load LISE++ Simulation Data
    // ---------------------------------------------------------------------
    vector<double> BeamEnergies = {13.0, 15.0, 15.5, 16.0, 17.5};

    // std::vector<int> ICSegment13_0MeV, ICSegment15_0MeV, ICSegment15_5MeV, ICSegment16_0MeV, ICSegment17_5MeV;
    // std::vector<double> EnergyLoss13_0MeV, EnergyLoss15_0MeV, EnergyLoss15_5MeV, EnergyLoss16_0MeV, EnergyLoss17_5MeV;

    std::map<std::string, std::vector<int>> ICSegmentMap;
    std::map<std::string, std::vector<float>> EnergyLossMap;

    std::vector<std::vector<int>> ICSegment2DVec;
    std::vector<std::vector<float>> EnergyLoss2DVec;

    std::vector<std::string> LISEFiles;

    for (int i = 0; i < BeamEnergies.size(); i++)
    {
        std::string liseFile = Form("%0.1fMeV.txt", BeamEnergies[i]); // + "MeV.txt";
        LISEFiles.push_back(liseFile.c_str());
    }

    // for (int i = 0; i < BeamEnergies.size(); i++)
    // {
    //     std::string ICSegment = Form("ICSegment%.0f_%0.0f", BeamEnergies[i], (BeamEnergies[i] - int(BeamEnergies[i])) * 10);
    //     cout << ICSegment.c_str() << endl;
    // }

    for (int k = 0; k < BeamEnergies.size(); k++)
    {
        std::string liseFile = "/u/ddas/software/work/artemis-oedo/output/Analysis/" + LISEFiles[k];
        ifstream inFile(liseFile.c_str());

        // // define the EnergyLoss and ICSegment vectors
        // std::string ICSegmentStr = Form("ICSegment%.0f_%0.0f", BeamEnergies[k], (BeamEnergies[k] - int(BeamEnergies[k])) * 10);
        // std::string EnergyLossStr = Form("EnergyLoss%.0f_%0.0f", BeamEnergies[k], (BeamEnergies[k] - int(BeamEnergies[k])) * 10);

        // // Create a new vector and add it to the map with the generated name as the key
        // ICSegmentMap[ICSegmentStr] = std::vector<double>();   // Empty vector
        // EnergyLossMap[EnergyLossStr] = std::vector<double>(); // Empty vector
        // ICSegmentMap.emplace(ICSegmentStr, std::vector<double>());
        // EnergyLossMap.emplace(EnergyLossStr, std::vector<double>());

        std::vector<int> ICSegmentsVector;
        std::vector<float> EnergyLossVector;

        if (inFile.is_open())
        {
            std::string SegmentNbr_str, EnergyLoss_str;

            int SegmentNbr;
            float EnergyLoss;

            string firstline;
            string line;

            getline(inFile, firstline);

            while (getline(inFile, line))
            {
                stringstream ss(line);

                getline(ss, SegmentNbr_str, '\t');
                getline(ss, EnergyLoss_str, '\t');

                SegmentNbr = stoi(SegmentNbr_str);
                EnergyLoss = stof(EnergyLoss_str);

                // push/append the data into the vectors
                ICSegmentsVector.push_back(SegmentNbr);
                EnergyLossVector.push_back(EnergyLoss);

                // cout << "SetmentNbr: " << SegmentNbr << " EnergyLoss: " << EnergyLoss << endl;
            }
        }
        inFile.close();

        ICSegment2DVec.push_back(ICSegmentsVector);
        EnergyLoss2DVec.push_back(EnergyLossVector);
    }

    for (int idx = 0; idx < ICSegment2DVec[0].size(); idx++)
    {
        cout << "SegmentNbr: " << ICSegment2DVec[0][idx] << " EnergyLoss: " << EnergyLoss2DVec[0][idx] << endl;
    }

    // // ---------------------------------------------------------------------

    // TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    // TTree *tree = (TTree *)file->Get("tree_new");

    // // TFile *graphFile = new TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/graphCuts/FE9_X_PID_cut.root");
    // // TCutG *cut = (TCutG *)graphFile->Get("myCut");

    // Double_t F3_T_0;  // Time at F3 Diamond (array).
    // Double_t PID_T_0; // Time of Flight between F3 and FE9.
    // Double_t FE9_X_0; // Horizontal position at FE9 focal plane.
    // Double_t FE9_Y_0; // Vertical position at FE9 focal plane.
    // Double_t FE9_A_0; // Horizontal angle at FE9 focal plane.
    // Double_t FE9_B_0; // Vertical angle at FE9 focal plane.

    // Double_t Beam_T_0; // Time of Flight between FE9 and FE12.
    // Double_t S0_X_0;   // Horizontal position at FE12/S0 focal plane.
    // Double_t S0_Y_0;   // Vertical position at FE12/S0 focal plane.
    // Double_t S0_A_0;   // Horizontal angle at FE12/S0 focal plane.
    // Double_t S0_B_0;   // Vertical angle at FE12/S0 focal plane.

    // Double_t S1PID_T_0; // Time at S1 (average between anodes).
    // Double_t S1_X_0;    // Horizontal position at S1 focal plane.
    // Double_t S1_Y_0;    // Vertical position at S1 focal plane.
    // Double_t S1_A_0;    // Horizontal angle at S1 focal plane.
    // Double_t S1_B_0;    // Vertical angle at S1 focal plane.

    // tree->SetBranchAddress("F3_T_0", &F3_T_0);
    // tree->SetBranchAddress("PID_T_0", &PID_T_0);
    // tree->SetBranchAddress("FE9_X_0", &FE9_X_0);
    // tree->SetBranchAddress("FE9_Y_0", &FE9_Y_0);
    // tree->SetBranchAddress("FE9_A_0", &FE9_A_0);
    // tree->SetBranchAddress("FE9_B_0", &FE9_B_0);

    // tree->SetBranchAddress("Beam_T_0", &Beam_T_0);
    // tree->SetBranchAddress("S0_X_0", &S0_X_0);
    // tree->SetBranchAddress("S0_Y_0", &S0_Y_0);
    // tree->SetBranchAddress("S0_A_0", &S0_A_0);
    // tree->SetBranchAddress("S0_B_0", &S0_B_0);

    // tree->SetBranchAddress("S1PID_T_0", &S1PID_T_0);
    // tree->SetBranchAddress("S1_X_0", &S1_X_0);
    // tree->SetBranchAddress("S1_Y_0", &S1_Y_0);
    // tree->SetBranchAddress("S1_A_0", &S1_A_0);
    // tree->SetBranchAddress("S1_B_0", &S1_B_0);

    // TCanvas *c1 = new TCanvas("c1", "IC Canvas", 2500, 750);
    // c1->Divide(2, 1);
    // c1->cd(1);

    // TList *hlist_E_Cal = new TList();

    // TH2F *h_IC_E_Cal = new TH2F("IC_E_Cal", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 50); // IC Energy Data
    // hlist_E_Cal->Add(h_IC_E_Cal);

    // vector<vector<double>> IC_E_Cal_vec;
    // Double_t IC_E_Cal[30];

    // for (int k = 0; k < 30; k++)
    // {
    //     string IC_E_SetBranch = "IC_E_Cal_" + to_string(k);
    //     tree->SetBranchAddress(IC_E_SetBranch.c_str(), &IC_E_Cal[k]); // .c_str() converts any string to CONSTANT
    // }

    // int entries = tree->GetEntries();
    // cout << "Number of Entries: " << entries << endl;

    // int count = 0;

    // tree->GetEntry(363409);

    // for (int ix = 0; ix < 30; ix++)
    // {
    //     cout << "IC Segment: " << ix << " >> EnergyLoss: " << IC_E_Cal[ix] << endl;

    //     h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);
    // }

    // count = count + 1;

    // cout << count << endl;

    // // gStyle->SetPalette(kRed);
    // // h_IC_E_Cal->SetFillColor(kRed);
    // // h_IC_E_Cal->SetLineColor(kBlue);
    // // h_IC_E_Cal->Draw("colz");
    // h_IC_E_Cal->Draw("BOX");

    // // Predefined Colours in ROOT++
    // // kRed, kBlue, kGreen, kYellow, kOrange, kMagenta, kCyan, kBlack, etc.

    // // Use .data() to pass raw pointers (like const Double_t*) to TGraph
    // // BeamEnergy 13MeV/u
    // TGraph *graph13_0 = new TGraph(n, Segment.data(), EnergyLoss13_0MeV.data());
    // graph13_0->SetMarkerStyle(20);
    // graph13_0->SetMarkerColor(kBlack);
    // graph13_0->Draw("P SAME");

    // // BeamEnergy 15MeV/u
    // TGraph *graph15_0 = new TGraph(n, Segment.data(), EnergyLoss15_0MeV.data());
    // graph15_0->SetMarkerStyle(20);
    // graph15_0->SetMarkerColor(kBlue);
    // graph15_0->Draw("P SAME");

    // // BeamEnergy 15.5MeV/u
    // TGraph *graph15_5 = new TGraph(n, Segment.data(), EnergyLoss15_5MeV.data());
    // graph15_5->SetMarkerStyle(20);
    // graph15_5->SetMarkerColor(kRed);
    // graph15_5->Draw("P SAME");

    // // BeamEnergy 16MeV/u
    // TGraph *graph16_0 = new TGraph(n, Segment.data(), EnergyLoss16_0MeV.data());
    // graph16_0->SetMarkerStyle(20);
    // graph16_0->SetMarkerColor(kOrange);
    // graph16_0->Draw("P SAME");

    // // BeamEnergy 17.5MeV/u
    // TGraph *graph17_5 = new TGraph(n, Segment.data(), EnergyLoss17_5MeV.data());
    // graph17_5->SetMarkerStyle(20);
    // graph17_5->SetMarkerColor(kMagenta);
    // graph17_5->Draw("P SAME");

    // // Define ChiSqr Map
    // std::map<std::string, double> ChiSqr;

    // for (int k = 0; k < LISEFiles.size(); k++)
    // {
    //     if (LISEFiles[k] == "13.0MeV.txt")
    //     {
    //         double chisqr = 0;
    //         for (int i = 0; i < n; i++)
    //         {
    //             if (IC_E_Cal[i] >= 10)
    //             {
    //                 chisqr += pow((EnergyLoss13_0MeV.at(i) - IC_E_Cal[i]), 2);
    //             }
    //         }
    //         ChiSqr["13.0"] = chisqr;
    //         cout << "Chi Square = " << ChiSqr["13.0"] << endl;
    //     }
    //     else if (LISEFiles[k] == "15.0MeV.txt")
    //     {
    //         double chisqr = 0;
    //         for (int i = 0; i < n; i++)
    //         {
    //             if (IC_E_Cal[i] >= 10)
    //             {
    //                 chisqr += pow((EnergyLoss15_0MeV.at(i) - IC_E_Cal[i]), 2);
    //             }
    //         }
    //         ChiSqr["15.0"] = chisqr;
    //         cout << "Chi Square = " << ChiSqr["15.0"] << endl;
    //     }
    //     else if (LISEFiles[k] == "15.5MeV.txt")
    //     {
    //         double chisqr = 0;
    //         for (int i = 0; i < n; i++)
    //         {
    //             if (IC_E_Cal[i] >= 10)
    //             {
    //                 chisqr += pow((EnergyLoss15_5MeV.at(i) - IC_E_Cal[i]), 2);
    //             }
    //         }
    //         ChiSqr["15.5"] = chisqr;
    //         cout << "Chi Square = " << ChiSqr["15.5"] << endl;
    //     }
    //     else if (LISEFiles[k] == "16.0MeV.txt")
    //     {
    //         double chisqr = 0;
    //         for (int i = 0; i < n; i++)
    //         {
    //             if (IC_E_Cal[i] >= 10)
    //             {
    //                 chisqr += pow((EnergyLoss16_0MeV.at(i) - IC_E_Cal[i]), 2);
    //             }
    //         }
    //         ChiSqr["16.0"] = chisqr;
    //         cout << "Chi Square = " << ChiSqr["16.0"] << endl;
    //     }
    //     else if (LISEFiles[k] == "17.5MeV.txt")
    //     {
    //         double chisqr = 0;
    //         for (int i = 0; i < n; i++)
    //         {
    //             if (IC_E_Cal[i] >= 10)
    //             {
    //                 chisqr += pow((EnergyLoss17_5MeV.at(i) - IC_E_Cal[i]), 2);
    //             }
    //         }
    //         ChiSqr["17.5"] = chisqr;
    //         cout << "Chi Square = " << ChiSqr["17.5"] << endl;
    //     }
    // }

    // vector<double> energyArr, chiSqrArr;

    // for (auto iterator = ChiSqr.begin(); iterator != ChiSqr.end(); iterator++)
    // {
    //     cout << "Key: " << iterator->first << " Value: " << iterator->second << endl;
    //     // double energy = stof(iterator->first);
    //     // double chisq = iterator->second;

    //     energyArr.push_back(stof(iterator->first));
    //     chiSqrArr.push_back(iterator->second);
    // }

    // c1->cd(2);

    // TGraph *graphChiSqrvsE = new TGraph(energyArr.size(), energyArr.data(), chiSqrArr.data());

    // graphChiSqrvsE->SetTitle("ChiSqr vs Energy; Energy in MeV; ChiSqr");
    // graphChiSqrvsE->SetLineColor(kRed);
    // graphChiSqrvsE->SetMarkerStyle(20);
    // graphChiSqrvsE->SetLineWidth(2);

    // // graphChiSqrvsE->Draw("AP");
    // graphChiSqrvsE->Draw("AL");

    // c1->Update();
}