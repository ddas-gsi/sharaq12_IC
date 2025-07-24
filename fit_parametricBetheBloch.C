#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

double parametric_bethe_bloch(double x, double Z, double a, double b, double c)
{
    double beta = 0.6; // Place holder value, use experimental beta later
    double gamma = 1 / (std::sqrt(1 - pow(beta, 2)));

    // Parametric Bragg Curve
    double stopping_power = (a * (Z * *2 / beta * *2) * exp(-b * x) * (1 + c * log(x + 1)));

    return stopping_power;
}

double calculateSlope(double X1, double Y1, double X2, double Y2)
{
    double slope = (Y2 - Y1) / (X2 - X1);
    return slope;
}

double calculateDy(double Y1, double Y2)
{
    double dY = Y2 - Y1;
    return dY;
}

template <typename K, typename V>
K getEnergyValue(const map<K, V> &ChiSqrMap)
{
    // Create a vector of pairs to sort by value
    vector<pair<K, V>> ChiSqrMapVector(ChiSqrMap.begin(), ChiSqrMap.end());

    // Sort the vector based on the second element of the pair (i.e., the value)
    sort(ChiSqrMapVector.begin(), ChiSqrMapVector.end(), [](const pair<K, V> &a, const pair<K, V> &b)
         { return a.second < b.second; });
    // Return the key with the lowest value
    return ChiSqrMapVector.front().first;
}

int getSegmentIDOfBraggCurve(vector<double> &EnergyLossExVector, double Val)
{
    for (int i = 0; i < EnergyLossExVector.size(); i++)
    {
        if (EnergyLossExVector[i] == Val)
        {
            return i;
        }
    }
    return -1;
}

vector<double> selectCorrectEnergyLossVec(vector<double> &EnergyLoss12_0MeV, vector<double> &EnergyLoss12_5MeV,
                                          vector<double> &EnergyLoss13_0MeV, vector<double> &EnergyLoss13_5MeV,
                                          vector<double> &EnergyLoss14_0MeV, vector<double> &EnergyLoss14_5MeV,
                                          vector<double> &EnergyLoss15_0MeV, vector<double> &EnergyLoss15_5MeV,
                                          vector<double> &EnergyLoss16_0MeV, vector<double> &EnergyLoss16_5MeV,
                                          vector<double> &EnergyLoss17_0MeV, vector<double> &EnergyLoss17_5MeV, string energyValKey)
{
    if (energyValKey == "12.0")
    {
        return EnergyLoss12_0MeV;
    }
    else if (energyValKey == "12.5")
    {
        return EnergyLoss12_5MeV;
    }
    else if (energyValKey == "13.0")
    {
        return EnergyLoss13_0MeV;
    }
    else if (energyValKey == "13.5")
    {
        return EnergyLoss13_5MeV;
    }
    else if (energyValKey == "14.0")
    {
        return EnergyLoss14_0MeV;
    }
    else if (energyValKey == "14.5")
    {
        return EnergyLoss14_5MeV;
    }
    else if (energyValKey == "15.0")
    {
        return EnergyLoss15_0MeV;
    }
    else if (energyValKey == "15.5")
    {
        return EnergyLoss15_5MeV;
    }
    else if (energyValKey == "16.0")
    {
        return EnergyLoss16_0MeV;
    }
    else if (energyValKey == "16.5")
    {
        return EnergyLoss16_5MeV;
    }
    else if (energyValKey == "17.0")
    {
        return EnergyLoss17_0MeV;
    }
    else if (energyValKey == "17.5")
    {
        return EnergyLoss17_5MeV;
    }
}

// // Energy Loss for BeamEnergy 15.5 MeV
// vector<double> EnergyLoss15_5MeV = {19.307, 19.777, 20.276, 20.82, 21.413, 22.06, 22.77,
//                                     23.554, 24.425, 25.4, 26.502, 27.759, 29.211, 30.905,
//                                     32.892, 35.247, 37.76, 35.479, 12.255, 0, 0, 0, 0, 0,
//                                     0, 0, 0, 0, 0, 0};

vector<double> Segment = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                          17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

const int n = 30;

void fit_one()
{

    // ---------------------------------------------------------------------
    // Load LISE++ Simulation Data
    // ---------------------------------------------------------------------
    vector<double> BeamEnergies = {12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5};

    vector<std::string> LISEFiles;

    for (int i = 0; i < BeamEnergies.size(); i++)
    {
        std::string liseFile = Form("%0.1fMeV.txt", BeamEnergies[i]); // + "MeV.txt";
        // cout << liseFile << endl;

        LISEFiles.push_back(liseFile.c_str());
    }

    std::vector<int> ICSegment12_0MeV, ICSegment12_5MeV, ICSegment13_0MeV, ICSegment13_5MeV,
        ICSegment14_0MeV, ICSegment14_5MeV, ICSegment15_0MeV, ICSegment15_5MeV,
        ICSegment16_0MeV, ICSegment16_5MeV, ICSegment17_0MeV, ICSegment17_5MeV;

    std::vector<double> EnergyLoss12_0MeV, EnergyLoss12_5MeV, EnergyLoss13_0MeV, EnergyLoss13_5MeV,
        EnergyLoss14_0MeV, EnergyLoss14_5MeV, EnergyLoss15_0MeV, EnergyLoss15_5MeV,
        EnergyLoss16_0MeV, EnergyLoss16_5MeV, EnergyLoss17_0MeV, EnergyLoss17_5MeV;

    for (int k = 0; k < LISEFiles.size(); k++)
    {
        std::string liseFile = "/u/ddas/software/work/artemis-oedo/output/Analysis/" + LISEFiles[k];
        ifstream inFile(liseFile.c_str());

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
                if (LISEFiles[k] == "12.0MeV.txt")
                {
                    ICSegment12_0MeV.push_back(SegmentNbr);
                    EnergyLoss12_0MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "12.5MeV.txt")
                {
                    ICSegment12_5MeV.push_back(SegmentNbr);
                    EnergyLoss12_5MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "13.0MeV.txt")
                {
                    ICSegment13_0MeV.push_back(SegmentNbr);
                    EnergyLoss13_0MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "13.5MeV.txt")
                {
                    ICSegment13_5MeV.push_back(SegmentNbr);
                    EnergyLoss13_5MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "14.0MeV.txt")
                {
                    ICSegment14_0MeV.push_back(SegmentNbr);
                    EnergyLoss14_0MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "14.5MeV.txt")
                {
                    ICSegment14_5MeV.push_back(SegmentNbr);
                    EnergyLoss14_5MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "15.0MeV.txt")
                {
                    ICSegment15_0MeV.push_back(SegmentNbr);
                    EnergyLoss15_0MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "15.5MeV.txt")
                {
                    ICSegment15_5MeV.push_back(SegmentNbr);
                    EnergyLoss15_5MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "16.0MeV.txt")
                {
                    ICSegment16_0MeV.push_back(SegmentNbr);
                    EnergyLoss16_0MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "16.5MeV.txt")
                {
                    ICSegment16_5MeV.push_back(SegmentNbr);
                    EnergyLoss16_5MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "17.0MeV.txt")
                {
                    ICSegment17_0MeV.push_back(SegmentNbr);
                    EnergyLoss17_0MeV.push_back(EnergyLoss);
                }
                else if (LISEFiles[k] == "17.5MeV.txt")
                {
                    ICSegment17_5MeV.push_back(SegmentNbr);
                    EnergyLoss17_5MeV.push_back(EnergyLoss);
                }

                // cout << "SetmentNbr: " << SegmentNbr << " EnergyLoss: " << EnergyLoss << endl;
            }
        }
        inFile.close();
    }

    // for (int idx = 0; idx < ICSegment15_5MeV.size(); idx++)
    // {
    //     cout << "SegmentNbr: " << ICSegment17_0MeV[idx] << " EnergyLoss: " << EnergyLoss17_0MeV[idx] << endl;
    // }

    // ---------------------------------------------------------------------

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

    TCanvas *c1 = new TCanvas("c1", "IC Canvas", 2500, 750);
    c1->Divide(2, 1);
    c1->cd(1);

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

    tree->GetEntry(363409); // 15.5 MeV
    // tree->GetEntry(358);    // 14.0 MeV
    // tree->GetEntry(152); // 14.5 MeV

    for (int ix = 0; ix < 30; ix++)
    {
        cout << "IC Segment: " << ix << " >> EnergyLoss: " << IC_E_Cal[ix] << endl;

        h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);
        // h_IC_E_Cal->Fill(ix, IC_E_Cal[ix] - 6);
    }

    count = count + 1;

    cout << count << endl;

    // gStyle->SetPalette(kRed);
    // h_IC_E_Cal->SetFillColor(kRed);
    // h_IC_E_Cal->SetLineColor(kBlue);
    // h_IC_E_Cal->Draw("colz");
    h_IC_E_Cal->Draw("BOX");

    // Predefined Colours in ROOT++
    // kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan, kTeal, kGreen, kSpring, kYellow, kOrange, kBlack, kGray, kWhite etc.

    // Use .data() to pass raw pointers (like const Double_t*) to TGraph
    // BeamEnergy 12.0MeV/u
    TGraph *graph12_0 = new TGraph(n, Segment.data(), EnergyLoss12_0MeV.data());
    graph12_0->SetMarkerStyle(20);
    graph12_0->SetMarkerColor(kPink + 10);
    graph12_0->Draw("P SAME");

    // BeamEnergy 12.5MeV/u
    TGraph *graph12_5 = new TGraph(n, Segment.data(), EnergyLoss12_5MeV.data());
    graph12_5->SetMarkerStyle(20);
    graph12_5->SetMarkerColor(kYellow + 1);
    graph12_5->Draw("P SAME");

    // BeamEnergy 13.0MeV/u
    TGraph *graph13_0 = new TGraph(n, Segment.data(), EnergyLoss13_0MeV.data());
    graph13_0->SetMarkerStyle(20);
    graph13_0->SetMarkerColor(kBlack);
    graph13_0->Draw("P SAME");

    // BeamEnergy 13.5MeV/u
    TGraph *graph13_5 = new TGraph(n, Segment.data(), EnergyLoss13_5MeV.data());
    graph13_5->SetMarkerStyle(20);
    graph13_5->SetMarkerColor(kViolet - 3);
    graph13_5->Draw("P SAME");

    // BeamEnergy 14.0MeV/u
    TGraph *graph14_0 = new TGraph(n, Segment.data(), EnergyLoss14_0MeV.data());
    graph14_0->SetMarkerStyle(20);
    graph14_0->SetMarkerColor(kRed - 4);
    graph14_0->Draw("P SAME");

    // BeamEnergy 14.5MeV/u
    TGraph *graph14_5 = new TGraph(n, Segment.data(), EnergyLoss14_5MeV.data());
    graph14_5->SetMarkerStyle(20);
    graph14_5->SetMarkerColor(kCyan + 2);
    graph14_5->Draw("P SAME");

    // BeamEnergy 15.0MeV/u
    TGraph *graph15_0 = new TGraph(n, Segment.data(), EnergyLoss15_0MeV.data());
    graph15_0->SetMarkerStyle(20);
    graph15_0->SetMarkerColor(kBlue);
    graph15_0->Draw("P SAME");

    // BeamEnergy 15.5MeV/u
    TGraph *graph15_5 = new TGraph(n, Segment.data(), EnergyLoss15_5MeV.data());
    graph15_5->SetMarkerStyle(20);
    graph15_5->SetMarkerColor(kRed);
    graph15_5->Draw("P SAME");

    // BeamEnergy 16.0MeV/u
    TGraph *graph16_0 = new TGraph(n, Segment.data(), EnergyLoss16_0MeV.data());
    graph16_0->SetMarkerStyle(20);
    graph16_0->SetMarkerColor(kOrange - 3);
    graph16_0->Draw("P SAME");

    // BeamEnergy 16.5MeV/u
    TGraph *graph16_5 = new TGraph(n, Segment.data(), EnergyLoss16_5MeV.data());
    graph16_5->SetMarkerStyle(20);
    graph16_5->SetMarkerColor(kSpring + 9);
    graph16_5->Draw("P SAME");

    // BeamEnergy 17.0MeV/u
    TGraph *graph17_0 = new TGraph(n, Segment.data(), EnergyLoss17_0MeV.data());
    graph17_0->SetMarkerStyle(20);
    graph17_0->SetMarkerColor(kAzure - 3);
    graph17_0->Draw("P SAME");

    // BeamEnergy 17.5MeV/u
    TGraph *graph17_5 = new TGraph(n, Segment.data(), EnergyLoss17_5MeV.data());
    graph17_5->SetMarkerStyle(20);
    graph17_5->SetMarkerColor(kMagenta);
    graph17_5->Draw("P SAME");

    // Define ChiSqrMap
    std::map<std::string, double> ChiSqrMap;

    for (int k = 0; k < LISEFiles.size(); k++)
    {
        if (LISEFiles[k] == "12.0MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                // double slopeExptData = calculateSlope(Segment[i], IC_E_Cal[i], Segment[i + 1], IC_E_Cal[i + 1]);
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss12_0MeV[i], Segment[i + 1], EnergyLoss12_0MeV[i + 1]);

                // cout << "Expt Slope: " << slopeExptData << " Expt dY: " << dYExptData << " LISESim Slope: " << slopeLISESim << endl;
                // cout << "Expt dY: " << dYExptData << " LISESim Slope: " << slopeLISESim << endl;

                // if (IC_E_Cal[i] >= 10)
                if (slopeLISESim > 0 || dYExptData > -5.0)
                {
                    chisqr += pow((EnergyLoss12_0MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["12.0"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["12.0"] << endl;
        }
        else if (LISEFiles[k] == "12.5MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss12_5MeV[i], Segment[i + 1], EnergyLoss12_5MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss12_5MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["12.5"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["12.5"] << endl;
        }
        else if (LISEFiles[k] == "13.0MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss13_0MeV[i], Segment[i + 1], EnergyLoss13_0MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss13_0MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["13.0"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["13.0"] << endl;
        }
        else if (LISEFiles[k] == "13.5MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss13_5MeV[i], Segment[i + 1], EnergyLoss13_5MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss13_5MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["13.5"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["13.5"] << endl;
        }
        else if (LISEFiles[k] == "14.0MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss14_0MeV[i], Segment[i + 1], EnergyLoss14_0MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss14_0MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["14.0"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["14.0"] << endl;
        }
        else if (LISEFiles[k] == "14.5MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss14_5MeV[i], Segment[i + 1], EnergyLoss14_5MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss14_5MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["14.5"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["14.5"] << endl;
        }
        else if (LISEFiles[k] == "15.0MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss15_0MeV[i], Segment[i + 1], EnergyLoss15_0MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss15_0MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["15.0"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["15.0"] << endl;
        }
        else if (LISEFiles[k] == "15.5MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss15_5MeV[i], Segment[i + 1], EnergyLoss15_5MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss15_5MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["15.5"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["15.5"] << endl;
        }
        else if (LISEFiles[k] == "16.0MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss16_0MeV[i], Segment[i + 1], EnergyLoss16_0MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss16_0MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["16.0"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["16.0"] << endl;
        }
        else if (LISEFiles[k] == "16.5MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss16_5MeV[i], Segment[i + 1], EnergyLoss16_5MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss16_5MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["16.5"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["16.5"] << endl;
        }
        else if (LISEFiles[k] == "17.0MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss17_0MeV[i], Segment[i + 1], EnergyLoss17_0MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss17_0MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["17.0"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["17.0"] << endl;
        }
        else if (LISEFiles[k] == "17.5MeV.txt")
        {
            double chisqr = 0;
            for (int i = 0; i < n; i++)
            {
                double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
                double slopeLISESim = calculateSlope(Segment[i], EnergyLoss17_5MeV[i], Segment[i + 1], EnergyLoss17_5MeV[i + 1]);
                if (slopeLISESim > 0 || dYExptData > -5.0)
                // if (IC_E_Cal[i] >= 10)
                {
                    chisqr += pow((EnergyLoss17_5MeV.at(i) - IC_E_Cal[i]), 2);
                }
            }
            ChiSqrMap["17.5"] = chisqr;
            cout << "Chi Square = " << ChiSqrMap["17.5"] << endl;
        }
    }

    vector<double> energyArr, chiSqrArr;
    for (auto iterator = ChiSqrMap.begin(); iterator != ChiSqrMap.end(); iterator++)
    {
        cout << "Key: " << iterator->first << " Value: " << iterator->second << endl;
        // double energy = stof(iterator->first);
        // double chisq = iterator->second;

        energyArr.push_back(stof(iterator->first));
        chiSqrArr.push_back(iterator->second);
    }

    c1->cd(2);

    TGraph *graphChiSqrvsE = new TGraph(energyArr.size(), energyArr.data(), chiSqrArr.data());

    graphChiSqrvsE->SetTitle("ChiSqr vs Energy; Energy in MeV; ChiSqr");
    graphChiSqrvsE->SetLineColor(kRed);
    graphChiSqrvsE->SetMarkerStyle(20);
    graphChiSqrvsE->SetLineWidth(2);

    // graphChiSqrvsE->Draw("AP");
    graphChiSqrvsE->Draw("AL");

    c1->Update();

    // cout << EnergyLoss12_0MeV.size() << endl;

    // Print out the EnergyValue with Lowest chiSquareValue
    std::string energyValKey = getEnergyValue(ChiSqrMap);
    cout << "Energy of Bragg Curve with Lowest ChiSqr Value: " << energyValKey << " and Value of ChiSqr: " << ChiSqrMap[energyValKey] << endl;

    // Select the correct Bragg Curve Vector / Energy Loss Vector
    vector<double> correctEnergyLossVec = selectCorrectEnergyLossVec(EnergyLoss12_0MeV, EnergyLoss12_5MeV, EnergyLoss13_0MeV, EnergyLoss13_5MeV,
                                                                     EnergyLoss14_0MeV, EnergyLoss14_5MeV, EnergyLoss15_0MeV, EnergyLoss15_5MeV,
                                                                     EnergyLoss16_0MeV, EnergyLoss16_5MeV, EnergyLoss17_0MeV, EnergyLoss17_5MeV, energyValKey);

    cout << "Ion Stopped at SegmentID: " << getSegmentIDOfBraggCurve(EnergyLoss15_5MeV, 0) << endl;
}
