#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>

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

struct BeamEnergyData {
    std::vector<int> ICSegment;
    std::vector<double> EnergyLoss;
};

double calculateSlope(double X1, double Y1, double X2, double Y2) {
    return (Y2 - Y1) / (X2 - X1);
}

double calculateDy(double Y1, double Y2) {
    return Y2 - Y1;
}

vector<double> Segment = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                          17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
const int n = 30;

void loadLISEData(const vector<double> &BeamEnergies, map<string, BeamEnergyData> &LISEDataMap) {
    for (auto &energy : BeamEnergies) {
        // string filename = "/u/ddas/software/work/artemis-oedo/output/Analysis/" + to_string(energy) + "MeV.txt";
        string liseFile = Form("%0.1fMeV.txt", energy);
        string filename = "/u/ddas/software/work/artemis-oedo/output/Analysis/" + liseFile;
        ifstream inFile(filename.c_str());
        
        if (!inFile) {
            cerr << "Unable to open file: " << filename << endl;
            continue;
        }
        
        BeamEnergyData data;
        string line, SegmentNbr_str, EnergyLoss_str;
        getline(inFile, line);  // Skip header
        
        while (getline(inFile, line)) {
            stringstream ss(line);
            getline(ss, SegmentNbr_str, '\t');
            getline(ss, EnergyLoss_str, '\t');
            
            data.ICSegment.push_back(stoi(SegmentNbr_str));
            data.EnergyLoss.push_back(stod(EnergyLoss_str));
        }
        
        LISEDataMap[to_string(energy) + "MeV"] = data;
        inFile.close();
    }
}

// void calculateChiSquare(const map<string, BeamEnergyData> &LISEDataMap, const vector<double> &IC_E_Cal) {
//     for (const auto &entry : LISEDataMap) {
//         const auto &energyLoss = entry.second.EnergyLoss;
//         double chiSquare = 0.0;
        
//         for (int i = 0; i < n - 1; ++i) {
//             double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
//             double slopeLISESim = calculateSlope(Segment[i], energyLoss[i], Segment[i + 1], energyLoss[i + 1]);
            
//             if (slopeLISESim > 0 || dYExptData > -5.0) {
//                 chiSquare += pow((energyLoss[i] - IC_E_Cal[i]), 2);
//             }
//         }
        
//         cout << "Chi Square for " << entry.first << " = " << chiSquare << endl;
//     }
// }

// void calculateChiSquare(const map<string, BeamEnergyData> &LISEDataMap, const vector<double> &IC_E_Cal) {
//     for (const auto &entry : LISEDataMap) {
//         const auto &energyLoss = entry.second.EnergyLoss;
//         double chiSquare = 0.0;  // Reset chi-square for each beam energy
        
//         for (int i = 0; i < n - 1; ++i) {
//             double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
//             double slopeLISESim = calculateSlope(Segment[i], energyLoss[i], Segment[i + 1], energyLoss[i + 1]);
            
//             if (slopeLISESim > 0 || dYExptData > -5.0) {
//                 chiSquare += pow((energyLoss[i] - IC_E_Cal[i]), 2);
//             }
//         }
        
//         // Print chi-square for this particular beam energy
//         cout << "Chi Square for " << entry.first << " = " << chiSquare << endl;
//     }
// }

void calculateChiSquare(const map<string, BeamEnergyData> &LISEDataMap, const vector<double> &IC_E_Cal, map<string, double> &chiSquareResults) {
    for (const auto &entry : LISEDataMap) {
        const auto &energyLoss = entry.second.EnergyLoss;
        double chiSquare = 0.0;  // Reset chi-square for each beam energy
        
        for (int i = 0; i < n - 1; ++i) {
            double dYExptData = calculateDy(IC_E_Cal[i], IC_E_Cal[i + 1]);
            double slopeLISESim = calculateSlope(Segment[i], energyLoss[i], Segment[i + 1], energyLoss[i + 1]);
            
            if (slopeLISESim > 0 || dYExptData > -5.0) {
                chiSquare += pow((energyLoss[i] - IC_E_Cal[i]), 2);
            }
        }
        
        // Store chi-square in the map
        chiSquareResults[entry.first] = chiSquare;
        cout << "Chi Square for " << entry.first << " = " << chiSquare << endl;
    }
}


void plotData(TTree *tree, const vector<double> &BeamEnergies, const map<string, BeamEnergyData> &LISEDataMap) {
    TCanvas *canvas = new TCanvas("canvas", "IC Canvas", 2500, 750);
    canvas->Divide(2, 1);
    canvas->cd(1);
    
    TH2F *h_IC_E_Cal = new TH2F("IC_E_Cal", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 50);
    h_IC_E_Cal->Draw("BOX");
    
    Double_t IC_E_Cal[30];
    for (int k = 0; k < 30; k++) {
        tree->SetBranchAddress(("IC_E_Cal_" + to_string(k)).c_str(), &IC_E_Cal[k]);
    }
    
    int entries = tree->GetEntries();
    tree->GetEntry(363409);  // Sample entry for demonstration
    
    for (int ix = 0; ix < 30; ix++) {
        h_IC_E_Cal->Fill(ix, IC_E_Cal[ix]);
    }
    
    for (size_t i = 0; i < BeamEnergies.size(); ++i) {
        string key = to_string(BeamEnergies[i]) + "MeV";
        TGraph *graph = new TGraph(n, Segment.data(), LISEDataMap.at(key).EnergyLoss.data());
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(i + 2); // Unique color for each graph
        graph->Draw("P SAME");
    }
}

void fit_one_v3() {
    vector<double> BeamEnergies = {12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5};
    map<string, BeamEnergyData> LISEDataMap;
    
    loadLISEData(BeamEnergies, LISEDataMap);
    
    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    TTree *tree = (TTree *)file->Get("tree_new");
    
    vector<double> IC_E_Cal(30, 0.0); // Adjust as needed for actual IC_E_Cal values
    
    plotData(tree, BeamEnergies, LISEDataMap);

    map<string, double> chiSquareResults;
    calculateChiSquare(LISEDataMap, IC_E_Cal, chiSquareResults);

    // calculateChiSquare(LISEDataMap, IC_E_Cal);
}
