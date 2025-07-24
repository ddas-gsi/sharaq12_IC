#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <ctime>   // For timestamp
#include <sstream> // For string formatting

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"

// Energy Loss for BeamEnergy 15 MeV
vector<double> EnergyLoss15MeV = {20.142, 20.676, 21.257, 21.889, 22.582, 23.346, 24.192,
                                  25.138, 26.203, 27.416, 28.812, 30.437, 32.351, 34.599,
                                  37.211, 37.504, 19.631, 0.494, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0};

// Energy Loss for BeamEnergy 15.5 MeV
vector<double> EnergyLoss15_5MeV = {19.307, 19.777, 20.276, 20.82, 21.413, 22.06, 22.77,
                                    23.554, 24.425, 25.4, 26.502, 27.759, 29.211, 30.905,
                                    32.892, 35.247, 37.76, 35.479, 12.255, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0};

vector<double> Segment = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                          17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

const int n = 30;

std::string beamEnergyStr = "15.5 MeV";

void sharaq12_BraggFit()
{
    TFile *file = TFile::Open("/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/902_icHistos.root", "READ");
    if (!file || !file->IsOpen())
    {
        std::cerr << "Error: Couldn't open the input file." << std::endl;
        return;
    }

    // Open the output ROOT file in "UPDATE" mode so that it keeps all previous canvases
    TFile *outputFile = TFile::Open("/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/902_icBraggPeak.root", "UPDATE");
    if (!outputFile || !outputFile->IsOpen())
    {
        std::cerr << "Error: Couldn't create or open the output file." << std::endl;
        file->Close();
        return;
    }

    // Create title of the Bragg Peak with Beam Energy value
    std::string hist2DTitle = "Bragg Peak for Beam Energy " + beamEnergyStr;

    // Create a canvas with a unique name
    std::ostringstream oss;
    std::time_t t = std::time(nullptr);
    std::tm *tm = std::localtime(&t);
    oss << "c1_" << (tm->tm_year + 1900) << "_" // Year
        << (tm->tm_mon + 1) << "_"              // Month
        << tm->tm_mday << "_"                   // Day
        << tm->tm_hour << "_"                   // Hour
        << tm->tm_min << "_"                    // Minute
        << tm->tm_sec;                          // Second
    std::string canvasName = oss.str();

    TCanvas *c1 = new TCanvas(canvasName.c_str(), "Bragg Peak", 800, 600);

    TH2F *hist2D = (TH2F *)file->Get("IC_E_Cal_2Dcut");

    if (hist2D)
    {

        hist2D->SetTitle(hist2DTitle.c_str()); // Set the title of the histogram

        hist2D->Draw("COLZ");

        // Use .data() to pass raw pointers (like const Double_t*) to TGraph
        // TGraph *graph = new TGraph(n, Segment.data(), EnergyLoss15MeV.data());
        TGraph *graph = new TGraph(n, Segment.data(), EnergyLoss15_5MeV.data());
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kRed);
        graph->Draw("P SAME");

        c1->Update();

        // Write the canvas with a unique name to the output file
        outputFile->cd(); // Switch to the output file directory
        c1->Write();      // Save the canvas
        std::cout << "Canvas saved to output file with name: " << canvasName << std::endl;
    }
    else
    {
        std::cerr << "Error: Couldn't retrieve the histogram 'IC_E_Cal_2Dcut' from the file." << std::endl;
    }

    // Clean up
    outputFile->Close();
    file->Close();
    delete c1;
}
