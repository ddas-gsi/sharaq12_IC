// Use this C++ code fastBeamChargePlot.C .. Not the fastBeamChargePlot.py file.

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGraph.h> // Include for TGraph
#include <iostream>

int fastBeamChargePlot()
{

    // Eloss value for beam Energy 78.11307 MeV/u and BRho = 3.2437 Tm
    // This array has 30 values, corresponding to pads 0-29.
    double Eloss[] = {
        4.429, 4.433, 4.437, 4.441, 4.445, 4.449, 4.453, 4.457,
        4.461, 4.465, 4.469, 4.473, 4.477, 4.481, 4.485, 4.489,
        4.493, 4.493, 4.480, 4.485, 4.490, 4.495, 4.499, 4.504,
        4.509, 4.514, 4.519, 4.523, 4.528, 4.533};
    const int nElossPoints = sizeof(Eloss) / sizeof(Eloss[0]); // Determine size dynamically

    double scalingFactors[] = {
        0.01884680851, 0.02141545894, 0.02035321101, 0.02009502262,
        0.01785140562, 0.01815918367, 0.02095529412, 0.02102358491,
        0.02074883721, 0.02157004831, 0.02013063063, 0.02412621359,
        0.01905106383, 0.02333854167, 0.02166666667, 0.02290306122,
        0.02455191257, 0.02269191919, 0.02395721925, 0.02231343284,
        0.02388297872, 0.02281725888, 0.02380423280, 0.02286294416,
        0.02437297297, 0.02351041667, 0.02378421053, 0.02368062827,
        0.02460869565, 0.02233004926};

    // Define pad IDs (should match the size of Eloss and scalingFactors)
    const int nPads = 30; // Assuming 30 pads, consistent with array sizes
    TString branchPrefix = "IC_C_";
    TString branchNames[nPads];
    Double_t IC_C_Values[nPads];

    for (int i = 0; i < nPads; ++i)
    {
        branchNames[i] = branchPrefix + std::to_string(i);
    }

    // Open ROOT file and get TTree
    TString inputFile = "/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/opticsOutput/FastBeam_0089_2024.root";
    TFile *file = TFile::Open(inputFile, "UPDATE");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return 1;
    }

    TTree *tree = (TTree *)file->Get("tree_new");
    if (!tree)
    {
        std::cerr << "TTree 'tree_new' not found in file: " << inputFile << std::endl;
        file->Close();
        return 1;
    }

    // Set branch addresses
    for (int i = 0; i < nPads; ++i)
    {
        tree->SetBranchAddress(branchNames[i], &IC_C_Values[i]);
    }

    Long64_t nEntries = tree->GetEntries();
    std::cout << "Number of entries: " << nEntries << std::endl;

    // Create canvas and histogram
    TCanvas *c1 = new TCanvas("c1", "IC Plots", 1500, 600);
    c1->Divide(2, 1);

    TList *hlist = new TList();

    TH2F *h2D_FAST_BEAM_IC_Charge = new TH2F("h2D_FAST_BEAM_IC_Charge",
                                             "FAST BEAM IC Charge vs IC Pad; Pad; IC Charge [au]",
                                             31, -0.5, 30.5, 1000, 0, 3000);
    hlist->Add(h2D_FAST_BEAM_IC_Charge);

    TH2F *h2D_FAST_BEAM_IC_Energy_Calibrated = new TH2F("h2D_FAST_BEAM_IC_Energy_Calibrated",
                                                        "FAST BEAM IC Energy vs IC Pad; Pad; IC Energy [MeV/u]",
                                                        31, -0.5, 30.5, 1000, 0, 30);
    hlist->Add(h2D_FAST_BEAM_IC_Energy_Calibrated);

    // Loop over entries
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        if (i % 100000 == 0)
        {
            std::cout << "Processing event: " << i << std::endl;
        }

        tree->GetEntry(i);

        for (int padID = 0; padID < nPads; ++padID)
        {
            double val = IC_C_Values[padID];
            double energy = val * scalingFactors[padID];
            if (val > 0)
            {
                h2D_FAST_BEAM_IC_Charge->Fill(padID, val);
                h2D_FAST_BEAM_IC_Energy_Calibrated->Fill(padID, energy);
            }
        }
    }

    // Create arrays for TGraph
    double padIDs[nPads];
    for (int i = 0; i < nPads; ++i)
    {
        padIDs[i] = i;
    }

    // Create TGraph for Eloss values
    TGraph *gEloss = new TGraph(nElossPoints, padIDs, Eloss);
    gEloss->SetName("gEloss");
    gEloss->SetTitle("LISE++ Eloss for 78.11307 MeV/u;Pad;Eloss [MeV/u]");
    gEloss->SetMarkerStyle(20); // Circle marker
    gEloss->SetMarkerSize(1.2);
    gEloss->SetMarkerColor(kRed); // Red color
    gEloss->SetLineColor(kRed);   // Red line
    gEloss->SetLineWidth(2);
    // Add the graph to the list of objects to be saved (optional, but good practice)
    hlist->Add(gEloss);

    // Draw histograms and graph
    c1->cd(1);
    h2D_FAST_BEAM_IC_Charge->SetStats(0);
    h2D_FAST_BEAM_IC_Charge->Draw("COLZ");

    c1->cd(2);
    h2D_FAST_BEAM_IC_Energy_Calibrated->SetStats(0);
    h2D_FAST_BEAM_IC_Energy_Calibrated->Draw("COLZ");
    // Draw the Eloss graph on top of the 2D histogram
    // "LP" draws points and a line. "SAME" draws it on the current pad without clearing.
    gEloss->Draw("LP SAME");
    c1->Update();

    // Save histogram to a new ROOT file
    TString outputFileName = "/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/opticsOutput/FastBeam_0089_2024_hist1107_withEloss.root";
    TFile *outFile = new TFile(outputFileName, "RECREATE");
    if (!outFile || outFile->IsZombie())
    {
        std::cerr << "Error creating output file: " << outputFileName << std::endl;
        file->Close();
        return 1;
    }

    c1->Write();                                                 // Write the canvas to the output file
    hlist->Write("HistogramAndGraphLists", TObject::kSingleKey); // Write the histogram and graph to the output file
    outFile->Close();                                            // Close the output file
    std::cout << "Histograms and graph saved to '" << outputFileName << "'" << std::endl;

    // Wait for user to close canvas
    std::cout << "Press enter to exit..." << std::endl;
    std::cin.get();

    file->Close();
    // Clean up
    delete c1;
    delete hlist;
    // TGraph and TH2F are owned by hlist and will be deleted when hlist is deleted.
    // However, it's good practice to delete dynamically allocated objects if not owned by a list.
    // In this case, since they are in hlist, deleting hlist will handle them.

    return 0;
}