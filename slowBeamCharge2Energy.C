#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGraph.h> // Include for TGraph
#include <iostream>

int slowBeamCharge2Energy(int RUN_Nbr)
{
    // Eloss value for beam Energy 78.11307 MeV/u and BRho = 3.2437 Tm
    // This array has 30 values, corresponding to pads 0-29.
    double Eloss[] = {
        4.429, 4.433, 4.437, 4.441, 4.445, 4.449, 4.453, 4.457,
        4.461, 4.465, 4.469, 4.473, 4.477, 4.481, 4.485, 4.489,
        4.493, 4.493, 4.480, 4.485, 4.490, 4.495, 4.499, 4.504,
        4.509, 4.514, 4.519, 4.523, 4.528, 4.533};
    const int nElossPoints = sizeof(Eloss) / sizeof(Eloss[0]); // Determine size dynamically

    // Scaling factors for each pad, corresponding to the Eloss values for beam Energy 78.11307 MeV/u and BRho = 3.2437 Tm
    double scalingFactors[] = {
        0.01884680851, 0.02141545894, 0.02035321101, 0.02009502262,
        0.01785140562, 0.01815918367, 0.02095529412, 0.02102358491,
        0.02074883721, 0.02157004831, 0.02013063063, 0.02412621359,
        0.01905106383, 0.02333854167, 0.02166666667, 0.02290306122,
        0.02455191257, 0.02269191919, 0.02395721925, 0.02231343284,
        0.02388297872, 0.02281725888, 0.02380423280, 0.02286294416,
        0.02437297297, 0.02351041667, 0.02378421053, 0.02368062827,
        0.02460869565, 0.02233004926};

    // Modifying the scaling factors to better readjust pad {6,7,10,12,13,14,15,16}
    double scalingFactorsMod[] = {
        1, 1, 1, 1,
        1, 1, 1.0312246, 1.080747,
        1, 1, 1.072643, 1,
        1.239050, 1.088723, 1.13045, 1.056176,
        1, 1.0806722, 1, 1,
        1, 1, 1, 1,
        1, 1, 1, 1,
        1, 1};

    // Define pad IDs (should match the size of Eloss and scalingFactors)
    const int nPads = 30; // Assuming 30 pads, consistent with array sizes
    TString branchPrefix = "IC_C_Cor_";
    TString branchNames[nPads];
    Double_t IC_C_Values[nPads];

    // New arrays for calibrated energy branches
    TString calEnergyBranchPrefix = "IC_E_Cal_";
    TString calEnergyBranchNames[nPads];
    Double_t IC_E_Cal_Values[nPads]; // To store the calculated energyMod values

    for (int i = 0; i < nPads; ++i)
    {
        branchNames[i] = branchPrefix + std::to_string(i);
        calEnergyBranchNames[i] = calEnergyBranchPrefix + std::to_string(i);
    }

    // Open ROOT file and get TTree
    // TString inputFile = "/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_1053new100725.root";
    TString inputFile = "/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_" + std::to_string(RUN_Nbr) + "newSpline.root";
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

    Long64_t nEntries = tree->GetEntries();
    std::cout << "Number of entries: " << nEntries << std::endl;

    // Create histograms
    TCanvas *c1 = new TCanvas("c1", "IC Plots", 1500, 1400);
    c1->Divide(2, 2);

    TList *hlist = new TList();

    TH2F *h2D_SLOW_BEAM_IC_Charge = new TH2F("h2D_SLOW_BEAM_IC_Charge",
                                             "SLOW BEAM IC Charge vs IC Pad; Pad; IC Charge [au]",
                                             31, -0.5, 30.5, 1000, 0, 3000);
    hlist->Add(h2D_SLOW_BEAM_IC_Charge);

    TH2F *h2D_SLOW_BEAM_IC_Energy_Calibrated = new TH2F("h2D_SLOW_BEAM_IC_Energy_Calibrated",
                                                        "SLOW BEAM IC Energy vs IC Pad; Pad; IC Energy [MeV/u]",
                                                        31, -0.5, 30.5, 1000, 0, 100);
    hlist->Add(h2D_SLOW_BEAM_IC_Energy_Calibrated);

    TH2F *h2D_SLOW_BEAM_IC_Energy_Calibrated_2ndMod = new TH2F("h2D_SLOW_BEAM_IC_Energy_Calibrated_2ndMod",
                                                               "SLOW BEAM IC Energy vs IC Pad 2nd Modification; Pad; IC Energy [MeV/u]",
                                                               31, -0.5, 30.5, 1000, 0, 100);
    hlist->Add(h2D_SLOW_BEAM_IC_Energy_Calibrated_2ndMod);

    // Clone tree structure
    TTree *newTree = tree->CloneTree(0);

    for (int i = 0; i < nPads; ++i)
    {
        tree->SetBranchAddress(branchNames[i], &IC_C_Values[i]);
        newTree->SetBranchAddress(calEnergyBranchNames[i], &IC_E_Cal_Values[i]);
    }

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        if (i % 100000 == 0)
            std::cout << "Processing event: " << i << std::endl;

        tree->GetEntry(i);

        for (int padID = 0; padID < nPads; ++padID)
        {
            double val = IC_C_Values[padID];
            double energy = val * scalingFactors[padID];
            double energyMod = val * scalingFactors[padID] * scalingFactorsMod[padID];

            if (val > 0)
            {
                IC_E_Cal_Values[padID] = energyMod;

                h2D_SLOW_BEAM_IC_Charge->Fill(padID, val);
                h2D_SLOW_BEAM_IC_Energy_Calibrated->Fill(padID, energy);
                h2D_SLOW_BEAM_IC_Energy_Calibrated_2ndMod->Fill(padID, energyMod);
            }
            else
            {
                IC_E_Cal_Values[padID] = -1.0e6;
            }
        }

        newTree->Fill();
    }

    // Draw plots
    c1->cd(1);
    h2D_SLOW_BEAM_IC_Charge->SetStats(0);
    h2D_SLOW_BEAM_IC_Charge->Draw("COLZ");

    c1->cd(2);
    h2D_SLOW_BEAM_IC_Energy_Calibrated->SetStats(0);
    h2D_SLOW_BEAM_IC_Energy_Calibrated->Draw("COLZ");

    c1->cd(3);
    h2D_SLOW_BEAM_IC_Energy_Calibrated_2ndMod->SetStats(0);
    h2D_SLOW_BEAM_IC_Energy_Calibrated_2ndMod->Draw("COLZ");

    c1->Update();

    // Save histogram to a new ROOT file
    TString outputFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/SLOW_Beam_1053_2024_hist1107.root";
    TFile *outFile = new TFile(outputFileName, "RECREATE");
    if (!outFile || outFile->IsZombie())
    {
        std::cerr << "Error creating output file: " << outputFileName << std::endl;
        file->Close();
        return 1;
    }

    c1->Write();
    hlist->Write("HistogramAndGraphLists", TObject::kSingleKey);
    outFile->Close();

    // Replace old tree in original file
    file->cd();
    tree->Delete(); // Delete old tree
    // newTree->Write("tree_new"); // Write updated tree TObject::kOverwrite
    newTree->Write("tree_new", TObject::kOverwrite); // Write updated tree
    std::cout << "Updated TTree written to '" << inputFile << "'" << std::endl;

    // std::cout << "Press enter to exit..." << std::endl;
    // std::cin.get();

    file->Close();
    delete c1;
    delete hlist;

    return 0;
}
