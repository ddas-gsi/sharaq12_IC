// Use this C++ code icChargeCorrection.C .. Not the icChargeCorrection.py file.

#include <iostream>
#include <vector>
#include <string>
#include <array>     // For fixed-size arrays like pad_ID
#include <cmath>     // For pow, full_like equivalent
#include <algorithm> // For std::fill

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObject.h" // For kOverwrite

// Pad correction parameters (using std::array for fixed size)
// The first element is the mean charge, followed by the polynomial coefficients.
std::array<double, 7> pad_0 = {982.2, 1023.8, -2.29364, -0.0365386, -0.000195611, 4.83869e-06, 9.0028e-08};
std::array<double, 7> pad_1 = {934.4, 973.409, -1.99299, -0.0336818, 6.20176e-06, 4.23026e-06, 1.71449e-08};
std::array<double, 7> pad_2 = {1018.0, 1054.38, -2.2641, -0.0295328, 1.74352e-05, 3.62653e-06, 9.09791e-09};
std::array<double, 7> pad_3 = {1060.0, 1095.58, -2.51861, -0.0272889, 0.000103151, 3.33841e-06, -8.39731e-09};
std::array<double, 7> pad_4 = {1343.0, 1390.09, -3.5073, -0.0219948, 0.000384003, 1.37212e-06, -6.33711e-08};
std::array<double, 7> pad_5 = {1394.0, 1436.05, -3.29489, -0.0276866, 7.1602e-05, 2.76825e-06, -6.84038e-09};
std::array<double, 7> pad_6 = {1170.0, 1207.79, -2.49011, -0.028619, -0.000186802, 3.69802e-06, 4.02745e-08};
std::array<double, 7> pad_7 = {1108.0, 1145.95, -2.37942, -0.0254606, -8.36485e-05, 2.46088e-06, 2.3824e-08};
std::array<double, 7> pad_8 = {1273.0, 1310.89, -2.51221, -0.0253169, -0.000191775, 2.87244e-06, 3.76846e-08};
std::array<double, 7> pad_9 = {1279.0, 1314.51, -2.45355, -0.0248844, -0.000129802, 3.21088e-06, 2.62168e-08};
std::array<double, 7> pad_10 = {1288.0, 1320.78, -2.49497, -0.0263923, -1.98699e-05, 4.13981e-06, 8.82406e-09};
std::array<double, 7> pad_11 = {1247.0, 1284.18, -2.19985, -0.0233099, -0.000112906, 2.45282e-06, 2.56592e-08};
std::array<double, 7> pad_12 = {1281.0, 1322.75, -1.83127, -0.0287894, -0.000283749, 3.49972e-06, 5.85403e-08};
std::array<double, 7> pad_13 = {1262.0, 1305.65, -1.78447, -0.0244827, -0.000382759, 2.27385e-06, 7.822e-08};
std::array<double, 7> pad_14 = {1330.0, 1372.36, -1.61442, -0.0247793, -0.000650052, 3.17676e-06, 1.31614e-07};
std::array<double, 7> pad_15 = {1349.0, 1388.48, -1.28014, -0.0144987, -0.000875432, 9.93095e-07, 1.61077e-07};
std::array<double, 7> pad_16 = {1338.0, 1389.58, -1.78708, -0.0257677, -0.000447052, 3.63395e-06, 8.38101e-08};
std::array<double, 7> pad_17 = {1359.0, 1428.52, -2.41925, -0.0276335, -0.000235782, 3.85709e-06, 6.92493e-08};
std::array<double, 7> pad_18 = {1336.0, 1406.13, -1.41996, -0.0314169, -0.000548039, 5.19766e-06, 4.17918e-08};
std::array<double, 7> pad_19 = {1530.0, 1571.08, -2.40544, 0.0125569, -0.000505481, 1.81144e-06, 1.58637e-07};
std::array<double, 7> pad_20 = {1447.0, 1637.04, 0.276491, -0.0644718, -0.00299767, -1.3076e-06, 5.94385e-07};
std::array<double, 7> pad_21 = {1024.0, 1090.41, -0.829893, -0.0846117, 0.0, 0.0, 0.0};
std::array<double, 7> pad_22 = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 7> pad_23 = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 7> pad_24 = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 7> pad_25 = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 7> pad_26 = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 7> pad_27 = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 7> pad_28 = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 7> pad_29 = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

// Array of all pad parameters for easier access
std::array<std::array<double, 7>, 30> pad_params = {
    pad_0, pad_1, pad_2, pad_3, pad_4, pad_5, pad_6, pad_7, pad_8, pad_9,
    pad_10, pad_11, pad_12, pad_13, pad_14, pad_15, pad_16, pad_17, pad_18, pad_19,
    pad_20, pad_21, pad_22, pad_23, pad_24, pad_25, pad_26, pad_27, pad_28, pad_29};

// Function for the polynomial fit (S1Y correction function)
double padIDCharge_S1Y_correction_function(double x, const std::array<double, 7> &pad_ID)
{
    // Coefficients start from index 1 in the pad_ID array
    double y = pad_ID[1] + pad_ID[2] * x + pad_ID[3] * std::pow(x, 2) + pad_ID[4] * std::pow(x, 3) + pad_ID[5] * std::pow(x, 4) + pad_ID[6] * std::pow(x, 5);
    return y;
}

// Function to get the mean charge (the first element in pad_ID)
double fn_meanCharge(double x, const std::array<double, 7> &pad_ID)
{
    // In C++, for a single value, you just return it. No need for full_like.
    return pad_ID[0];
}

// Function to calculate the correction factor
double padIDCharge_S1Y_correction_factor(double x, const std::array<double, 7> &pad_ID)
{
    double dy = fn_meanCharge(x, pad_ID) - padIDCharge_S1Y_correction_function(x, pad_ID);
    return dy;
}

void icChargeCorrection(int RUN_Nbr)
{
    // int RUN_Nbr = 1052;

    // std::cout << "Copying the root file 24sharaq12phys_" << RUN_Nbr << "new100725.root to the physOutput folder ..." << std::endl;
    std::cout << "Copying the root file 24sharaq12phys_" << RUN_Nbr << "newSpline.root to the physOutput folder ..." << std::endl;
    std::string copy_command = "cp /u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_" + std::to_string(RUN_Nbr) + "new.root /u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_" + std::to_string(RUN_Nbr) + "newSpline.root";
    std::system(copy_command.c_str());

    std::string inputFile = "/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_" + std::to_string(RUN_Nbr) + "newSpline.root";
    TFile *rootFile = TFile::Open(inputFile.c_str(), "UPDATE");
    if (!rootFile || rootFile->IsZombie())
    {
        std::cerr << "Error opening ROOT file: " << inputFile << std::endl;
        return;
    }

    TTree *tree = static_cast<TTree *>(rootFile->Get("tree_new"));
    if (!tree)
    {
        std::cerr << "Error getting TTree 'tree_new' from file." << std::endl;
        rootFile->Close();
        return;
    }

    Long64_t num_entries = tree->GetEntries();
    // num_entries = 10000; // For testing
    std::cout << "Number of entries in the TTree: " << num_entries << std::endl;

    // Declare variables for existing branches (input)
    Double_t S1_Y_0;
    std::array<Double_t, 30> IC_C_raw; // Array to hold raw IC charges

    // Declare variables for new branches (output)
    std::array<Double_t, 30> IC_C_Cor; // Array to hold corrected IC charges

    // Set branch addresses for existing branches
    tree->SetBranchAddress("S1_Y_0", &S1_Y_0);
    for (int i = 0; i < 30; ++i)
    {
        std::string branch_name = "IC_C_" + std::to_string(i);
        tree->SetBranchAddress(branch_name.c_str(), &IC_C_raw[i]);
    }

    // Create new branches for corrected charges
    std::array<TBranch *, 30> IC_C_Cor_br;
    for (int i = 0; i < 30; ++i)
    {
        std::string branch_name = "IC_C_Cor_" + std::to_string(i);
        std::string branch_type = "IC_C_Cor_" + std::to_string(i) + "/D"; // /D for Double_t
        IC_C_Cor_br[i] = tree->Branch(branch_name.c_str(), &IC_C_Cor[i], branch_type.c_str());
    }

    // HISTOGRAMS
    TCanvas *c1 = new TCanvas("c1", "IC Charge Correction", 1500, 600);
    c1->Divide(2, 1);
    TH2F *h2D_IC_Charge = new TH2F("h2D_IC_Charge", "IC Charge vs IC Pad; Pad; IC Charge ", 31, -0.5, 30.5, 1000, 0, 3000);
    TH2F *h2D_IC_Charge_Corrected = new TH2F("h2D_IC_Charge_Corrected", "IC Charge Corrected vs IC Pad; Pad; IC Charge ", 31, -0.5, 30.5, 1000, 0, 3000);

    // Create a TList
    TList *hlist = new TList();
    hlist->Add(h2D_IC_Charge);
    hlist->Add(h2D_IC_Charge_Corrected);

    // Loop over entries
    for (Long64_t i_Event = 0; i_Event < num_entries; ++i_Event)
    {
        tree->GetEntry(i_Event);

        if (i_Event % 100000 == 0)
        {
            std::cout << "Event Number: " << i_Event << std::endl;
        }

        // Initialize new branches to -1.0e6
        std::fill(IC_C_Cor.begin(), IC_C_Cor.end(), -1.0e6);

        if (S1_Y_0 > -1.0e5)
        { // Check if S1_Y_0 is valid
            for (int i = 0; i < 22; ++i)
            { // Pads 0 to 21 get correction
                if (IC_C_raw[i] > 350.0)
                {
                    IC_C_Cor[i] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_params[i]) + IC_C_raw[i];
                }
                else
                {
                    IC_C_Cor[i] = 0;
                }
            }
            // Pads 22 to 29 do not get correction (just threshold)
            for (int i = 22; i < 30; ++i)
            {
                if (IC_C_raw[i] > 350.0)
                {
                    IC_C_Cor[i] = IC_C_raw[i];
                }
                else
                {
                    IC_C_Cor[i] = 0;
                }
            }
        }

        // Fill the 2D Histograms
        for (int i = 0; i < 30; ++i)
        {
            h2D_IC_Charge->Fill(i, IC_C_raw[i]);
            h2D_IC_Charge_Corrected->Fill(i, IC_C_Cor[i]);
        }

        // Fill the new branches
        for (int i = 0; i < 30; ++i)
        {
            IC_C_Cor_br[i]->Fill();
        }
    }

    // Draw the histograms
    c1->cd(1);
    h2D_IC_Charge->Draw("colz");

    c1->cd(2);
    h2D_IC_Charge_Corrected->Draw("colz");

    c1->Update();
    // c1->Draw(); // Not strictly needed to see in C++ interactive mode

    // Save the Histograms to the root file
    std::string histFileName = "/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_" + std::to_string(RUN_Nbr) + "histFileSpline.root";
    TFile *histFile = TFile::Open(histFileName.c_str(), "RECREATE");
    if (!histFile || histFile->IsZombie())
    {
        std::cerr << "Error creating histogram file: " << histFileName << std::endl;
    }
    else
    {
        c1->Write();
        hlist->Write("HistogramLists", TObject::kSingleKey);
        histFile->Close();
        delete histFile;
    }

    // Save the values to the root file (write the updated tree)
    rootFile->cd();
    tree->Write("", TObject::kOverwrite); // Overwrite the existing tree

    // Clean up
    rootFile->Close();
    delete rootFile;
    delete c1;
    delete h2D_IC_Charge;
    delete h2D_IC_Charge_Corrected;
    delete hlist;

    std::cout << "Processing complete." << std::endl;
}

// // Main function to execute the processing
// int main()
// {
//     process_root_file();
//     return 0;
// }