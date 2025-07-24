#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <stdio.h>
#include <math.h>
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
const int nave = 5;
void test_fit()
{
    TFile *file = TFile::Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root");
    TTree *tree = (TTree *)file->Get("tree_new");
    TH2F *h = new TH2F("h", "h", 30, 0, 30, 100000, -100, 100);
    vector<vector<double>> IC_E_Cal_vec;
    Double_t IC_E_Cal[30];

    for (int k = 0; k < 30; k++)
    {
        string IC_E_SetBranch = "IC_E_Cal_" + to_string(k);
        tree->SetBranchAddress(IC_E_SetBranch.c_str(), &IC_E_Cal[k]); // .c_str() converts any string to CONSTANT
    }

    int entries = tree->GetEntries();
    cout << "Number of Entries: " << entries << endl;
    // entries = 1000;
    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        if (IC_E_Cal[0] < 10)
            continue;
        for (int k = 0; k < 30 - nave; k++)
        {
            double ave = 0;
            for (int j = 0; j < nave; j++)
            {
                ave += IC_E_Cal[k + j];
            }
            ave /= nave;
            h->Fill(k, IC_E_Cal[k] - ave);
        }
        if (i % 10000 == 0)
        {
            cout << i << "/" << entries << " done ..." << flush;
        }
    }
    TFile *fout = new TFile("fout.root", "RECREATE");
    fout->cd();
    h->Write("", TObject::kOverwrite);
    fout->Close();
}