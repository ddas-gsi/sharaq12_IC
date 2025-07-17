# // Use this C++ code fastBeamChargePlot.C .. Not the fastBeamChargePlot.py file.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import product
from datetime import datetime
from multiprocessing import Pool
from functools import partial
import os
import csv
import ROOT


def generateICBranchNameList(ICSegmentList, branchNameStr):
    """Returns the IC Branch Names as a List of strings"""
    ICBranchNameList = []
    for segmentNbr in ICSegmentList:
        ICBranchName = branchNameStr + str(segmentNbr)
        ICBranchNameList.append(ICBranchName)
        # print(ICBranchName)
    return ICBranchNameList


# Define the IC Branch Names
ICSegmentList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
IC_C_BranchNameList = generateICBranchNameList(ICSegmentList, "IC_C_")                              # Raw IC Charge Branch Names

print(IC_C_BranchNameList)

# HISTOGRAMS
c1 = ROOT.TCanvas("c1", "IC Charge Correction", 1500, 600)
# c1.Divide(2, 1)
h2D_FAST_BEAM_IC_Charge = ROOT.TH2F("h2D_FAST_BEAM_IC_Charge", "FAST BEAM IC Charge vs IC Pad; Pad; IC Charge ", 31, -0.5, 30.5, 1000, 0, 3000)


inputFile = f"/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/opticsOutput/FastBeam_0089_2024.root"
rootFile = ROOT.TFile(inputFile, "UPDATE")
tree = rootFile.Get("tree_new")
num_entries = tree.GetEntries()
# num_entries = 10000
print(f"Number of entries in the TTree: {num_entries}")


for i_Event in range(num_entries):        
    tree.GetEntry(i_Event)

    if (i_Event % 100000 == 0):
        print(f"Event Number: {i_Event}")

    for padID in ICSegmentList:
        IC_C_BranchName = IC_C_BranchNameList[padID]
        IC_C_padID_Charge = getattr(tree, IC_C_BranchName)

        if IC_C_padID_Charge > 0:
            h2D_FAST_BEAM_IC_Charge.Fill(padID, IC_C_padID_Charge)

# Draw the 2D Histogram
c1.cd()
h2D_FAST_BEAM_IC_Charge.SetStats(0)
h2D_FAST_BEAM_IC_Charge.Draw("COLZ")
