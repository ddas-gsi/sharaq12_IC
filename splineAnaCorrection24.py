#!/usr/bin/env python3

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
import sys



avg3padE_param = [20, 166.591, -7863.94, 97677.6, 0, 0]   # avg3padE vs betaS1_0**2
# splinePeakE_param = [38, 38.4398, -0.0303254, -0.000167409, -6.7234e-07, 1.03274e-08]  # splinePeakE vs S1_X_0
splinePeakE_param = [38, 38.5555, -0.00679176, 0.000213653, -8.43577e-08, -7.59006e-09]  # splinePeakE vs S1_X_0
splineRange_param = [150.0, -810.934, 5401.56, 0, 0, 0]  # splineRange/betaS1_0**2/100 vs betaS1_0


def correction_factor(x, param):
    par = param[1:]
    y = (par[0] + par[1] * x + par[2] * x**2 + par[3] * x**3 + par[4] * x**4)
    meanY = np.full_like(x, param[0])
    dy = meanY -y
    return dy

def main(RUN_Nbr):

    print(f"Processing RUN Number: {RUN_Nbr} ...")
    print(f"Copying the root file {RUN_Nbr}_Spline_2024.root to the splineCorrection folder ...")

    os.system(f"cp /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/{RUN_Nbr}_Spline_2024.root \
              /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/{RUN_Nbr}_Spline_2024.root")
    
    inputFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/" + str(RUN_Nbr) + "_Spline_2024.root"
    rootFile = ROOT.TFile.Open(inputFileName, "UPDATE")
    tree = rootFile.Get("tree_new")
    num_entries = tree.GetEntries()
    # num_entries = 1000  # For testing purposes, limit to 1000 entries
    print(f"Number of entries in the TTree: {num_entries}")

    # Enable all branches first
    tree.SetBranchStatus("*", 1)

    # define new variables
    a3EBetaS1 = np.zeros(1, dtype=np.float64)
    sPeakES1X = np.zeros(1, dtype=np.float64)
    sRangeBetaS1 = np.zeros(1, dtype=np.float64)

    # Create new branches for the new variables
    a3EBetaS1_branch = tree.Branch("a3EBetaS1", a3EBetaS1, "a3EBetaS1/D")
    sPeakES1X_branch = tree.Branch("sPeakES1X", sPeakES1X, "sPeakES1X/D")
    sRangeBetaS1_branch = tree.Branch("sRangeBetaS1", sRangeBetaS1, "sRangeBetaS1/D")

    countGoodEvent = 0

    # Loop over the entries in the tree
    for i_Event in range(num_entries):
        tree.GetEntry(i_Event)

        if (i_Event % 100000 == 0):
            print(f"Event Number: {i_Event}")

        # Define the branches
        FE9_X_0 = getattr(tree, "FE9_X_0")
        PID_T_0 = getattr(tree, "PID_T_0")
        S1_X_0 = getattr(tree, "S1_X_0")
        S1PID_T_0 = getattr(tree, "S1PID_T_0")
        S1_E_0 = getattr(tree, "S1_E_0")
        S1_Y_0 = getattr(tree, "S1_Y_0")
        S1_A_0 = getattr(tree, "S1_A_0")
        S1_B_0 = getattr(tree, "S1_B_0")
        betaS1_0 = getattr(tree, "betaS1_0")
        avg3padE = getattr(tree, "avg3padE")
        splinePeakE = getattr(tree, "splinePeakE")
        splineRange = getattr(tree, "splineRange")
        


        # Initialize new variables to -1e6
        a3EBetaS1[0] = -1e6
        sPeakES1X[0] = -1e6
        sRangeBetaS1[0] = -1e6

        if betaS1_0 > 0:
            countGoodEvent += 1

            # Calculate the correction factors
            a3EBetaS1[0] = correction_factor(betaS1_0**2, avg3padE_param) + avg3padE
            sPeakES1X[0] = correction_factor(S1_X_0, splinePeakE_param) + splinePeakE
            sRangeBetaS1[0] = correction_factor(betaS1_0, splineRange_param) + splineRange/ betaS1_0**2 / 100.0

        
        
        # ::::: Fill the branches with the calculated values :::::
        a3EBetaS1_branch.Fill()
        sPeakES1X_branch.Fill()
        sRangeBetaS1_branch.Fill()  

    # Write the updated tree to the file    
    rootFile.cd()
    rootFile.Write("", ROOT.TObject.kOverwrite)
    rootFile.Close()
    print(f"Total Number of Events processed: {countGoodEvent}")


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python3 splineCorrection24.py RUN_Nbr")
        sys.exit(1)

    RUN_Nbr = int(sys.argv[1])

    # RUN_Nbr = 1052
    start_time = datetime.now()
    main(RUN_Nbr)
    end_time = datetime.now()
    print(f"Time taken for processing: {end_time - start_time}")


    
