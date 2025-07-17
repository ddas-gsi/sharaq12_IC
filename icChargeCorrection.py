# Use this C++ code icChargeCorrection.C .. Not the icChargeCorrection.py file.

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


pad_0 = [982.2, 1023.8, -2.29364,  -0.0365386, -0.000195611,  4.83869e-06 , 9.0028e-08]
pad_1 = [934.4, 973.409, -1.99299, -0.0336818, 6.20176e-06, 4.23026e-06 , 1.71449e-08]
pad_2 = [1018.0, 1054.38, -2.2641, -0.0295328, 1.74352e-05, 3.62653e-06, 9.09791e-09]
pad_3 = [1060.0, 1095.58, -2.51861, -0.0272889, 0.000103151, 3.33841e-06, -8.39731e-09]
pad_4 = [1343.0, 1390.09, -3.5073, -0.0219948, 0.000384003, 1.37212e-06, -6.33711e-08]
pad_5 = [1394.0, 1436.05, -3.29489, -0.0276866, 7.1602e-05, 2.76825e-06, -6.84038e-09]
pad_6 = [1170.0, 1207.79, -2.49011, -0.028619, -0.000186802, 3.69802e-06, 4.02745e-08]
pad_7 = [1108.0, 1145.95, -2.37942, -0.0254606, -8.36485e-05, 2.46088e-06, 2.3824e-08]
pad_8 = [1273.0, 1310.89, -2.51221, -0.0253169, -0.000191775, 2.87244e-06, 3.76846e-08]
pad_9 = [1279.0, 1314.51, -2.45355, -0.0248844, -0.000129802, 3.21088e-06, 2.62168e-08]
pad_10 = [1288.0, 1320.78, -2.49497, -0.0263923, -1.98699e-05, 4.13981e-06, 8.82406e-09]
pad_11 = [1247.0, 1284.18, -2.19985, -0.0233099, -0.000112906, 2.45282e-06, 2.56592e-08]
pad_12 = [1281.0, 1322.75, -1.83127, -0.0287894, -0.000283749, 3.49972e-06, 5.85403e-08]
pad_13 = [1262.0, 1305.65, -1.78447, -0.0244827, -0.000382759, 2.27385e-06, 7.822e-08]
pad_14 = [1330.0, 1372.36, -1.61442, -0.0247793, -0.000650052, 3.17676e-06, 1.31614e-07]
pad_15 = [1349.0, 1388.48, -1.28014, -0.0144987, -0.000875432, 9.93095e-07, 1.61077e-07]
pad_16 = [1338.0, 1389.58, -1.78708, -0.0257677, -0.000447052, 3.63395e-06, 8.38101e-08]
pad_17 = [1359.0, 1428.52, -2.41925, -0.0276335, -0.000235782, 3.85709e-06, 6.92493e-08]
pad_18 = [1336.0, 1406.13, -1.41996, -0.0314169, -0.000548039, 5.19766e-06, 4.17918e-08]
pad_19 = [1530.0, 1571.08, -2.40544, 0.0125569, -0.000505481, 1.81144e-06, 1.58637e-07]
pad_20 = [1447.0, 1637.04, 0.276491, -0.0644718, -0.00299767, -1.3076e-06, 5.94385e-07]
pad_21 = [1024.0, 1090.41, -0.829893, -0.0846117, 0.0, 0.0, 0.0]
# pad_22 = [183.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# pad_23 = [179.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# pad_24 = [167.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# pad_25 = [171.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# pad_26 = [169.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# pad_27 = [169.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# pad_28 = [162.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# pad_29 = [168.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pad_22 = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pad_23 = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pad_24 = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pad_25 = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pad_26 = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pad_27 = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pad_28 = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pad_29 = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


def generateICBranchNameList(ICSegmentList, branchNameStr):
    """Returns the IC Branch Names as a List of strings"""
    ICBranchNameList = []
    for segmentNbr in ICSegmentList:
        ICBranchName = branchNameStr + str(segmentNbr)
        ICBranchNameList.append(ICBranchName)
        # print(ICBranchName)
    return ICBranchNameList

def padIDCharge_S1Y_correction_function(x, pad_ID):
    par = pad_ID[1:]
    y = (par[0] + par[1] * x + par[2] * x**2 + par[3] * x**3 + par[4] * x**4 + par[5] * x**5)
    return y

def fn_meanCharge(x, pad_ID):
    y = np.full_like(x, pad_ID[0])
    return y

def padIDCharge_S1Y_correction_factor(x, pad_ID):
    dy = fn_meanCharge(x, pad_ID) - padIDCharge_S1Y_correction_function(x, pad_ID)
    return dy    
    
# # HISTOGRAMS
c1 = ROOT.TCanvas("c1", "IC Charge Correction", 1500, 600)
c1.Divide(2, 1)
h2D_IC_Charge = ROOT.TH2F("h2D_IC_Charge", "IC Charge vs IC Pad; Pad; IC Charge ", 31, -0.5, 30.5, 1000, 0, 3000)
h2D_IC_Charge_Corrected = ROOT.TH2F("h2D_IC_Charge_Corrected", "IC Charge Corrected vs IC Pad; Pad; IC Charge ", 31, -0.5, 30.5, 1000, 0, 3000)

# Create a TList
hlist = ROOT.TList()
hlist.Add(h2D_IC_Charge)
hlist.Add(h2D_IC_Charge_Corrected)



# Define the IC Branch Names
ICSegmentList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
IC_C_BranchNameList = generateICBranchNameList(ICSegmentList, "IC_C_")                              # Raw IC Charge Branch Names
IC_C_Corrected_BranchNameList = generateICBranchNameList(ICSegmentList, "IC_C_Cor_")          # Corrected IC Charge Branch Names


# Define variables for the new branches
IC_C_Cor_0 = np.zeros(1, dtype=np.float64)
IC_C_Cor_1 = np.zeros(1, dtype=np.float64)
IC_C_Cor_2 = np.zeros(1, dtype=np.float64)
IC_C_Cor_3 = np.zeros(1, dtype=np.float64)
IC_C_Cor_4 = np.zeros(1, dtype=np.float64)
IC_C_Cor_5 = np.zeros(1, dtype=np.float64)
IC_C_Cor_6 = np.zeros(1, dtype=np.float64)
IC_C_Cor_7 = np.zeros(1, dtype=np.float64)
IC_C_Cor_8 = np.zeros(1, dtype=np.float64)
IC_C_Cor_9 = np.zeros(1, dtype=np.float64)
IC_C_Cor_10 = np.zeros(1, dtype=np.float64)
IC_C_Cor_11 = np.zeros(1, dtype=np.float64)
IC_C_Cor_12 = np.zeros(1, dtype=np.float64)
IC_C_Cor_13 = np.zeros(1, dtype=np.float64)
IC_C_Cor_14 = np.zeros(1, dtype=np.float64)
IC_C_Cor_15 = np.zeros(1, dtype=np.float64)
IC_C_Cor_16 = np.zeros(1, dtype=np.float64)
IC_C_Cor_17 = np.zeros(1, dtype=np.float64)
IC_C_Cor_18 = np.zeros(1, dtype=np.float64)
IC_C_Cor_19 = np.zeros(1, dtype=np.float64)
IC_C_Cor_20 = np.zeros(1, dtype=np.float64)
IC_C_Cor_21 = np.zeros(1, dtype=np.float64)
IC_C_Cor_22 = np.zeros(1, dtype=np.float64)
IC_C_Cor_23 = np.zeros(1, dtype=np.float64)
IC_C_Cor_24 = np.zeros(1, dtype=np.float64)
IC_C_Cor_25 = np.zeros(1, dtype=np.float64)
IC_C_Cor_26 = np.zeros(1, dtype=np.float64)
IC_C_Cor_27 = np.zeros(1, dtype=np.float64)
IC_C_Cor_28 = np.zeros(1, dtype=np.float64)
IC_C_Cor_29 = np.zeros(1, dtype=np.float64)


# # Read the root file
RUN_Nbr = 1053

print(f"Copying the root file 24sharaq12phys_{RUN_Nbr}new100725.root to the physOutput folder ...")
os.system(f"cp /u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_{RUN_Nbr}new.root /u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_{RUN_Nbr}new100725.root")

inputFile = f"/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_{RUN_Nbr}new100725.root"
rootFile = ROOT.TFile(inputFile, "UPDATE")
tree = rootFile.Get("tree_new")
num_entries = tree.GetEntries()
# num_entries = 10000
print(f"Number of entries in the TTree: {num_entries}")


# Add the new branches to the output tree
IC_C_Cor_0_br = tree.Branch("IC_C_Cor_0", IC_C_Cor_0, "IC_C_Cor_0/D")
IC_C_Cor_1_br = tree.Branch("IC_C_Cor_1", IC_C_Cor_1, "IC_C_Cor_1/D")
IC_C_Cor_2_br = tree.Branch("IC_C_Cor_2", IC_C_Cor_2, "IC_C_Cor_2/D")
IC_C_Cor_3_br = tree.Branch("IC_C_Cor_3", IC_C_Cor_3, "IC_C_Cor_3/D")
IC_C_Cor_4_br = tree.Branch("IC_C_Cor_4", IC_C_Cor_4, "IC_C_Cor_4/D")
IC_C_Cor_5_br = tree.Branch("IC_C_Cor_5", IC_C_Cor_5, "IC_C_Cor_5/D")
IC_C_Cor_6_br = tree.Branch("IC_C_Cor_6", IC_C_Cor_6, "IC_C_Cor_6/D")
IC_C_Cor_7_br = tree.Branch("IC_C_Cor_7", IC_C_Cor_7, "IC_C_Cor_7/D")
IC_C_Cor_8_br = tree.Branch("IC_C_Cor_8", IC_C_Cor_8, "IC_C_Cor_8/D")
IC_C_Cor_9_br = tree.Branch("IC_C_Cor_9", IC_C_Cor_9, "IC_C_Cor_9/D")
IC_C_Cor_10_br = tree.Branch("IC_C_Cor_10", IC_C_Cor_10, "IC_C_Cor_10/D")
IC_C_Cor_11_br = tree.Branch("IC_C_Cor_11", IC_C_Cor_11, "IC_C_Cor_11/D")
IC_C_Cor_12_br = tree.Branch("IC_C_Cor_12", IC_C_Cor_12, "IC_C_Cor_12/D")
IC_C_Cor_13_br = tree.Branch("IC_C_Cor_13", IC_C_Cor_13, "IC_C_Cor_13/D")
IC_C_Cor_14_br = tree.Branch("IC_C_Cor_14", IC_C_Cor_14, "IC_C_Cor_14/D")
IC_C_Cor_15_br = tree.Branch("IC_C_Cor_15", IC_C_Cor_15, "IC_C_Cor_15/D")
IC_C_Cor_16_br = tree.Branch("IC_C_Cor_16", IC_C_Cor_16, "IC_C_Cor_16/D")
IC_C_Cor_17_br = tree.Branch("IC_C_Cor_17", IC_C_Cor_17, "IC_C_Cor_17/D")
IC_C_Cor_18_br = tree.Branch("IC_C_Cor_18", IC_C_Cor_18, "IC_C_Cor_18/D")
IC_C_Cor_19_br = tree.Branch("IC_C_Cor_19", IC_C_Cor_19, "IC_C_Cor_19/D")
IC_C_Cor_20_br = tree.Branch("IC_C_Cor_20", IC_C_Cor_20, "IC_C_Cor_20/D")
IC_C_Cor_21_br = tree.Branch("IC_C_Cor_21", IC_C_Cor_21, "IC_C_Cor_21/D")
IC_C_Cor_22_br = tree.Branch("IC_C_Cor_22", IC_C_Cor_22, "IC_C_Cor_22/D")
IC_C_Cor_23_br = tree.Branch("IC_C_Cor_23", IC_C_Cor_23, "IC_C_Cor_23/D")
IC_C_Cor_24_br = tree.Branch("IC_C_Cor_24", IC_C_Cor_24, "IC_C_Cor_24/D")
IC_C_Cor_25_br = tree.Branch("IC_C_Cor_25", IC_C_Cor_25, "IC_C_Cor_25/D")
IC_C_Cor_26_br = tree.Branch("IC_C_Cor_26", IC_C_Cor_26, "IC_C_Cor_26/D")
IC_C_Cor_27_br = tree.Branch("IC_C_Cor_27", IC_C_Cor_27, "IC_C_Cor_27/D")
IC_C_Cor_28_br = tree.Branch("IC_C_Cor_28", IC_C_Cor_28, "IC_C_Cor_28/D")
IC_C_Cor_29_br = tree.Branch("IC_C_Cor_29", IC_C_Cor_29, "IC_C_Cor_29/D")


for i_Event in range(num_entries):        
    tree.GetEntry(i_Event)

    if (i_Event % 10000 == 0):
        print(f"Event Number: {i_Event}")

    # Define the branches
    S1_Y_0 = getattr(tree, "S1_Y_0")

    IC_C_0 = getattr(tree, "IC_C_0")
    IC_C_1 = getattr(tree, "IC_C_1")
    IC_C_2 = getattr(tree, "IC_C_2")
    IC_C_3 = getattr(tree, "IC_C_3")
    IC_C_4 = getattr(tree, "IC_C_4")
    IC_C_5 = getattr(tree, "IC_C_5")
    IC_C_6 = getattr(tree, "IC_C_6")
    IC_C_7 = getattr(tree, "IC_C_7")
    IC_C_8 = getattr(tree, "IC_C_8")
    IC_C_9 = getattr(tree, "IC_C_9")
    IC_C_10 = getattr(tree, "IC_C_10")
    IC_C_11 = getattr(tree, "IC_C_11")
    IC_C_12 = getattr(tree, "IC_C_12")
    IC_C_13 = getattr(tree, "IC_C_13")
    IC_C_14 = getattr(tree, "IC_C_14")
    IC_C_15 = getattr(tree, "IC_C_15")
    IC_C_16 = getattr(tree, "IC_C_16")
    IC_C_17 = getattr(tree, "IC_C_17")
    IC_C_18 = getattr(tree, "IC_C_18")
    IC_C_19 = getattr(tree, "IC_C_19")
    IC_C_20 = getattr(tree, "IC_C_20")
    IC_C_21 = getattr(tree, "IC_C_21")
    IC_C_22 = getattr(tree, "IC_C_22")
    IC_C_23 = getattr(tree, "IC_C_23")
    IC_C_24 = getattr(tree, "IC_C_24")
    IC_C_25 = getattr(tree, "IC_C_25")
    IC_C_26 = getattr(tree, "IC_C_26")
    IC_C_27 = getattr(tree, "IC_C_27")
    IC_C_28 = getattr(tree, "IC_C_28")
    IC_C_29 = getattr(tree, "IC_C_29")

    # Initialize new branches to -1.0e6
    IC_C_Cor_0[0] = -1.0e6
    IC_C_Cor_1[0] = -1.0e6
    IC_C_Cor_2[0] = -1.0e6
    IC_C_Cor_3[0] = -1.0e6
    IC_C_Cor_4[0] = -1.0e6
    IC_C_Cor_5[0] = -1.0e6
    IC_C_Cor_6[0] = -1.0e6
    IC_C_Cor_7[0] = -1.0e6
    IC_C_Cor_8[0] = -1.0e6
    IC_C_Cor_9[0] = -1.0e6
    IC_C_Cor_10[0] = -1.0e6
    IC_C_Cor_11[0] = -1.0e6
    IC_C_Cor_12[0] = -1.0e6
    IC_C_Cor_13[0] = -1.0e6
    IC_C_Cor_14[0] = -1.0e6
    IC_C_Cor_15[0] = -1.0e6
    IC_C_Cor_16[0] = -1.0e6
    IC_C_Cor_17[0] = -1.0e6
    IC_C_Cor_18[0] = -1.0e6
    IC_C_Cor_19[0] = -1.0e6
    IC_C_Cor_20[0] = -1.0e6
    IC_C_Cor_21[0] = -1.0e6
    IC_C_Cor_22[0] = -1.0e6
    IC_C_Cor_23[0] = -1.0e6
    IC_C_Cor_24[0] = -1.0e6
    IC_C_Cor_25[0] = -1.0e6
    IC_C_Cor_26[0] = -1.0e6
    IC_C_Cor_27[0] = -1.0e6
    IC_C_Cor_28[0] = -1.0e6
    IC_C_Cor_29[0] = -1.0e6

    if S1_Y_0 > -1.0e5:  # Check if S1_Y_0 is valid

        # Calculate the corrected IC charges using the correction factor for each pad ID and the corresponding raw IC charge:
        if IC_C_0 > 350.0:
            IC_C_Cor_0[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_0) + IC_C_0
        else:
            # IC_C_Cor_0[0] = IC_C_0
            IC_C_Cor_0[0] = 0
        if IC_C_1 > 350.0:
            IC_C_Cor_1[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_1) + IC_C_1
        else:
            # IC_C_Cor_1[0] = IC_C_1
            IC_C_Cor_1[0] = 0            
        if IC_C_2 > 350.0:
            IC_C_Cor_2[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_2) + IC_C_2
        else:
            # IC_C_Cor_2[0] = IC_C_2
            IC_C_Cor_2[0] = 0
        if IC_C_3 > 350.0:
            IC_C_Cor_3[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_3) + IC_C_3
        else:
            # IC_C_Cor_3[0] = IC_C_3
            IC_C_Cor_3[0] = 0
        if IC_C_4 > 350.0:
            IC_C_Cor_4[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_4) + IC_C_4
        else:
            # IC_C_Cor_4[0] = IC_C_4
            IC_C_Cor_4[0] = 0
        if IC_C_5 > 350.0:
            IC_C_Cor_5[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_5) + IC_C_5
        else:
            # IC_C_Cor_5[0] = IC_C_5
            IC_C_Cor_5[0] = 0
        if IC_C_6 > 350.0:
            IC_C_Cor_6[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_6) + IC_C_6
        else:
            # IC_C_Cor_6[0] = IC_C_6
            IC_C_Cor_6[0] = 0
        if IC_C_7 > 350.0:
            IC_C_Cor_7[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_7) + IC_C_7
        else:
            # IC_C_Cor_7[0] = IC_C_7
            IC_C_Cor_7[0] = 0
        if IC_C_8 > 350.0:
            IC_C_Cor_8[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_8) + IC_C_8
        else:
            # IC_C_Cor_8[0] = IC_C_8
            IC_C_Cor_8[0] = 0
        if IC_C_9 > 350.0:
            IC_C_Cor_9[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_9) + IC_C_9
        else:
            # IC_C_Cor_9[0] = IC_C_9
            IC_C_Cor_9[0] = 0
        if IC_C_10 > 350.0:
            IC_C_Cor_10[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_10) + IC_C_10
        else:
            # IC_C_Cor_10[0] = IC_C_10
            IC_C_Cor_10[0] = 0
        if IC_C_11 > 350.0:
            IC_C_Cor_11[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_11) + IC_C_11
        else:
            # IC_C_Cor_11[0] = IC_C_11
            IC_C_Cor_11[0] = 0
        if IC_C_12 > 350.0:
            IC_C_Cor_12[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_12) + IC_C_12
        else:
            # IC_C_Cor_12[0] = IC_C_12
            IC_C_Cor_12[0] = 0
        if IC_C_13 > 350.0:
            IC_C_Cor_13[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_13) + IC_C_13
        else:
            # IC_C_Cor_13[0] = IC_C_13
            IC_C_Cor_13[0] = 0
        if IC_C_14 > 350.0:
            IC_C_Cor_14[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_14) + IC_C_14
        else:
            # IC_C_Cor_14[0] = IC_C_14
            IC_C_Cor_14[0] = 0
        if IC_C_15 > 350.0:
            IC_C_Cor_15[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_15) + IC_C_15
        else:
            # IC_C_Cor_15[0] = IC_C_15
            IC_C_Cor_15[0] = 0
        if IC_C_16 > 350.0:
            IC_C_Cor_16[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_16) + IC_C_16
        else:
            # IC_C_Cor_16[0] = IC_C_16
            IC_C_Cor_16[0] = 0
        if IC_C_17 > 350.0:
            IC_C_Cor_17[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_17) + IC_C_17
        else:
            # IC_C_Cor_17[0] = IC_C_17
            IC_C_Cor_17[0] = 0
        if IC_C_18 > 350.0:
            IC_C_Cor_18[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_18) + IC_C_18
        else:
            # IC_C_Cor_18[0] = IC_C_18
            IC_C_Cor_18[0] = 0
        if IC_C_19 > 350.0:
            IC_C_Cor_19[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_19) + IC_C_19
        else:
            # IC_C_Cor_19[0] = IC_C_19
            IC_C_Cor_19[0] = 0
        if IC_C_20 > 350.0:
            IC_C_Cor_20[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_20) + IC_C_20
        else:
            # IC_C_Cor_20[0] = IC_C_20
            IC_C_Cor_20[0] = 0
        if IC_C_21 > 350.0:
            IC_C_Cor_21[0] = padIDCharge_S1Y_correction_factor(S1_Y_0, pad_21) + IC_C_21
        else:
            # IC_C_Cor_21[0] = IC_C_21
            IC_C_Cor_21[0] = 0
            
        # From pad 22 to pad 29, the correction factor is not applied
        if IC_C_22 > 350.0:
            IC_C_Cor_22[0] = IC_C_22
        else:
            IC_C_Cor_22[0] = 0
        if IC_C_23 > 350.0:
            IC_C_Cor_23[0] = IC_C_23
        else:
            IC_C_Cor_23[0] = 0
        if IC_C_24 > 350.0:
            IC_C_Cor_24[0] = IC_C_24
        else:
            IC_C_Cor_24[0] = 0
        if IC_C_25 > 350.0:
            IC_C_Cor_25[0] = IC_C_25
        else:
            IC_C_Cor_25[0] = 0
        if IC_C_26 > 350.0:
            IC_C_Cor_26[0] = IC_C_26
        else:
            IC_C_Cor_26[0] = 0
        if IC_C_27 > 350.0:
            IC_C_Cor_27[0] = IC_C_27
        else:
            IC_C_Cor_27[0] = 0
        if IC_C_28 > 350.0:
            IC_C_Cor_28[0] = IC_C_28
        else:
            IC_C_Cor_28[0] = 0
        if IC_C_29 > 350.0:
            IC_C_Cor_29[0] = IC_C_29
        else:
            IC_C_Cor_29[0] = 0

        # Fill the 2D Histogram with the pad ID and the IC charge
        # For Raw IC Charge
        h2D_IC_Charge.Fill(0, IC_C_0)
        h2D_IC_Charge.Fill(1, IC_C_1)
        h2D_IC_Charge.Fill(2, IC_C_2)
        h2D_IC_Charge.Fill(3, IC_C_3)
        h2D_IC_Charge.Fill(4, IC_C_4)
        h2D_IC_Charge.Fill(5, IC_C_5)
        h2D_IC_Charge.Fill(6, IC_C_6)
        h2D_IC_Charge.Fill(7, IC_C_7)
        h2D_IC_Charge.Fill(8, IC_C_8)
        h2D_IC_Charge.Fill(9, IC_C_9)
        h2D_IC_Charge.Fill(10, IC_C_10)
        h2D_IC_Charge.Fill(11, IC_C_11)
        h2D_IC_Charge.Fill(12, IC_C_12)
        h2D_IC_Charge.Fill(13, IC_C_13)
        h2D_IC_Charge.Fill(14, IC_C_14)
        h2D_IC_Charge.Fill(15, IC_C_15)
        h2D_IC_Charge.Fill(16, IC_C_16)
        h2D_IC_Charge.Fill(17, IC_C_17)
        h2D_IC_Charge.Fill(18, IC_C_18)
        h2D_IC_Charge.Fill(19, IC_C_19)
        h2D_IC_Charge.Fill(20, IC_C_20)
        h2D_IC_Charge.Fill(21, IC_C_21)
        h2D_IC_Charge.Fill(22, IC_C_22)
        h2D_IC_Charge.Fill(23, IC_C_23)
        h2D_IC_Charge.Fill(24, IC_C_24)
        h2D_IC_Charge.Fill(25, IC_C_25)
        h2D_IC_Charge.Fill(26, IC_C_26)
        h2D_IC_Charge.Fill(27, IC_C_27)
        h2D_IC_Charge.Fill(28, IC_C_28)
        h2D_IC_Charge.Fill(29, IC_C_29)

        # For Corrected IC Charge
        h2D_IC_Charge_Corrected.Fill(0, IC_C_Cor_0[0])
        h2D_IC_Charge_Corrected.Fill(1, IC_C_Cor_1[0])
        h2D_IC_Charge_Corrected.Fill(2, IC_C_Cor_2[0])
        h2D_IC_Charge_Corrected.Fill(3, IC_C_Cor_3[0])
        h2D_IC_Charge_Corrected.Fill(4, IC_C_Cor_4[0])
        h2D_IC_Charge_Corrected.Fill(5, IC_C_Cor_5[0])
        h2D_IC_Charge_Corrected.Fill(6, IC_C_Cor_6[0])
        h2D_IC_Charge_Corrected.Fill(7, IC_C_Cor_7[0])
        h2D_IC_Charge_Corrected.Fill(8, IC_C_Cor_8[0])
        h2D_IC_Charge_Corrected.Fill(9, IC_C_Cor_9[0])
        h2D_IC_Charge_Corrected.Fill(10, IC_C_Cor_10[0])
        h2D_IC_Charge_Corrected.Fill(11, IC_C_Cor_11[0])
        h2D_IC_Charge_Corrected.Fill(12, IC_C_Cor_12[0])
        h2D_IC_Charge_Corrected.Fill(13, IC_C_Cor_13[0])
        h2D_IC_Charge_Corrected.Fill(14, IC_C_Cor_14[0])
        h2D_IC_Charge_Corrected.Fill(15, IC_C_Cor_15[0])
        h2D_IC_Charge_Corrected.Fill(16, IC_C_Cor_16[0])
        h2D_IC_Charge_Corrected.Fill(17, IC_C_Cor_17[0])
        h2D_IC_Charge_Corrected.Fill(18, IC_C_Cor_18[0])
        h2D_IC_Charge_Corrected.Fill(19, IC_C_Cor_19[0])
        h2D_IC_Charge_Corrected.Fill(20, IC_C_Cor_20[0])
        h2D_IC_Charge_Corrected.Fill(21, IC_C_Cor_21[0])
        h2D_IC_Charge_Corrected.Fill(22, IC_C_Cor_22[0])
        h2D_IC_Charge_Corrected.Fill(23, IC_C_Cor_23[0])
        h2D_IC_Charge_Corrected.Fill(24, IC_C_Cor_24[0])
        h2D_IC_Charge_Corrected.Fill(25, IC_C_Cor_25[0])
        h2D_IC_Charge_Corrected.Fill(26, IC_C_Cor_26[0])
        h2D_IC_Charge_Corrected.Fill(27, IC_C_Cor_27[0])
        h2D_IC_Charge_Corrected.Fill(28, IC_C_Cor_28[0])
        h2D_IC_Charge_Corrected.Fill(29, IC_C_Cor_29[0])


    # Fill the new branches
    IC_C_Cor_0_br.Fill()
    IC_C_Cor_1_br.Fill()
    IC_C_Cor_2_br.Fill()
    IC_C_Cor_3_br.Fill()
    IC_C_Cor_4_br.Fill()
    IC_C_Cor_5_br.Fill()
    IC_C_Cor_6_br.Fill()
    IC_C_Cor_7_br.Fill()
    IC_C_Cor_8_br.Fill()
    IC_C_Cor_9_br.Fill()
    IC_C_Cor_10_br.Fill()
    IC_C_Cor_11_br.Fill()
    IC_C_Cor_12_br.Fill()
    IC_C_Cor_13_br.Fill()
    IC_C_Cor_14_br.Fill()
    IC_C_Cor_15_br.Fill()
    IC_C_Cor_16_br.Fill()
    IC_C_Cor_17_br.Fill()
    IC_C_Cor_18_br.Fill()
    IC_C_Cor_19_br.Fill()
    IC_C_Cor_20_br.Fill()
    IC_C_Cor_21_br.Fill()
    IC_C_Cor_22_br.Fill()
    IC_C_Cor_23_br.Fill()
    IC_C_Cor_24_br.Fill()
    IC_C_Cor_25_br.Fill()
    IC_C_Cor_26_br.Fill()
    IC_C_Cor_27_br.Fill()
    IC_C_Cor_28_br.Fill()
    IC_C_Cor_29_br.Fill()


# Draw the histograms
c1.cd(1)
h2D_IC_Charge.Draw("colz")

c1.cd(2)
h2D_IC_Charge_Corrected.Draw("colz")

c1.Update()
c1.Draw()

# input("Press Enter to Save File and exit...")

# Save the Histograms to the root file
histFileName = f"/u/ddas/Lustre/gamma/ddas/SHARAQ12/output/Analysis/physOutput/24sharaq12phys_{RUN_Nbr}histFile100725.root" 
histFile = ROOT.TFile(histFileName, "RECREATE")
c1.Write()
hlist.Write("HistogramLists", ROOT.TObject.kSingleKey)
histFile.Close()


# Save the values to the root file
rootFile.cd()
tree.Write("", ROOT.TObject.kOverwrite)
rootFile.Close()