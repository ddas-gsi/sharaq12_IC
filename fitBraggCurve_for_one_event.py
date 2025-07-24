import ROOT
import numpy as np
import scipy 
from scipy.optimize import curve_fit
import json
# import uproot     # not there in lxpool python installation
# import pandas     # not there in lxpool python installation

def parametric_bethe_bloch(x, Z, a, b, c):
    """
    Parametric Bethe-Bloch formula for fitting experimental Bragg curves.

    Parameters:
    x : array-like
        Position along the path (e.g., range in cm).
    Z : float
        Charge of the particle (fitting parameter).
    a, b, c : floats
        Empirical fitting parameters to scale and shape the curve.

    Returns:
    dE_dx : array-like
        Differential energy loss (-dE/dx) at each x.
    """
    # Experimental velocity (beta = v/c), assumed constant for simplicity
    beta = 0.6  # Placeholder value; use experimental beta if given.
    gamma = 1 / np.sqrt(1 - beta**2)

    # Parametric Bragg curve
    stopping_power = (
        a * (Z**2 / beta**2)
        * np.exp(-b * x)  # Approximation of energy loss shape
        * (1 + c * np.log(x + 1))  # Empirical correction term
    )
    return stopping_power


def fit_bragg_curve(x_exp, dE_dx_exp):
    """
    Fits the experimental Bragg curve using the parametric Bethe-Bloch formula.

    Parameters:
    x_exp : array-like
        Experimental positions (e.g., range in cm).
    dE_dx_exp : array-like
        Experimental differential energy loss values (-dE/dx).

    Returns:
    popt : array
        Optimal parameters for the fit (Z, a, b, c).
    pcov : 2D array
        Covariance of the fit parameters.
    """
    # Initial guesses for Z, a, b, c
    initial_guess = [20, 10, 0.5, 0.01]

    # Fit the parametric function to experimental data
    popt, pcov = curve_fit(parametric_bethe_bloch, x_exp, dE_dx_exp, p0=initial_guess)

    return popt, pcov



# # Example experimental data (replace with actual data)
# x_exp = np.linspace(0, 30, 30)  # Range positions in cm
# dE_dx_exp = np.random.uniform(15, 35, 30)  # Replace with experimental values

# # Fit the experimental data
# popt, pcov = fit_bragg_curve(x_exp, dE_dx_exp)

# # Extract parameters
# Z_fitted, a_fitted, b_fitted, c_fitted = popt
# print(f"Fitted Parameters: Z = {Z_fitted}, a = {a_fitted}, b = {b_fitted}, c = {c_fitted}")

# --------------------------------------------------------------

# --------------------------------------------------------------

# rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root")
# tree = rootFile.Get("tree_new")     

# def generateICBranchNameList(ICSegmentList, branchNameStr):
#     """Returns the IC Branch Names as a List of strings"""
#     ICBranchNameList = []
#     for segmentNbr in ICSegmentList:
#         ICBranchName = branchNameStr + str(segmentNbr)
#         ICBranchNameList.append(ICBranchName)
#         # print(ICBranchName)
#     return ICBranchNameList

# ICSegmentList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
# ICBranchNameList = generateICBranchNameList(ICSegmentList, "IC_E_Cal_")

# i_Event = 363409
# tree.GetEntry(i_Event)
# IC_E_Cal_ExptArr_iEvent = []
# for idx in range(len(ICBranchNameList)):
#     # print(getattr(tree, ICBranchNameList[idx]))
#     ic_E_Cal_idx = getattr(tree, ICBranchNameList[idx])
#     IC_E_Cal_ExptArr_iEvent.append(ic_E_Cal_idx)

# print(IC_E_Cal_ExptArr_iEvent)



# Create a TCanvas
canvas = ROOT.TCanvas("c1", "IC Canvas", 2500, 750)
# Divide the canvas into 2 parts (1 row, 2 columns)
canvas.Divide(2, 1)  # (columns, rows) 
# Create a 2D histogram
hist2d_IC_E_Cal = ROOT.TH2F("IC_E_Cal", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 50)

IC_E_Cal_ExptArr_iEvent = [19.30883442006729, 20.07757334620357, 18.961192169495646, 20.863828861326464, 20.48170201163091, 21.36438847886429, 21.63018097893729, 21.78063574955428, 23.27558740839977, 24.935359534156735, 23.598545133143375, 26.10346947974618, 25.56415320656543, 28.99915234738814, 32.45695704909536, 36.3346556802946, 41.47629481050695, 42.45185952712907, 33.99725735284845, 8.003822782169488, 6.399985336097896, 6.164303604087812, 6.3044716518416575, 6.421929781555445, 6.211521150003354, 6.161777665149174, 6.5173214041458625, 6.061974586086006, 6.132044328979943, 5.933165787595918]

for idx in range(len(IC_E_Cal_ExptArr_iEvent)):
    hist2d_IC_E_Cal.Fill(idx, IC_E_Cal_ExptArr_iEvent[idx])

canvas.cd(1)
hist2d_IC_E_Cal.Draw("BOX")




# # Example experimental data (replace with actual data)
# x_exp = np.linspace(0, 30, 30)  # Range positions in cm
# dE_dx_exp = np.random.uniform(15, 35, 30)  # Replace with experimental values

# # Fit the experimental data
x_exp = np.linspace(0,30,30)
popt, pcov = fit_bragg_curve(x_exp, IC_E_Cal_ExptArr_iEvent)

# # Extract parameters
Z_fitted, a_fitted, b_fitted, c_fitted = popt
print(f"Fitted Parameters: Z = {Z_fitted}, a = {a_fitted}, b = {b_fitted}, c = {c_fitted}")


fitted_dE_dX = parametric_bethe_bloch(x_exp, *popt)

canvas.cd(2)
gr = ROOT.TGraph(30, x_exp, fitted_dE_dX)
hist2d_IC_E_Cal.Draw("BOX")
gr.SetMarkerStyle(20)
gr.SetMarkerColor(2)
gr.Draw('P SAME')





# Update the canvas
canvas.Update()
# Close the root file
# rootFile.Close()
# Keep the Canvas open in different ways
# Keep the canvas open
input("Press Enter to exit...")
# Keep the canvas open with the ROOT event loop
# ROOT.gApplication.Run()