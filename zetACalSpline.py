#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from itertools import product
from datetime import datetime
from multiprocessing import Pool
from functools import partial
import os
import csv
import ROOT
import sys


# Parameters for Z correction for a3EBetaS1 from spline calculation
# [mean, sigma, Z]
Ca5120_par_a3EBetaS1 = [20.3690, 0.32591, 20]
Ca5020_par_a3EBetaS1 = [20.0405, 1.06117, 20]
Ca5019_par_a3EBetaS1 = [19.4142, 0.9853, 20]
Sc5121_par_a3EBetaS1 = [21.9879, 1.09855, 21]
Sc5120_par_a3EBetaS1 = [22.3151, 1.99474, 21]

# Parameters for Z correction for sPeakES1X from spline calculation
# [mean, sigma, Z]
Ca5020_par_sPeakES1X = [37.9078, 1.2168, 20]
Ca5019_par_sPeakES1X = [38.3266, 1.35241, 20]
Sc5121_par_sPeakES1X = [39.641, 0.964208, 21]

# Parameters for A correction for sRangeBetaS1 from spline calculation
# [mean, sigma, A]
Ca5020_par_sRangeBetaS1 = [149.196, 7.73661, 50]
Ca5120_par_sRangeBetaS1 = [144.054, 9.93274, 51]
Sc5121_par_sRangeBetaS1 = [134.681, 6.3936, 51]



def Z_correction_a3EBetaS1():
    data = np.array([
    Ca5020_par_a3EBetaS1,
    Ca5120_par_a3EBetaS1,
    Sc5121_par_a3EBetaS1,
    ])  # Use only Ca5020 and Sc5121 for the fit

    X = data[:, 0]         # mean (horizontal)
    sigma_X = data[:, 1]   # uncertainty in mean
    Y = data[:, 2]         # Z (vertical)

    # Fit: Y = m*X + c
    coeffs = np.polyfit(X, Y, deg=1)
    m, c = coeffs

    # Fit line for plotting
    x_fit = np.linspace(min(X) - 1, max(X) + 1, 100)
    y_fit = m * x_fit + c

    print(f"Z correction coefficients: m = {m:.6f}, c = {c:.6f}")
    return m, c, x_fit, y_fit, X, Y, sigma_X


def Z_correction_sPeakES1X():
    data = np.array([
    Ca5020_par_sPeakES1X,
    Ca5019_par_sPeakES1X,
    Sc5121_par_sPeakES1X,
    ])  # Use only Ca5020 and Sc5121 for the fit

    X = data[:, 0]         # mean (horizontal)
    sigma_X = data[:, 1]   # uncertainty in mean
    Y = data[:, 2]         # Z (vertical)

    # Fit: Y = m*X + c
    coeffs = np.polyfit(X, Y, deg=1)
    m, c = coeffs

    # Fit line for plotting
    x_fit = np.linspace(min(X) - 1, max(X) + 1, 100)
    y_fit = m * x_fit + c

    print(f"Z correction coefficients: m = {m:.6f}, c = {c:.6f}")
    return m, c, x_fit, y_fit, X, Y, sigma_X


def A_correction_sRangeBetaS1():
    data = np.array([
    Ca5020_par_sRangeBetaS1,
    Ca5120_par_sRangeBetaS1,
    Sc5121_par_sRangeBetaS1,
    ])  # Use only Ca5020 and Sc5121 for the fit

    X = data[:, 0]         # mean (horizontal)
    sigma_X = data[:, 1]   # uncertainty in mean
    Y = data[:, 2]         # Z (vertical)

    # Fit: Y = m*X + c
    coeffs = np.polyfit(X, Y, deg=1)
    m, c = coeffs

    # Fit line for plotting
    x_fit = np.linspace(min(X) - 1, max(X) + 1, 100)
    y_fit = m * x_fit + c

    print(f"A correction coefficients: m = {m:.6f}, c = {c:.6f}")
    return m, c, x_fit, y_fit, X, Y, sigma_X


def plot_fit(m, c, x_fit, y_fit, X, Y, sigma_X):
    plt.figure(figsize=(10, 6))
    plt.errorbar(X, Y, xerr=sigma_X, fmt='o', label='Data with σ_X', capsize=5)
    plt.plot(x_fit, y_fit, 'r-', label=f'Fit: Z = {m:.3f}·mean + {c:.3f}')
    plt.xlabel('Mean')
    plt.ylabel('Z')
    plt.title('Linear Fit: Z vs Mean')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# Get the coefficients for Z correction
# m, c, x_fit, y_fit, X, Y, sigma_X = Z_correction_a3EBetaS1()

# Plotting
# plot_fit(m, c, x_fit, y_fit, X, Y, sigma_X)




def main(RUN_Nbr):

    print(f"Processing RUN Number: {RUN_Nbr} ...")
    print(f"Copying the root file {RUN_Nbr}_Spline_2024.root to the zetCorr folder ...")

    os.system(f"cp /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/AoQ/{RUN_Nbr}_Spline_2024.root \
              /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/AoQ/zetCorr/{RUN_Nbr}_Spline_2024.root")

    inputFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/splineCorrection/AoQ/zetCorr/" + str(RUN_Nbr) + "_Spline_2024.root"
    rootFile = ROOT.TFile.Open(inputFileName, "UPDATE")
    tree = rootFile.Get("tree_new")
    num_entries = tree.GetEntries()
    # num_entries = 1000  # For testing purposes, limit to 1000 entries
    print(f"Number of entries in the TTree: {num_entries}")

    # Enable all branches first
    tree.SetBranchStatus("*", 1)

    # define new variables
    Z_a3EBetaS1 = np.zeros(1, dtype=np.float64)
    Z_sPeakES1X = np.zeros(1, dtype=np.float64)
    A_sRangeBetaS1 = np.zeros(1, dtype=np.float64)

    # Create new branches for the new variables
    Z_a3EBetaS1_branch = tree.Branch("Z_a3EBetaS1", Z_a3EBetaS1, "Z_a3EBetaS1/D")
    Z_sPeakES1X_branch = tree.Branch("Z_sPeakES1X", Z_sPeakES1X, "Z_sPeakES1X/D")
    A_sRangeBetaS1_branch = tree.Branch("A_sRangeBetaS1", A_sRangeBetaS1, "A_sRangeBetaS1/D")

    countGoodEvent = 0

    print("Calculation Z & A correction coefficients ...")
    # Get the coefficients for Z and A correction
    m_a3EBetaS1, c_a3EBetaS1, _, _, _, _, _ = Z_correction_a3EBetaS1()
    m_sPeakES1X, c_sPeakES1X, _, _, _, _, _ = Z_correction_sPeakES1X()
    m_sRangeBetaS1, c_sRangeBetaS1, _, _, _, _, _ = A_correction_sRangeBetaS1()

    # Loop over the entries in the tree
    for i_Event in range(num_entries):
        tree.GetEntry(i_Event)

        if (i_Event % 100000 == 0):
            print(f"Event Number: {i_Event}")

        # Define the branches
        a3EBetaS1 = getattr(tree, "a3EBetaS1")
        sPeakES1X = getattr(tree, "sPeakES1X")
        sRangeBetaS1 = getattr(tree, "sRangeBetaS1")

        # Initialize new variables to -1e6
        Z_a3EBetaS1[0] = -1e6
        Z_sPeakES1X[0] = -1e6
        A_sRangeBetaS1[0] = -1e6

        # Run Calculation
        if a3EBetaS1 > 0:
            # Calculate Z_a3EBetaS1 using the linear fit coefficients
            Z_a3EBetaS1[0] = m_a3EBetaS1 * a3EBetaS1 + c_a3EBetaS1
            Z_sPeakES1X[0] = m_sPeakES1X * sPeakES1X + c_sPeakES1X
            A_sRangeBetaS1[0] = m_sRangeBetaS1 * sRangeBetaS1 + c_sRangeBetaS1
            
            countGoodEvent += 1



        Z_a3EBetaS1_branch.Fill()
        Z_sPeakES1X_branch.Fill()
        A_sRangeBetaS1_branch.Fill()


    # Write the new branches to the output file
    rootFile.cd()
    rootFile.Write("", ROOT.TObject.kOverwrite)
    rootFile.Close()
    print(f"Total Number of Events processed: {countGoodEvent}")
    print(f"Finished processing RUN Number: {RUN_Nbr}")


if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("Usage: python3 zetACalSpline.py RUN_Nbr")
        sys.exit(1)

    RUN_Nbr = int(sys.argv[1])


    # RUN_Nbr = 1053          # Enter RUN number
    start_time = datetime.now()
    main(RUN_Nbr)
    end_time = datetime.now()
    print(f"Time taken for processing: {end_time - start_time}")









































# ----------------------------------------------------------------------
# # Example of using Orthogonal Distance Regression (ODR) for a more robust fit
# # Your data: [mean, sigma, Z]
# data = np.array([
#     # [20.3690, 0.32591, 20],  # Ca5120
#     [20.0405, 1.06117, 20],    # Ca5020
#     # [19.4142, 0.9853, 20],   # Ca5019
#     [21.9879, 1.09855, 21],    # Sc5121
#     # [22.3151, 1.99474, 21]   # Sc5120
# ])

# mean = data[:, 0]       # X
# sigma = data[:, 1]      # σ_X
# Z = data[:, 2]          # Y

# # Define linear model: Z = m*mean + c
# def linear_func(B, x):
#     return B[0] * x + B[1]  # B[0]=m, B[1]=c

# # Prepare data for ODR
# model = Model(linear_func)
# data_odr = RealData(mean, Z, sx=sigma, sy=None)  # sx: uncertainty in X
# odr = ODR(data_odr, model, beta0=[1., 0.])       # initial guess m=1, c=0
# output = odr.run()

# # Extract results
# m, c = output.beta
# m_err, c_err = output.sd_beta

# # Generate fit line
# x_fit = np.linspace(min(mean) - 1, max(mean) + 1, 100)
# y_fit = m * x_fit + c

# # Plot
# plt.errorbar(mean, Z, xerr=sigma, fmt='o', capsize=5, label='Data with σ_X')
# plt.plot(x_fit, y_fit, 'r-', label=f'ODR Fit: Z = {m:.3f}±{m_err:.3f} · mean + {c:.3f}±{c_err:.3f}')
# plt.xlabel('Mean')
# plt.ylabel('Z')
# plt.title('Orthogonal Distance Regression: Z vs Mean')
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()