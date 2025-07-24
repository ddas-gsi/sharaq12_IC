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

# Define the constants

c = 300000000;	# velocity in m/s
amu_E = 931.494; # in MeV/u

distance_F3_FE9 = 68.54307;  # distance in meter
distance_FE9_FE12 = 14.87168; # distance in meter   
distance_FE12_S1 = 9.4892;	 # distance in meter


def calculateBeta(TOF, pathLength):
    if TOF>0:
        velocity = (pathLength / TOF) * pow(10, 9);
        beta = velocity / c;
    elif TOF <= 0:
        beta = -1000000;    # maybe put -1000000 later
    return beta;

def fn_peak38(x):
    y = np.full_like(x, 38)  # Create an array of 38 with the same shape as x
    return y

def peakS1Y_correction_function(x, *par):
    y = par[0] + par[1]*x + par[2]*x**2 + par[3]*x**3
    return y

def peakS1Y_correction_factor(x, *par):
    dy = fn_peak38(x) - peakS1Y_correction_function(x, *par)
    return dy

def peakS1X_correction_function(x, *par):
    y = par[0] + par[1]*x + par[2]*x**2 + par[3]*x**3 + par[4]*x**4 + par[5]*x**5 + par[6]*x**6 + par[7]*x**7 + par[8]*x**8
    return y

def peakS1X_correction_factor(x, *par):
    dy = fn_peak38(x) - peakS1X_correction_function(x, *par)
    return dy

def fn_range147(x):
    y = np.full_like(x, 147.6647)   # Create an array of 147.6647 with the same shape as x
    return y

def fn_rangeBetaS1LISE(x):
    betaS1 = [0.17, 0.18, 0.19]
    # range_betaS1_2_100 = [117.870, 155.759, 185.228]
    range_betaS1_2_100 = [122.23598, 155.759, 188.72326]

    # Fit a straight line (degree 1 polynomial)
    m, c = np.polyfit(betaS1, range_betaS1_2_100, 1)

    y = m * x + c
    return y


def rangeBetaS1_correction_function(x, *par):
    y = par[0] + par[1]*x
    return y

def rangeBetaS1_correction_factor(x, *par):
    # dy = fn_range147(x) - rangeBetaS1_correction_function(x, *par)
    dy = fn_rangeBetaS1LISE(x) - rangeBetaS1_correction_function(x, *par)
    return dy


def calculateBetaFromEnergy(energy):
    '''
    totalE = energy + amu_E    in MeV/u  
    '''
    amu_E = 931.494
    totalE = energy + amu_E
    gamma = totalE / amu_E
    beta = np.sqrt(1 - ((1 / gamma)**2))
    return beta

def calculateNormalizedRange(energy):
    
    if energy == 15.5:
        Range = 19 * 25.233 
    elif energy == 15.0:
        Range = 17.5 * 25.233
    elif energy == 14.5:
        Range = 16 * 25.233
    elif energy == 14.0:
        Range = 14 * 25.233
    
    beta = calculateBetaFromEnergy(energy)
    normalizedR = (Range/beta**2)/100
    print(f"Beta: {beta}, Normalized Range: {normalizedR}")
    return beta, normalizedR

# calculateNormalizedRange(14)
# plt.scatter(calculateNormalizedRange(14.0)[0], calculateNormalizedRange(14.0)[1])
# plt.scatter(calculateNormalizedRange(14.5)[0], calculateNormalizedRange(14.5)[1])
# plt.scatter(calculateNormalizedRange(15.0)[0], calculateNormalizedRange(15.0)[1])
# plt.scatter(calculateNormalizedRange(15.5)[0], calculateNormalizedRange(15.5)[1])
# paramsBetaS1 = (-694.053, 4648.12)
# betaS1 = np.linspace(0.16, 0.22, 10000)
# range_Beta2_100 = rangeBetaS1_correction_function(betaS1, *paramsBetaS1)
# plt.plot(betaS1, range_Beta2_100, label="range/betaS1^2/100 vs betaS1")
# # plt.plot(betaS1, fn_range147(betaS1), label="y = 147")
# plt.plot(betaS1, fn_rangeBetaS1LISE(betaS1), label="Range-BetaS1 fit from LISE++")
# plt.plot(betaS1, rangeBetaS1_correction_factor(betaS1, *paramsBetaS1), label="Correction Factor")
# plt.scatter(0.17, 117.870)
# plt.scatter(0.18, 155.759)
# plt.scatter(0.19, 185.228)



# paramsS1X = (38.2941, -0.0316204, 0.000472052, 1.61642e-05, -2.66285e-07, -1.73942e-09, 3.16725e-11, 5.40495e-14, -1.08301e-15)
# # paramsS1Y = (36.6426, -0.0471184, -0.000388751, -9.97883e-07)
# x = np.linspace(-200, 200, 10000)
# y = peakS1Y_correction_function(x, *paramsS1Y)
# y = peakS1X_correction_function(x, *paramsS1X)
# plt.plot(x, y, label='peak @S1X correction function')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.plot(x, fn_peak38(x), label='y = 38')
# plt.plot(x, peakS1X_correction_factor(x, *paramsS1X), label="Correction Factor")
# plt.title('Peak @S1X correction')
# plt.legend()
# plt.grid()
# plt.show()

# Read the root file
# rootFile = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/Ca_1005_Analysis.root", "UPDATE")
rootFile = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/Sc_1005_Analysis.root", "UPDATE")
tree = rootFile.Get("tree_new")     # this method also works...
# Get the number of entries in the tree
num_entries = tree.GetEntries()
# Print the number of entries
print(f"Number of entries in the TTree: {num_entries}")

# # Create a canvas 3x3 for the plots
c1 = ROOT.TCanvas("c1", "c1", 1800, 1350)
c1.Divide(3,3)

# define the params of the polynomial
paramsS1Y = (36.7351, -0.0457725, -0.00043151, -1.53105e-06)
paramsS1X = (38.2941, -0.0316204, 0.000472052, 1.61642e-05, -2.66285e-07, -1.73942e-09, 3.16725e-11, 5.40495e-14, -1.08301e-15)
paramsBetaS1 = (-694.053, 4648.12)

# Define the histograms
h_peakS1Y = ROOT.TH2F("peak_expt", "peak_expt:S1_Y_0; S1_Y_0; peak_expt", 1000, -200, 200, 1000, 1, 70)
h_peakS1Y_corrected = ROOT.TH2F("peak_S1Y_corrected", "peakS1Y:S1_Y_0; S1_Y_0; peakS1Y", 1000, -200, 200, 1000, 1, 70)
h_peakS1Y_X = ROOT.TH2F("peakS1Y", "peakS1Y:S1_X_0; S1_X_0; peakS1Y", 1000, -200, 200, 1000, 1, 70)
h_peakS1YX_corrected = ROOT.TH2F("peak_S1YX_corrected", "peakS1YX:S1_X_0; S1_X_0; peakS1YX", 1000, -200, 200, 1000, 1, 70)

h_peakS1Y_betaS1 = ROOT.TH2F("peakS1Y_betaS1", "peakS1Y:betaS1; betaS1; peakS1Y", 1000, 0.16, 0.22, 1000, 1, 70)
h_peakS1YX_betaS1 = ROOT.TH2F("peakS1YX_betaS1", "peakS1YX:betaS1; betaS1; peakS1YX", 1000, 0.16, 0.22, 1000, 1, 70)

h_rangeBeta2S1 = ROOT.TH2F("range/betaS1^2/100", "range/betaS1^2/100:betaS1; betaS1; range/betaS1^2/100", 1000, 0.16, 0.22, 1000, 0, 250)
h_rangeBeta2S1_corrected = ROOT.TH2F("range/betaS1^2/100 corrected", "Corrected range/betaS1^2/100:betaS1; betaS1; range/betaS1^2/100", 1000, 0.16, 0.22, 1000, 0, 250)
h_rangeBetaS1_corrected = ROOT.TH2F("rangeb corrected", "Corrected rangeb:betaS1; betaS1; rangeb", 1000, 0.16, 0.22, 1000, 0, 800)

# Create a TList
hlist = ROOT.TList()
# Add all the histogram to the histogram TList
hlist.Add(h_peakS1Y)
hlist.Add(h_peakS1Y_corrected)
hlist.Add(h_peakS1Y_X)
hlist.Add(h_peakS1YX_corrected)
hlist.Add(h_peakS1Y_betaS1)
hlist.Add(h_peakS1YX_betaS1)
hlist.Add(h_rangeBeta2S1)
hlist.Add(h_rangeBeta2S1_corrected)
hlist.Add(h_rangeBetaS1_corrected)


# # Create a new ROOT file and clone the structure and data of the old tree
# output_file = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/1005_Analysis.root", "RECREATE")
# output_tree = tree.CloneTree(-1)  # Clone structure and all data from the input tree

# Define variables for the new branches
peakS1Y = np.zeros(1, dtype=np.float64)     # peakS1Y is the corrected peak_expt for S1_Y_0
peakS1YX = np.zeros(1, dtype=np.float64)    # peakS1YX is the corrected peakS1Y for S1_X_0
betaFE9 = np.zeros(1, dtype=np.float64)
betaFE12 = np.zeros(1, dtype=np.float64)
betaS1 = np.zeros(1, dtype=np.float64)
rangeb = np.zeros(1, dtype=np.float64)

# Add the new branches to the output tree
# output_tree.Branch("peakS1Y", peakS1Y, "peakS1Y/D")  # /D for double.. peakS1Y is corrected the peak_expt for S1_Y_0
peakS1Y_branch = tree.Branch("peakS1Y", peakS1Y, "peakS1Y/D")  # /D for double.. peakS1Y is the corrected peak_expt for S1_Y_0
peakS1YX_branch = tree.Branch("peakS1YX", peakS1YX, "peakS1YX/D")  # /D for double.. peakS1YX is the corrected peakS1Y for S1_X_0
betaFE9_branch = tree.Branch("betaFE9", betaFE9, "betaFE9/D")  # /D for double..
betaFE12_branch = tree.Branch("betaFE12", betaFE12, "betaFE12/D")  # /D for double..
betaS1_branch = tree.Branch("betaS1", betaS1, "betaS1/D")  # /D for double..
rangeb_branch = tree.Branch("rangeb", rangeb, "rangeb/D")  # /D for double..

# Enable all branches
tree.SetBranchStatus("*", 1)    

count = 0

for i_Event in range(num_entries):        
    tree.GetEntry(i_Event)

    if (i_Event % 100000 == 0):
        print(f"Event Number: {i_Event}")


    # Define the branches
    peak_expt = getattr(tree, "peak_expt")
    S1_Y_0 = getattr(tree, "S1_Y_0")
    S1_X_0 = getattr(tree, "S1_X_0")
    PID_T_0 = getattr(tree, "PID_T_0")
    Beam_T_0 = getattr(tree, "Beam_T_0")
    S1PID_T_0 = getattr(tree, "S1PID_T_0")
    range_simFine = getattr(tree, "range_simFine")

    # Initialize new branches to -1000000
    peakS1Y[0] = -1000000
    peakS1YX[0] = -1000000
    betaFE9[0] = -1000000
    betaFE12[0] = -1000000
    betaS1[0] = -1000000
    rangeb[0] = -1000000

    # Calculate beta for different part of beamline
    betaFE9[0] = calculateBeta(PID_T_0, distance_F3_FE9)
    betaFE12[0] = calculateBeta(Beam_T_0 + 22.2, distance_FE9_FE12)      # 22.2 ns is the TOF offset for Beam_T_0
    betaS1[0] = calculateBeta(S1PID_T_0 - 467.679, distance_FE12_S1)     # 467.679 ns is the TOF offset for S1PID_T_0

    if (peak_expt>0):   

        peakS1Y[0] = peakS1Y_correction_factor(S1_Y_0, *paramsS1Y) + peak_expt
        h_peakS1Y_corrected.Fill(S1_Y_0, peakS1Y[0])

        peakS1YX[0] = peakS1X_correction_factor(S1_X_0, *paramsS1X) + peakS1Y[0]
        h_peakS1YX_corrected.Fill(S1_X_0, peakS1YX[0])

        h_rangeBeta2S1.Fill(betaS1[0], range_simFine/(betaS1[0]**2)/100)

        range_beta2_100 = rangeBetaS1_correction_factor(betaS1[0], *paramsBetaS1) + range_simFine/(betaS1[0]**2)/100
        rangeb[0] = range_beta2_100 * (betaS1[0]**2) * 100
        h_rangeBeta2S1_corrected.Fill(betaS1[0], rangeb[0]/(betaS1[0]**2)/100)
        h_rangeBetaS1_corrected.Fill(betaS1[0], rangeb[0])

        # print(f"EventNbr: {i_Event}, peak_expt: {peak_expt}, peak_expt_Y_corrected: {peakS1Y[0]}, peak_expt_YX_corrected: {peakS1YX[0]}")

        h_peakS1Y.Fill(S1_Y_0, peak_expt)
        h_peakS1Y_X.Fill(S1_X_0, peakS1Y[0])
        h_peakS1Y_betaS1.Fill(betaS1[0], peakS1Y[0])
        h_peakS1YX_betaS1.Fill(betaS1[0], peakS1YX[0])

        count = count + 1

    # output_tree.Fill()
    # tree.Fill()
    peakS1Y_branch.Fill()
    peakS1YX_branch.Fill()
    betaFE9_branch.Fill()
    betaFE12_branch.Fill()
    betaS1_branch.Fill()
    rangeb_branch.Fill()


c1.cd(1)
h_peakS1Y.Draw("colz")

c1.cd(2)
h_peakS1Y_corrected.Draw("colz")

c1.cd(3)
h_peakS1Y_X.Draw("colz")

c1.cd(4)
h_peakS1Y_betaS1.Draw("colz")

c1.cd(5)
h_peakS1YX_corrected.Draw("colz")

c1.cd(6)
h_peakS1YX_betaS1.Draw("colz")

c1.cd(7)
h_rangeBeta2S1.Draw("colz")

c1.cd(8)
h_rangeBeta2S1_corrected.Draw("colz")

c1.cd(9)
h_rangeBetaS1_corrected.Draw("colz")

c1.Update()
c1.Draw()

print(f"Total good events calculated: {count}")
input("Press Enter to exit...")

# Create a ROOT File for the Canvas and the histograms
# histFile = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/Ca_1005_histFile.root", "RECREATE")
histFile = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/Sc_1005_histFile.root", "RECREATE")
c1.Write()
hlist.Write("Histogram List", ROOT.TObject.kSingleKey)
histFile.Close()

# Save it in the Analysis File
rootFile.cd()
tree.Write("", ROOT.TObject.kOverwrite)
# rootFile.Write("", ROOT.TObject.kOverwrite)
rootFile.Close()