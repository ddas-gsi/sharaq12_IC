# ==================================================
# This code is to Fit Bragg Curve simulated from
# Analytical Function to the Experimental data 
# for all the events. 
# The Bragg curve is simulated considering the Peak
# and Range of the Experimental data
# 
# 
# 
# ++++ Next tasks ++++
# - Do scipy.optimize.curve_fit for more fine tuned Calculation
# - Do brute-force optimization for all the Good Events
# - Make the brute-force optimization parallel
# - Intelligent choice of brute-force parameters
# - in the brute_force_parallel() function return some best_chi2 value
# - Create a function that chooses the last_hit_channel in the best way for LOW Energy particles
# 
# 
# ==================================================


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


# ===================================================================================
available_cores = os.cpu_count()
print(f"Number of available CPU cores: {available_cores}")
use_cores = int(available_cores*0.13)
print(f"Running calculation on {use_cores} CPU cores. To increase number of cores, increase threshold.")
# ===================================================================================


def fncsvDictWriter(filename, list_of_dicts, field_names):
    with open(filename, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=field_names)
        writer.writeheader()
        writer.writerows(list_of_dicts)


def generateICBranchNameList(ICSegmentList, branchNameStr):
    """Returns the IC Branch Names as a List of strings"""
    ICBranchNameList = []
    for segmentNbr in ICSegmentList:
        ICBranchName = branchNameStr + str(segmentNbr)
        ICBranchNameList.append(ICBranchName)
        # print(ICBranchName)
    return ICBranchNameList

# Define parametric Bragg curve
def bragg_curve(x, peak, range_, c):
    """
    Parametric Bragg curve function.
    
    Parameters:
        x (array): x positions (material thickness or segments).
        peak (float): Maximum energy loss at the Bragg peak.
        range_ (float): Stopping range of the particle.
        c (float): Shape parameter controlling the logarithmic rise.
        
    Returns:
        dE_dx (array): Energy loss at different x positions.
    """
    Z = 20  # Atomic number of Calcium ion
    beta = 0.5  # Velocity parameter (adjustable if needed)
    normalized_x = x / range_  # Normalize x by range
    dE_dx = (
        peak * (Z**2 / beta**2)  # Base Bethe-Bloch rise
        * np.exp(-normalized_x)  # Exponential decay
        * (1 + c * np.log(1 + normalized_x))  # Logarithmic rise
    )
    # Ensure dE/dx is zero beyond range (optional for realism)
    dE_dx[normalized_x > 1] = 0.0
    return dE_dx

def fnExpLorenz(X0, X1, Xpeak, Aexp, Bexp, ALorentz, fwhm):

    x = np.arange(X0,X1,1)
    xExp = x[:Xpeak]
    xLorentz = x[Xpeak:]
    # xexp = np.linspace(X0, Xpeak, 100)
    # xLorentz = np.linspace(Xpeak, X1, 100)
    yExp = (Aexp + np.exp(Bexp * xExp)) 
    yLorentz = ALorentz/np.pi * ((fwhm/2) / ((xLorentz - Xpeak)**2 + (fwhm/2)**2))
    y = np.concatenate([yExp,yLorentz])
    # x = np.concatenate([xExp,xLorentz])
    return x,y

def fnExpGaussian(X0, X1, Xpeak, Aexp, Bexp, sigmaGauss, ampGauss):
    x = np.arange(X0,X1,1)
    xExp = x[:Xpeak]
    xGauss = x[Xpeak:]
    muGauss = Xpeak
    yExp = (Aexp + np.exp(Bexp * xExp)) 
    # yGauss = ampGauss * (1 / (sigmaGauss * np.sqrt(2 * np.pi)) * np.exp(-((xGauss - muGauss)**2 / (2 * sigmaGauss**2))))
    # yGauss = ampGauss * (1 / (np.sqrt(2 * np.pi)) * np.exp(-((xGauss - muGauss)**2 / (2 * sigmaGauss**2))))
    yGauss = ampGauss * np.exp(-((xGauss - muGauss)**2 / (2 * sigmaGauss**2)))
    y = np.concatenate([yExp, yGauss])
    return x, y

def fnExpGaussianSplined(X0, X1, NumOfPts, Xpeak, Aexp, Bexp, sigmaGauss, ampGauss):
    x = np.linspace(X0,X1,NumOfPts)
    index = np.where(x>=Xpeak)[0]
    index = index[0]
    xExp = x[:index]
    xGauss = x[index:]
    muGauss = Xpeak
    yExp = (Aexp + np.exp(Bexp * xExp)) 
    yGauss = ampGauss * np.exp(-((xGauss - muGauss)**2 / (2 * sigmaGauss**2)))
    y = np.concatenate([yExp, yGauss])
    return x, y

def fnExpGaussianFineCalculation(xArray, Xpeak, Aexp, Bexp, sigmaGauss, ampGauss):
    """
    This function is to be used by scipy.optimize.curve_fit()
    """
    index = np.where(xArray>=Xpeak)[0]
    index = index[0]
    xExp = xArray[:index]
    xGauss = xArray[index:]
    muGauss = Xpeak
    yExp = (Aexp + np.exp(Bexp * xExp)) 
    yGauss = ampGauss * np.exp(-((xGauss - muGauss)**2 / (2 * sigmaGauss**2)))
    y = np.concatenate([yExp, yGauss])
    return y

def fnCalculateRange(ElossArray):
    """
    Returns the range in mm
    """
    icStripWidth = 25.233    # in mm
    icStripTotalWidth = 756.990    # in mm
    index = np.where(ElossArray<0.1)[0]    # used threshold 0.1 from where we consider the last hit happened
    # index = np.where(ElossArray<0.1)[0] + 1
    if index.size == 0:
        # index = len(ElossArray)-1
        index = len(ElossArray)
    else:
        index = index[0]
    # print("index", index)
    dx = icStripTotalWidth / len(ElossArray)
    Range = dx * index
    return Range

def fnCalculateTotalEloss(ElossArray):
    """
    Returns total Energy loss in Ion-Chamber in MeV
    """
    index = np.where(ElossArray<0.1)[0]    # can change the threshold. Here I put 0.5
    if index.size == 0:
        index = len(ElossArray)-1
    else:
        index = index[0]
    # print("index", index)
    newElossArray = ElossArray[:index]
    totalEloss = np.sum(newElossArray)
    return totalEloss

def chi_square(params, X0, X1, expData):
    Xpeak, Aexp, Bexp, sigmaGauss, ampGauss = params
    x_, model_y = fnExpGaussian(X0, X1, Xpeak, Aexp, Bexp, sigmaGauss, ampGauss)
    chi2 = np.sum(((expData - model_y)) ** 2)
    return chi2

def parallel_bruteForce(params, X0, X1, expData):
    """
    Worker function to calculate chi-square for a given set of parameters.
    """
    return chi_square(params, X0, X1, expData)

def bruteForce_minmization(dE_dX_ExptArr_iEvent):
    # Parameter ranges for brute force
    Xpeak_range = np.arange(5,25,1)
    Aexp_range = np.linspace(7,30,25)
    Bexp_range = np.linspace(0.1,0.7,7)
    sigmaGauss_range = np.linspace(1,1.5,5)
    ampGauss_range = np.linspace(30,45,15)
    # Perform Brute-force search
    best_params = None
    best_chi2 = float('inf')
    # count = 0
    for params in product(Xpeak_range,Aexp_range,Bexp_range,sigmaGauss_range,ampGauss_range):
        chi2 = chi_square(params, 0, 30, dE_dX_ExptArr_iEvent)
        # count = count+1
        if chi2 < best_chi2:
            best_chi2 = chi2
            # print(best_chi2)
            best_params = params

    return best_params

def bruteForce_minimization_parallel(dE_dX_ExptArr_iEvent):
    """
    Perform brute-force optimization using parallelization.

    Parameters:
        dE_dX_ExptArr_iEvent (array): Experimental energy loss data.

    Returns:
        best_params (tuple): Best-fit parameters (Xpeak, Aexp, Bexp, sigmaGauss, ampGauss).
    """
    # # Parameter ranges
    # Xpeak_range = np.arange(5, 25, 1)
    # Aexp_range = np.linspace(7, 30, 25)
    # Bexp_range = np.linspace(0.1, 0.7, 7)
    # sigmaGauss_range = np.linspace(1, 1.5, 5)
    # ampGauss_range = np.linspace(30, 45, 15)

    #  Now it takes param_ranges generated in main function()

    # Generate all combinations of parameters
    param_combinations = list(
        product(Xpeak_range, Aexp_range, Bexp_range, sigmaGauss_range, ampGauss_range)
    )

    # Create a partial function to pass additional arguments to parallel_bruteForce
    worker_function = partial(parallel_bruteForce, X0=0, X1=30, expData=dE_dX_ExptArr_iEvent)

    # Use multiprocessing pool to compute in parallel. Define Pool(processes=20) to use 20 CPU cores and so on
    with Pool(processes=use_cores) as pool:
        chi2_values = pool.map(worker_function, param_combinations)

    # Find the best parameters
    best_idx = np.argmin(chi2_values)
    best_params = param_combinations[best_idx]
    return best_params

def moving_average(y, window_size):
    """
    Compute the moving average of the input data y using a specified window size.
    
    Parameters:
    y (numpy array): The data for which the moving average is to be calculated.
    window_size (int): The size of the moving window.
    
    Returns:
    numpy array: The moving average of the input data.
    """
    return np.convolve(y, np.ones(window_size)/window_size, mode='valid')

def fn_baseline_correction(y, window_size=2, threshold=0.2):
    dy = np.abs(np.diff(y))
    # print("Differences between y values (dy):", dy)
    max_diff = np.max(dy)
    # print("Max difference:", max_diff)
    max_diff_index = np.where(dy == max_diff)[0]
    # print("Max difference index:", max_diff_index)
    # print(y[max_diff_index[0]:])
    y_tail = y[max_diff_index[0]:]
    # Compute the moving average for the tail part of the data
    smoothed_y_tail = moving_average(y_tail, window_size)
    # print("Smoothed y tail:", smoothed_y_tail)
    # Compute the differences between adjacent values in the smoothed tail data
    diff = np.abs(np.diff(smoothed_y_tail))
    # print("Differences in smoothed y tail:", diff)
    if max_diff_index[0] > 25:
        threshold = 5.0
    else:
        threshold = threshold
    # Find regions where the difference is smaller than the threshold
    flat_indices = np.where(diff < threshold)[0]
    # print("Flat indices shape: ", flat_indices.shape)
    # for safety check, if flat_indices is empty, set it to 0
    if flat_indices.size == 0:
        flat_indices = np.array([0])
    else:
        pass 
    # print("Flat indices:", flat_indices)
    baseline_corrected_y = y
    baseline_corrected_y[max_diff_index[0]+flat_indices[0]:]= 0.0
    # print("Baseline corrected y:", baseline_corrected_y)
    return baseline_corrected_y


# ===================================================================================
# 
# ===================================================================================

#  ==================================================================================
# New IC Calibration Constants

# From charge to voltage:    y = A*x + B

IC_V_A = [0.000168032, 0.000155687, 0.000148696, 0.000148154, 0.000149055, 0.000156954,
                     0.000162803, 0.000174383, 0.000156753, 0.000152314, 0.000163803, 0.000156695,
                     0.000171786, 0.000165882, 0.000163502, 0.000157136, 0.000159309, 0.000151592,
                     0.000173288, 0.000160501, 0.00015524, 0.000164984, 0.000158545, 0.000144887,
                     0.000162171, 0.00016246, 0.000177322, 0.00015614, 0.000162964, 0.000161043]

IC_V_B = [0.040754, 0.0463789, 0.0584141, 0.0502596, 0.0565706, 0.054101, 0.0711771,
                     0.054378, 0.0681855, 0.049819, 0.109847, 0.00333231, 0.116113, 0.00149862,
                     0.0468207, 0.0291755, 0.0266396, 0.0293459, 0.028094, 0.0401911, 0.0266147,
                     0.0313477, 0.0269507, 0.0282567, 0.029073, 0.0322126, 0.0310689, 0.0282302,
                     0.0214893, 0.0223735]

# From voltage to energy:   y = A*x + B

IC_E_A = [185.76, 211.80, 204.22, 208.80, 209.99, 203.18, 193.55, 193.99, 193.04,
                     207.48, 163.56, 195.58, 158.65, 177.20, 189.61, 205.84, 214.86, 213.22,
                     196.01, 198.81, 216.87, 199.13, 213.05, 223.90, 213.76, 204.70, 188.67,
                     212.10, 211.80, 188.89]

IC_E_B = [-7.57, -9.82, -11.93, -10.49, -11.88, -10.99, -13.78, -10.55, -13.16,
                     -10.34, -17.97, -0.65, -18.42, -0.27, -8.88, -6.01, -5.72, -6.26, -5.51,
                     -7.99, -5.77, -6.24, -5.74, -6.33, -6.21, -6.59, -5.86, -5.99, -4.55, -4.23]

# ==================================================================================

script_start_time = datetime.now()

# Read the root file

RUN_Nbr = 1042

# rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root")
# rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/sharaq12phys_1005new.root", "READ")
# rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_1053new.root", "READ")
# rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_1056new.root", "READ")
inputFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_" + str(RUN_Nbr) + "new.root"
rootFile = ROOT.TFile.Open(inputFileName, "READ")
tree = rootFile.Get("tree_new")     # this method also works...
# Get the number of entries in the tree
num_entries = tree.GetEntries()
# num_entries = 10000
# Print the number of entries
print(f"Number of entries in the TTree: {num_entries}")

# # Open the ROOT file containing the graphical cut
# graphFileFE9_50Ca20 = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_50Ca20.root", "READ")
# FE9cut50Ca20 = graphFileFE9_50Ca20.Get("FE9pidCut_50Ca_1053_2024")
# graphFileS1_50Ca20 = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_50Ca20.root", "READ")
# S1cut50Ca20 = graphFileS1_50Ca20.Get("s1pidCut_50Ca20_1053_2024")
# graphFileFE9_51Sc21 = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/1005_FE9_PID_51Sc21_0502.root", "READ")
# FE9cut51Sc21 = graphFileFE9_51Sc21.Get("FE9pidCut_51Sc_1005")
# graphFileS1_51Sc21 = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/1005_S1_PID_51Sc21_0502.root", "READ")
# S1cut51Sc21 = graphFileS1_51Sc21.Get("s1pidCut_51Sc_1005")

# Load the C++ files containing the TCutG definitions
ROOT.gROOT.ProcessLine('.L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1035_FE9_PID_50Ca20.cxx')
ROOT.gROOT.ProcessLine('.L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1035_S1_PID_50Ca20.cxx')
FE9cut50Ca20 = ROOT.gROOT.FindObject("FE9pidCut_50Ca_1035_2024")
S1cut50Ca20 = ROOT.gROOT.FindObject("s1pidCut_50Ca20_1035_2024")


if not FE9cut50Ca20:
    print("Error: Could not load the FE9cut50Ca20 graphical cut!")
    exit()
if not S1cut50Ca20:
    print("Error: Could not load the S1cut50Ca20 graphical cut!")
    exit()

# if not FE9cut51Sc21:
#     print("Error: Could not load the FE9cut51Sc21 graphical cut!")
#     exit()
# if not S1cut51Sc21:
#     print("Error: Could not load the S1cut51Sc21 graphical cut!")
#     exit()


# Define the IC Branch Names
ICSegmentList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
IC_Calib_BranchNameList = generateICBranchNameList(ICSegmentList, "IC_E_Cal_")
IC_C_BranchNameList = generateICBranchNameList(ICSegmentList, "IC_C_")
IC_E_BranchNameList = generateICBranchNameList(ICSegmentList, "IC_E_")

# Create a list for all the good events, with calculated Range and Parameters
EventValues = []

print(f"Copying the TTree structure and {num_entries} Entries into the output TTree ...")
# Create a new ROOT file and clone the structure and data of the old tree
# output_file = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/sharaq12phys_1005new.root", "RECREATE")
# output_file = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/1005_Analysis.root", "RECREATE")
# output_file = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/Sc_1005_Analysis.root", "RECREATE")
# output_file = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/all_1053_Analysis_2024.root", "RECREATE")
# output_file = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/Ca_1053_Analysis_2024.root", "RECREATE")
# output_file = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/1053_Analysis_2024.root", "RECREATE")
# output_file = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/1056_Analysis_2024.root", "RECREATE")
outputFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/" + str(RUN_Nbr) + "_Analysis_2024.root"
output_file = ROOT.TFile(outputFileName, "RECREATE")
output_tree = tree.CloneTree(-1)  # Clone structure and all data from the input tree
# output_tree = tree.CloneTree(num_entries)  # Clone structure and 100000 data from the input tree

# Define variables for the new branches
range_expt = np.zeros(1, dtype=np.float64)
range_simCoarse = np.zeros(1, dtype=np.float64)
peak_sim = np.zeros(1, dtype=np.float64)
rangeAll = np.zeros(1, dtype=np.float64)
peakAll = np.zeros(1, dtype=np.float64)
rangeCa = np.zeros(1, dtype=np.float64)
peakCa = np.zeros(1, dtype=np.float64)

deltaE = np.zeros(1, dtype=np.float64)
totalE = np.zeros(1, dtype=np.float64)

# Add the new branches to the output tree
range_expt_branch = output_tree.Branch("range_expt", range_expt, "range_expt/D")  # /D for double
range_simCoarse_branch = output_tree.Branch("range_simCoarse", range_simCoarse, "range_simCoarse/D")  # /D for double
peak_sim_branch = output_tree.Branch("peak_sim", peak_sim, "peak_sim/D")  # /D for double
rangeAll_branch = output_tree.Branch("rangeAll", rangeAll, "rangeAll/D")  # /D for double
peakAll_branch = output_tree.Branch("peakAll", peakAll, "peakAll/D")  # /D for double
rangeCa_branch = output_tree.Branch("rangeCa", rangeCa, "rangeCa/D")  # /D for double
peakCa_branch = output_tree.Branch("peakCa", peakCa, "peakCa/D")  # /D for double

deltaE_branch = output_tree.Branch("deltaE", deltaE, "deltaE/D")  # /D for double
totalE_branch = output_tree.Branch("totalE", totalE, "totalE/D")  # /D for double

# Create a canvas 3x3 for the plots
c1 = ROOT.TCanvas("c1", "c1", 2000, 1400)
c1.Divide(4,4)

h_FE9_PID = ROOT.TH2F("FE9_PID", "FE9_X_0:PID_T_0; TOF (ns); X Position", 1000, 1050, 1150, 1000, -50, 50)
h_S1_PID = ROOT.TH2F("S1_PID", "S1_X_0:S1PID_T_0; TOF (ns); X Position", 1000, 600, 700, 1000, -200, 200)

h_FE9_PID_Ca = ROOT.TH2F("FE9_PID_Ca", "FE9_X_0:PID_T_0; TOF (ns); X Position", 1000, 1050, 1150, 1000, -50, 50)
h_S1_PID_Ca = ROOT.TH2F("S1_PID_Ca", "S1_X_0:S1PID_T_0; TOF (ns); X Position", 1000, 600, 700, 1000, -200, 200)

h_S1_E_0 = ROOT.TH1F("S1_E_0", "Beam Energy at S1; Energy (MeV/u); Counts", 1000, 1, 30)
h_S1_X_0 = ROOT.TH1F("S1_X_0", "S1 X Position; X Position; Counts", 1000, -200, 200)
h_S1_Y_0 = ROOT.TH1F("S1_Y_0", "S1 Y Position; Y Position; Counts", 1000, -200, 200)
h_S1_A_0 = ROOT.TH1F("S1_A_0", "S1 horizontal Angle; Horizontal Angle; Counts", 1000, -0.5, 0.5)
h_S1_B_0 = ROOT.TH1F("S1_B_0", "S1 vertical Angle; Vertical Angle; Counts", 1000, -0.5, 0.5)

h_S1_E_0_cut = ROOT.TH1F("S1_E_0_cut", "Beam Energy at S1; Energy (MeV/u); Counts", 1000, 1, 30)
h_S1_X_0_cut = ROOT.TH1F("S1_X_0_cut", "S1 X Position; X Position; Counts", 1000, -200, 200)
h_S1_Y_0_cut = ROOT.TH1F("S1_Y_0_cut", "S1 Y Position; Y Position; Counts", 1000, -200, 200)
h_S1_A_0_cut = ROOT.TH1F("S1_A_0_cut", "S1 horizontal Angle; Horizontal Angle; Counts", 1000, -0.5, 0.5)
h_S1_B_0_cut = ROOT.TH1F("S1_B_0_cut", "S1 vertical Angle; Vertical Angle; Counts", 1000, -0.5, 0.5)

h_IC_E_Cal_UpdatedC_2Dcut_all = ROOT.TH2F("IC_E_Cal_UpdatedC_2Dcut_all", "IC Calibrated Energy Data Updated Params; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60)

h_S1_E_0_cut_Ca = ROOT.TH1F("S1_E_0_cut_Ca", "Beam Energy at S1; Energy (MeV/u); Counts", 1000, 1, 30)
h_S1_X_0_cut_Ca = ROOT.TH1F("S1_X_0_cut_Ca", "S1 X Position; X Position; Counts", 1000, -200, 200)
h_S1_Y_0_cut_Ca = ROOT.TH1F("S1_Y_0_cut_Ca", "S1 Y Position; Y Position; Counts", 1000, -200, 200)
h_S1_A_0_cut_Ca = ROOT.TH1F("S1_A_0_cut_Ca", "S1 horizontal Angle; Horizontal Angle; Counts", 1000, -0.5, 0.5)
h_S1_B_0_cut_Ca = ROOT.TH1F("S1_B_0_cut_Ca", "S1 vertical Angle; Vertical Angle; Counts", 1000, -0.5, 0.5)

h_IC_E_Cal_UpdatedC_2Dcut_Ca = ROOT.TH2F("IC_E_Cal_UpdatedC_2Dcut_Ca", "IC Calibrated Energy Data Updated Params; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 60)

# Create a TList
hlist = ROOT.TList()

# Add all the histogram to the histogram TList
hlist.Add(h_FE9_PID)
hlist.Add(h_FE9_PID_Ca)
hlist.Add(h_S1_PID)
hlist.Add(h_S1_PID_Ca)
hlist.Add(FE9cut50Ca20)
hlist.Add(S1cut50Ca20)
# hlist.Add(FE9cut51Sc21)
# hlist.Add(S1cut51Sc21)
hlist.Add(h_S1_E_0)
hlist.Add(h_S1_E_0_cut)
hlist.Add(h_S1_E_0_cut_Ca)
hlist.Add(h_S1_X_0)
hlist.Add(h_S1_X_0_cut)
hlist.Add(h_S1_X_0_cut_Ca)
hlist.Add(h_S1_Y_0)
hlist.Add(h_S1_Y_0_cut)
hlist.Add(h_S1_Y_0_cut_Ca)
hlist.Add(h_S1_A_0)
hlist.Add(h_S1_A_0_cut)
hlist.Add(h_S1_A_0_cut_Ca)
hlist.Add(h_S1_B_0)
hlist.Add(h_S1_B_0_cut)
hlist.Add(h_S1_B_0_cut_Ca)
hlist.Add(h_IC_E_Cal_UpdatedC_2Dcut_all)
hlist.Add(h_IC_E_Cal_UpdatedC_2Dcut_Ca)

# Enable all branches
tree.SetBranchStatus("*", 1)    

countGoodEvent = 0
countCa = 0

# Read through all the Events
for i_Event in range(num_entries):        # Change it later 
    tree.GetEntry(i_Event)
 
    # Define the branches
    FE9_X_0 = getattr(tree, "FE9_X_0")
    PID_T_0 = getattr(tree, "PID_T_0")
    S1_X_0 = getattr(tree, "S1_X_0")
    S1PID_T_0 = getattr(tree, "S1PID_T_0")
    S1_E_0 = getattr(tree, "S1_E_0")
    S1_Y_0 = getattr(tree, "S1_Y_0")
    S1_A_0 = getattr(tree, "S1_A_0")
    S1_B_0 = getattr(tree, "S1_B_0")

    # Initialize new branches to -1000000.0
    range_expt[0] = -1000000.0
    range_simCoarse[0] = -1000000.0
    peak_sim[0] = -1000000.0
    rangeAll[0] = -1000000.0
    peakAll[0] = -1000000.0
    rangeCa[0] = -1000000.0
    peakCa[0] = -1000000.0
    deltaE[0] = -1000000.0
    totalE[0] = -1000000.0


    # Fill the whole histograms -------------------
    h_S1_E_0.Fill(S1_E_0)
    h_S1_X_0.Fill(S1_X_0)
    h_S1_Y_0.Fill(S1_Y_0)
    h_S1_A_0.Fill(S1_A_0)
    h_S1_B_0.Fill(S1_B_0)


    # Apply the graphical cut from FE9 and S1 focal planes. Also check if S1 Energy is greater than 0.1
    if S1_E_0 > 0.1:
    # if FE9cut50Ca20.IsInside(PID_T_0, FE9_X_0) and S1cut50Ca20.IsInside(S1PID_T_0, S1_X_0) and S1_E_0 > 0.1:
    # if FE9cut51Sc21.IsInside(PID_T_0, FE9_X_0) and S1cut51Sc21.IsInside(S1PID_T_0, S1_X_0) and S1_E_0 > 0.1:

        # Read the IC dE/dX values for the event
        IC_dE_dX_ExptArr_iEvent = []
        IC_C_ExptArr_iEvent = []
        IC_E_ExptArr_iEvent = []
        for idx in range(len(IC_C_BranchNameList)):

            IC_C_idx = getattr(tree, IC_C_BranchNameList[idx])    # DON'T comment this
            IC_E_idx = getattr(tree, IC_E_BranchNameList[idx])    # DON'T comment this
            
            IC_C_ExptArr_iEvent.append(IC_C_idx)                  # DON'T comment this
            IC_E_ExptArr_iEvent.append(IC_E_idx)                  # DON'T comment this

            IC_dE_dX_idx = getattr(tree, IC_Calib_BranchNameList[idx])
            IC_dE_dX_ExptArr_iEvent.append(IC_dE_dX_idx)

            # # Using Updated Params, do IC calibration
            # IC_E_Cal_UpdatedC_idx = (IC_C_idx * IC_V_A[idx] + IC_V_B[idx]) * IC_E_A[idx] + IC_E_B[idx]
            # IC_dE_dX_ExptArr_iEvent.append(IC_E_Cal_UpdatedC_idx)

        # added a threshold of 7 if the average of first 8 channel is more than 7
        if (np.sum(IC_dE_dX_ExptArr_iEvent[:8])/8.0) > 6:  # for 2024 data use this  

            start_time = datetime.now()

            countGoodEvent = countGoodEvent + 1
           
            IC_dE_dX_ExptArr_iEvent = np.array(IC_dE_dX_ExptArr_iEvent)
            print("RAW array: ")
            print(IC_dE_dX_ExptArr_iEvent)        

            # params for baseline correction
            window_size = 2
            threshold = 0.3  # You can adjust this to define what you consider "flat"
            # Perform Baseline correction
            IC_dE_dX_ExptArr_iEvent = fn_baseline_correction(IC_dE_dX_ExptArr_iEvent, window_size, threshold)

            print("Filtered array: ")
            print(IC_dE_dX_ExptArr_iEvent)
            print(f"Performing Brute-force optimization for {i_Event}-th event")
            
            # Define some parameters for the brute-force space
            dEdX_Max_i = np.max(IC_dE_dX_ExptArr_iEvent)
            Xpeak_i = np.where(IC_dE_dX_ExptArr_iEvent==dEdX_Max_i)[0]
            Xpeak_i = Xpeak_i[0]
            Aexp_avg = np.average(IC_dE_dX_ExptArr_iEvent[:4])
            # Parameter ranges
            Xpeak_range = np.arange(Xpeak_i-2,Xpeak_i+3, 1)
            Aexp_range = np.linspace(Aexp_avg-6, Aexp_avg+4, 11)
            Bexp_range = np.linspace(0.05,0.7,20)    # np.linspace(0.1, 0.7, 7)
            sigmaGauss_range = np.linspace(0.6, 1.5, 10)  # np.linspace(0.8, 1.5, 8)
            ampGauss_range = np.linspace(dEdX_Max_i-5, dEdX_Max_i+4, 10)
        
            # Perform Brute-force search for i-th Event
            # best_params_iEvent = bruteForce_minmization(IC_dE_dX_ExptArr_iEvent)    # brute-force optimization single core
            best_params_iEvent = bruteForce_minimization_parallel(IC_dE_dX_ExptArr_iEvent)    # brute-force optimization parallel processing
            Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best = best_params_iEvent
            print(f"Best Parameters: Xpeak = {Xpeak_best:.3f}, Aexp = {Aexp_best:.3f}, Bexp = {Bexp_best:.3f}, sigmaGauss = {sigmaGauss_best:.3f}, ampGauss = {ampGauss_best:.3f}")
        
            #  :::  Simulate the fitted curve: Coarse and Splined  :::  
            x_sim, fitted_curve = fnExpGaussian(0, 30, Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best)        
            x_simSplined, fitted_curveSplined = fnExpGaussianSplined(0, 30, 1000, Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best)
        
            #  :::  Calculate Range and Total Eloss  :::
            exptRange = fnCalculateRange(np.array(IC_dE_dX_ExptArr_iEvent))
            coarseSimRange = fnCalculateRange(fitted_curve)
            fineSimRange = fnCalculateRange(fitted_curveSplined)
            exptEloss = fnCalculateTotalEloss(np.array(IC_dE_dX_ExptArr_iEvent))
            simEloss = fnCalculateTotalEloss(fitted_curve)

            eventDict = {"EventNbr":i_Event, "exptRange":exptRange, "coarseSimRange":coarseSimRange, 
                         "fineSimRange":fineSimRange, "exptEloss":exptEloss, "simEloss":simEloss}

            print(eventDict)
            EventValues.append(eventDict)
            end_time = datetime.now()
            print(f"Opt time for {i_Event}-th event: {end_time-start_time}")


            if FE9cut50Ca20.IsInside(PID_T_0, FE9_X_0) and S1cut50Ca20.IsInside(S1PID_T_0, S1_X_0) and S1_E_0 > 0.1:

                #  :::  Write the Range and Peak for Ca events :::
                rangeCa[0] = fineSimRange
                peakCa[0] = dEdX_Max_i

                #  :::  Fill the Ca cut histograms  :::
                h_FE9_PID_Ca.Fill(PID_T_0, FE9_X_0)
                h_S1_PID_Ca.Fill(S1PID_T_0, S1_X_0)
                h_S1_E_0_cut_Ca.Fill(S1_E_0)
                h_S1_X_0_cut_Ca.Fill(S1_X_0)
                h_S1_Y_0_cut_Ca.Fill(S1_Y_0)
                h_S1_A_0_cut_Ca.Fill(S1_A_0)
                h_S1_B_0_cut_Ca.Fill(S1_B_0)

                #  :::  Fill the 2D histogram for Ca IC Calibrated Energy  :::
                for idx in range(len(IC_C_BranchNameList)):
                    h_IC_E_Cal_UpdatedC_2Dcut_Ca.Fill(idx, IC_dE_dX_ExptArr_iEvent[idx])

                # count Ca events only
                countCa = countCa + 1


            #  :::  Write the Range and Peak for All events :::
            range_expt[0] = exptRange
            range_simCoarse[0] = coarseSimRange
            peak_sim[0] = ampGauss_best
            rangeAll[0] = fineSimRange
            peakAll[0] = dEdX_Max_i

            #  ::: Calculate partial and Total ELoss for all events :::
            deltaE[0] = np.sum(IC_dE_dX_ExptArr_iEvent[:5])                # partial Eloss we calculate using firts 5 strips
            totalE[0] = np.sum(IC_dE_dX_ExptArr_iEvent)                    # total Eloss we calculate using all strips            


            #  :::  Fill the All cut histograms  :::
            h_FE9_PID.Fill(PID_T_0, FE9_X_0)
            h_S1_PID.Fill(S1PID_T_0, S1_X_0)
            h_S1_E_0_cut.Fill(S1_E_0)
            h_S1_X_0_cut.Fill(S1_X_0)
            h_S1_Y_0_cut.Fill(S1_Y_0)
            h_S1_A_0_cut.Fill(S1_A_0)
            h_S1_B_0_cut.Fill(S1_B_0)

            #   :::  Fill the 2D histogram for All IC Calibrated Energy  :::
            for idx in range(len(IC_C_BranchNameList)):
                h_IC_E_Cal_UpdatedC_2Dcut_all.Fill(idx, IC_dE_dX_ExptArr_iEvent[idx])

                  
    # print(f"EventNbr: {i_Event}, S1 Energy: {S1_E_0}, range_expt: {range_expt[0]}, range_simCoarse: {range_simCoarse[0]}, range_simFine: {range_simFine[0]}, peak_expt: {peak_expt[0]}, peak_sim: {peak_sim[0]}")

    #  :::  Fill the new Branches  :::
    range_expt_branch.Fill()
    range_simCoarse_branch.Fill()
    peak_sim_branch.Fill()
    rangeAll_branch.Fill()
    peakAll_branch.Fill()
    rangeCa_branch.Fill()
    peakCa_branch.Fill()
    deltaE_branch.Fill()
    totalE_branch.Fill()


script_end_time = datetime.now()
print(f"Script Run time for {countGoodEvent} good events: {script_end_time-script_start_time}")
print(f"Total Good Events: {countGoodEvent}")
print(f"Total Ca Events: {countCa}")


# plot the histograms
c1.cd(1)
h_FE9_PID_Ca.Draw("colz")
FE9cut50Ca20.Draw("SAME")
FE9cut50Ca20.SetLineColor(ROOT.kRed)
FE9cut50Ca20.SetLineWidth(2)
# FE9cut51Sc21.Draw("SAME")
# FE9cut51Sc21.SetLineColor(ROOT.kRed)
# FE9cut51Sc21.SetLineWidth(2)

c1.cd(2)
h_S1_PID_Ca.Draw("colz")
S1cut50Ca20.Draw("SAME")
S1cut50Ca20.SetLineColor(ROOT.kRed)
S1cut50Ca20.SetLineWidth(2)
# S1cut51Sc21.Draw("SAME")
# S1cut51Sc21.SetLineColor(ROOT.kRed)
# S1cut51Sc21.SetLineWidth(2)

c1.cd(3)
h_S1_E_0.Draw()
h_S1_E_0_cut_Ca.SetLineColor(ROOT.kRed)
h_S1_E_0_cut_Ca.Draw("SAME")

c1.cd(4)
h_S1_X_0.Draw()
h_S1_X_0_cut_Ca.SetLineColor(ROOT.kRed)
h_S1_X_0_cut_Ca.Draw("SAME")

c1.cd(5)
h_S1_Y_0.Draw()
h_S1_Y_0_cut_Ca.SetLineColor(ROOT.kRed)
h_S1_Y_0_cut_Ca.Draw("SAME")

c1.cd(6)
h_S1_A_0.Draw()
h_S1_A_0_cut_Ca.SetLineColor(ROOT.kRed)
h_S1_A_0_cut_Ca.Draw("SAME")

c1.cd(7)
h_S1_B_0.Draw()
h_S1_B_0_cut_Ca.SetLineColor(ROOT.kRed)
h_S1_B_0_cut_Ca.Draw("SAME") 

c1.cd(8)
h_IC_E_Cal_UpdatedC_2Dcut_Ca.Draw("colz")

c1.cd(9)
h_FE9_PID.Draw("colz")
FE9cut50Ca20.Draw("SAME")
FE9cut50Ca20.SetLineColor(ROOT.kRed)
FE9cut50Ca20.SetLineWidth(2)

c1.cd(10)
h_S1_PID.Draw("colz")
S1cut50Ca20.Draw("SAME")
S1cut50Ca20.SetLineColor(ROOT.kRed)
S1cut50Ca20.SetLineWidth(2)

c1.cd(11)
h_S1_E_0.Draw()
h_S1_E_0_cut.SetLineColor(ROOT.kRed)
h_S1_E_0_cut.Draw("SAME")

c1.cd(12)
h_S1_X_0.Draw()
h_S1_X_0_cut.SetLineColor(ROOT.kRed)
h_S1_X_0_cut.Draw("SAME")

c1.cd(13)
h_S1_Y_0.Draw()
h_S1_Y_0_cut.SetLineColor(ROOT.kRed)
h_S1_Y_0_cut.Draw("SAME")

c1.cd(14)
h_S1_A_0.Draw()
h_S1_A_0_cut.SetLineColor(ROOT.kRed)
h_S1_A_0_cut.Draw("SAME")

c1.cd(15)
h_S1_B_0.Draw()
h_S1_B_0_cut.SetLineColor(ROOT.kRed)
h_S1_B_0_cut.Draw("SAME")

c1.cd(16)
h_IC_E_Cal_UpdatedC_2Dcut_all.Draw("colz")


c1.Update()

# Create a ROOT File for the Canvas and the histograms
# histFile = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/Ca_1053_histFile_2024.root", "RECREATE")
# histFile = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/1053_histFile_2024.root", "RECREATE")
# histFile = ROOT.TFile("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/1056_histFile_2024.root", "RECREATE")
histFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/" + str(RUN_Nbr) + "_histFile_2024.root"
histFile = ROOT.TFile(histFileName, "RECREATE")
c1.Write()
hlist.Write("Histogram List", ROOT.TObject.kSingleKey)
histFile.Close()


output_file.cd()
output_tree.Write("", ROOT.TObject.kOverwrite)        # Overwrite the existing tree
output_file.Close()                                   # Close the output ROOT file

rootFile.Close()       # Close the input ROOT file


# filename = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/allEventsRangeEloss.csv"
# fieldnames = ["EventNbr","exptRange","coarseSimRange","fineSimRange","exptEloss","simEloss"]
# fncsvDictWriter(filename,EventValues,fieldnames)