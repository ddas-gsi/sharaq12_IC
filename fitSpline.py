#!/usr/bin/env python3

from scipy.interpolate import make_interp_spline
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
from scipy.stats import binned_statistic
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, interp1d
from scipy.optimize import brute



# EnergyLoss Spline from LISE++ =======================================================
# # for 14 MeV/u 50Ca20 ion
# EnergyLoss = np.array([
#     22.136, 22.854, 23.647, 24.529, 25.517, 26.636, 27.914, 29.39, 31.114, 33.135,
#     35.539, 37.947, 34.116, 9.291, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
# ])

# # for 14.5 MeV/u 50Ca20 ion
# EnergyLoss = np.array([
#     21.084, 21.698, 22.369, 23.11, 23.932, 24.848, 25.878, 27.047, 28.387, 29.929,
#     31.75, 33.901, 36.45, 38.153, 27.575, 3.083, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
# ])

# for 15.0 MeV/u 50Ca20 ion
EnergyLoss = np.array([
    20.142, 20.676, 21.257, 21.889, 22.582, 23.346, 24.192, 25.138, 26.203, 27.416,
    28.812, 30.437, 32.351, 34.599, 37.211, 37.504, 19.631, 0.494, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
])

# # for 15.5 MeV/u 50Ca20 ion
# EnergyLoss = np.array([
#     19.307, 19.777, 20.276, 20.82, 21.413, 22.06, 22.77, 23.554, 24.425, 25.4,
#     26.502, 27.759, 29.211, 30.905, 32.892, 35.247, 37.76, 35.479, 12.255, 0, 0, 0,
#     0, 0, 0, 0, 0, 0, 0, 0
# ])

# # for 16 MeV/u 50Ca20 ion
# EnergyLoss = np.array([
#     18.562, 18.961, 19.404, 19.881, 20.395, 20.945, 21.548, 22.208, 22.933, 23.734,
#     24.626, 25.628, 26.761, 28.058, 29.558, 31.303, 33.371, 35.812, 38.078, 32.5,
#     6.919, 0, 0, 0, 0, 0, 0, 0, 0, 0
# ])    

# # for 17.0 MeV/u 50Ca20 ion
# EnergyLoss=np.array([
#     17.217, 17.573, 17.884, 18.284, 18.689, 19.102, 19.551, 20.038, 20.565, 21.136,
#     21.758, 22.433, 23.181, 24.01, 24.935, 25.977, 27.16, 28.513, 30.08, 31.931,
#     34.112, 36.688, 38.093, 25.243, 2.078, 0, 0, 0, 0, 0
# ])   

# # for 17.5 MeV/u 50Ca20 ion
# EnergyLoss= np.array([
#     16.657, 16.909, 17.254, 17.61, 17.919, 18.327, 18.729, 19.15, 19.6, 20.091,
#     20.622, 21.197, 21.825, 22.512, 23.263, 24.102, 25.038, 26.093, 27.29, 28.661,
#     30.261, 32.144, 34.359, 36.958, 37.908, 22.37, 1.168, 0, 0, 0
# ])  

# # for 17.5 MeV/u 49K19 ion
# EnergyLoss= np.array([
#     14.874, 15.134, 15.364, 15.627, 15.936, 16.221, 16.526, 16.882, 17.229, 17.605, 
#     18.002, 18.432, 18.898, 19.402, 19.945, 20.534, 21.165, 21.903, 22.699, 23.588,
#     24.591, 25.732, 27.045, 28.572, 30.369, 32.487, 35.002, 35.602, 19.412, 0.588 
# ])

# # for 17.5 MeV/u 51Sc21 ion
# EnergyLoss= np.array([
#     18.490, 18.894, 19.274, 19.680, 20.144, 20.600, 21.102, 21.644, 22.232, 22.870, 
#     23.567, 24.332, 25.176, 26.118, 27.171, 28.365, 29.731, 31.310, 33.154, 35.316,
#     37.878, 40.146, 34.696, 7.662, 0, 0, 0, 0, 0, 0
# ])
# =======================================================

# # Experimental Data =======================================================
# dE_dX_Exp = np.array([
#     23.15828857, 26.25686581, 27.12598498, 27.88286197, 29.55377179, 30.43793649,
#     32.06677529, 33.60875744, 36.65823772, 39.9846389, 36.14198983, 33.06684032,
#     10.34674719, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
# ])

# dE_dX_Exp = np.array([
#     19.30883442006729, 20.07757334620357, 18.961192169495646, 20.863828861326464,
#     20.48170201163091, 21.36438847886429, 21.63018097893729, 21.78063574955428,
#     23.27558740839977, 24.935359534156735, 23.598545133143375, 26.10346947974618,
#     25.56415320656543, 28.99915234738814, 32.45695704909536, 36.3346556802946,
#     41.47629481050695, 42.45185952712907, 33.99725735284845, 8.003822782169488,
#     6.399985336097896, 0.0, 0.0, 0.0,
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0
# ])  # Add all points

# dE_dX_Exp =  np.array([ 9.168516,   10.01462886, 10.2900509,  11.87485822, 11.7350413 , 12.04700109,
#  11.8688846,  12.40485406, 13.97172702 ,13.91271879 ,13.48389765 ,14.27730995,
#  13.60621551 ,16.68141319, 17.82849901 ,20.74381663, 21.33160314, 22.05066379,
#  24.77273985, 26.24818362 ,31.27253815 ,31.73278556 ,36.89880655 ,35.82970397,
#  18.34126153 , 0.  ,        0.     ,     0.  ,        0. ,         0.        ])

# =======================================================




# Build the spline =======================================================

ICSegment = np.arange(30)

x_spline = np.linspace(ICSegment.min(), ICSegment.max(), 300)
spline_fn = make_interp_spline(ICSegment, EnergyLoss, k=3)
y_spline = spline_fn(x_spline)


def generateICBranchNameList(ICSegmentList, branchNameStr):
    """Returns the IC Branch Names as a List of strings"""
    ICBranchNameList = []
    for segmentNbr in ICSegmentList:
        ICBranchName = branchNameStr + str(segmentNbr)
        ICBranchNameList.append(ICBranchName)
        # print(ICBranchName)
    return ICBranchNameList

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
    flat_indices_y = np.where(diff < threshold)[0]
    # print("Flat indices_y shape: ", flat_indices_y.shape)
    # for safety check, if flat_indices_y is empty, set it to 0
    if flat_indices_y.size == 0:
        flat_indices_y = np.array([0])
    else:
        pass 
    # print("Flat indices_y:", flat_indices_y)
    baseline_corrected_y = y
    baseline_corrected_y[max_diff_index[0]+flat_indices_y[0]:]= 0.0
    # print("Baseline corrected y:", baseline_corrected_y)
    return baseline_corrected_y


def fit_spline_chi2_min(filtered_dE_dX_Exp, x_spline, y_spline):
    
    ICSegment = np.arange(30)
    mask_exp = filtered_dE_dX_Exp > 0
    x_data = ICSegment[mask_exp]
    y_target = filtered_dE_dX_Exp[mask_exp]

    indices_x = [np.abs(y_spline - min(y_target)).argmin()]
    x_scaling =  x_data[-1] / x_spline[indices_x]
    indices_y = [np.abs(x_spline*x_scaling - val).argmin() for val in x_data]

    matched_y_spline = y_spline[indices_y]
    y_target = np.array(y_target)

    # matched_y_spline = y_spline[indices_y[:-1]]
    # y_target = np.array(y_target[:-1])  # Exclude the last element to match lengths

    y_scaling_range = np.linspace(0.5, 2.0, 200)

    dict_chi2 = {}

    for y_scaling in y_scaling_range:
        scaled_spline = matched_y_spline * y_scaling
        chi2 = np.sum((y_target - scaled_spline) ** 2)

        dict_chi2[y_scaling] = chi2

    # Find the optimal scaling factor with the minimum chi-square
    optimal_y_scaling = min(dict_chi2, key=dict_chi2.get)
    chi2 = dict_chi2[optimal_y_scaling]
    y_scaling = optimal_y_scaling


    print("Optimal x_scaling:", x_scaling)
    print("Optimal y_scaling:", y_scaling)
    print("Chi-square:", chi2)

    return x_scaling, y_scaling 

def fit_spline(filtered_dE_dX_Exp, x_spline, y_spline):    
    ICSegment = np.arange(30)
    mask_exp = filtered_dE_dX_Exp > 0
    x_data = ICSegment[mask_exp]
    y_target = filtered_dE_dX_Exp[mask_exp]
    
    indices_x = [np.abs(y_spline - min(y_target)).argmin()]
    x_scaling =  x_data[-1] / x_spline[indices_x]
    indices_y = [np.abs(x_spline*x_scaling - val).argmin() for val in x_data]

    matched_y_spline = y_spline[indices_y[:-1]]

    y_target = np.array(y_target[:-1])  # Exclude the last element to match lengths

    # Compute optimal scaling factor to minimize chi-square
    numerator = np.sum(y_target * matched_y_spline)
    denominator = np.sum(matched_y_spline ** 2)
    y_scaling = numerator / denominator

    # Optionally scale and compute chi2
    scaled_spline = matched_y_spline * y_scaling
    chi2 = np.sum((y_target - scaled_spline) ** 2)
    reduced_chi2 = chi2 / (len(y_target) - 2)  # Reduced chi-square = chi2 / (N - p)

    print("Optimal x_scaling:", x_scaling)
    print("Optimal y_scaling:", y_scaling)
    print("Chi-square:", chi2)
    print("Reduced Chi-square:", reduced_chi2)

    return x_scaling, y_scaling, reduced_chi2

def calculateRange(x_scaling, y_scaling, x_spline, y_spline):

    ICActiveRegion = 756.990

    scaled_x_spline = x_spline * x_scaling
    scaled_y_spline = y_spline * y_scaling
    mask_y_spline = scaled_y_spline > 1
    range_x_spline = scaled_x_spline[mask_y_spline]

    if len(range_x_spline) >= 1:
        range_value = ICActiveRegion / 30 * range_x_spline[-1]
    else:
        range_value = ICActiveRegion / 30

    print("Range: ", range_value)
    
    return range_value

def calculateTotalEnergy(filtered_dE_dX_Exp, x_scaling, y_scaling, x_spline, y_spline):
    
    ICSegment = np.arange(30)
    mask_exp = filtered_dE_dX_Exp > 0
    x_data = ICSegment[mask_exp]

    indices_y = [np.abs(x_spline*x_scaling - val).argmin() for val in x_data]
    totalEnergy = np.sum(y_spline[indices_y] * y_scaling)
    print("Total Energy:", totalEnergy)
    print("Spline Points: ",y_spline[indices_y] * y_scaling)
    
    return totalEnergy

# def calculateTotalEnergy2(x_scaling, y_scaling, x_spline, y_spline):
#     scaled_x_spline = x_spline * x_scaling
#     scaled_y_spline = y_spline * y_scaling
#     xbin_means, xbin_edges, xbinnumber = binned_statistic(scaled_x_spline, scaled_x_spline, statistic='mean', bins=30)
#     ybin_means, ybin_edges, ybinnumber = binned_statistic(scaled_x_spline, scaled_y_spline, statistic='mean', bins=30)
#     mask = ybin_means > 1
#     print("Binned Y:", ybin_means[mask])
#     print("Total Energy: ", np.sum(ybin_means[mask]))

def calculatePeaKEnergy(y_scaling, y_spline):
    
    peak_energy = np.max(y_spline) * y_scaling
    print("Peak Energy:", peak_energy)

    return peak_energy

    

# fit_spline(dE_dX_Exp, x_spline, y_spline)



def main(RUN_Nbr):

    # inputFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/" + str(RUN_Nbr) + "_Analysis_2024_new1206.root"
    
    print(f"Processing RUN Number: {RUN_Nbr} ...")
    # print(f"Copying the root file 24sharaq12phys_{RUN_Nbr}new100725.root to the outputSpline folder as {RUN_Nbr}_Spline_2024.root ...") 
    print(f"Copying the root file 24sharaq12phys_{RUN_Nbr}newSpline.root to the outputSpline folder as {RUN_Nbr}_Spline_2024.root ...")

    os.system(f"cp /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_{RUN_Nbr}newSpline.root \
              /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/{RUN_Nbr}_Spline_2024.root")

    inputFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputSpline/" + str(RUN_Nbr) + "_Spline_2024.root"
    rootFile = ROOT.TFile.Open(inputFileName, "UPDATE")
    tree = rootFile.Get("tree_new")
    # Get the number of entries in the tree
    num_entries = tree.GetEntries()
    # num_entries = 1000  # For testing purposes, limit to 1000 entries
    print(f"Number of entries in the TTree: {num_entries}")
    count = 0

    # Load the files containing the TCutG definitions
    ROOT.gROOT.ProcessLine('.L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_50Ca20.cxx')
    ROOT.gROOT.ProcessLine('.L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_51Sc21.cxx')
    ROOT.gROOT.ProcessLine('.L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_FE9_PID_49K19.cxx')

    ROOT.gROOT.ProcessLine('.L /u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/graphCutsPhys/2024_1053_S1_PID_50Ca20.cxx')    
    
    FE9cut50Ca20 = ROOT.gROOT.FindObject("FE9pidCut_50Ca_1053_2024")
    FE9cut51Sc21 = ROOT.gROOT.FindObject("FE9pidCut_51Sc_1053_2024")
    FE9cut49K19 = ROOT.gROOT.FindObject("FE9pidCut_49K_1053_2024")
    S1cut50Ca20 = ROOT.gROOT.FindObject("s1pidCut_50Ca20_1053_2024")

    # Enable all branches first
    tree.SetBranchStatus("*", 1)

    ICSegmentList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
    IC_Calib_BranchNameList = generateICBranchNameList(ICSegmentList, "IC_E_Cal_")

    # define new splineRange, splineTotalE and splinePeakE variables
    splineRange = np.zeros(1, dtype=np.float64)
    splineTotalE = np.zeros(1, dtype=np.float64)
    splinePeakE = np.zeros(1, dtype=np.float64)
    x_scaling = np.zeros(1, dtype=np.float64)
    y_scaling = np.zeros(1, dtype=np.float64)
    chi2 = np.zeros(1, dtype=np.float64)
    avg3padE = np.zeros(1, dtype=np.float64)


    # define new splineRange, splineTotalE and splinePeakE branches
    splineRange_branch = tree.Branch("splineRange", splineRange, "splineRange/D")
    splineTotalE_branch = tree.Branch("splineTotalE", splineTotalE, "splineTotalE/D")
    splinePeakE_branch = tree.Branch("splinePeakE", splinePeakE, "splinePeakE/D")
    x_scaling_branch = tree.Branch("x_scaling", x_scaling, "x_scaling/D")
    y_scaling_branch = tree.Branch("y_scaling", y_scaling, "y_scaling/D")
    chi2_branch = tree.Branch("chi2", chi2, "chi2/D")
    avg3padE_branch = tree.Branch("avg3padE", avg3padE, "avg3padE/D")


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

        # Initialize the splineRange, splineTotalE and splinePeakE values to -1e6
        splineRange[0] = -1e6
        splineTotalE[0] = -1e6
        splinePeakE[0] = -1e6
        x_scaling[0] = -1e6
        y_scaling[0] = -1e6
        chi2[0] = -1e6
        avg3padE[0] = -1e6


        if S1_E_0 > 0.1:
        # if FE9cut50Ca20.IsInside(PID_T_0, FE9_X_0) and S1cut50Ca20.IsInside(S1PID_T_0, S1_X_0) and S1_E_0 > 0.1:
            # Read the IC dE/dX values for the event
            IC_dE_dX_ExptArr_iEvent = []
            
            for idx in range(len(IC_Calib_BranchNameList)):

                IC_dE_dX_idx = getattr(tree, IC_Calib_BranchNameList[idx])
                IC_dE_dX_ExptArr_iEvent.append(IC_dE_dX_idx)

            if (np.sum(IC_dE_dX_ExptArr_iEvent[:8])/8.0) > 6 and (np.sum(IC_dE_dX_ExptArr_iEvent[:8])/8.0) < 1.0e4:  # for 2024 data use this.. (np.sum(IC_dE_dX_ExptArr_iEvent[:8])/8.0) < 1.0e4 added for the problem caused by 557039th event of 1053_Spline_2024.root file

                countGoodEvent = countGoodEvent + 1
           
                IC_dE_dX_ExptArr_iEvent = np.array(IC_dE_dX_ExptArr_iEvent)
                print("RAW array: ")
                print(IC_dE_dX_ExptArr_iEvent)        

                # ::::: Perform Baseline correction :::::

                # Params for baseline correction
                window_size = 2
                threshold = 0.3  # You can adjust this to define what you consider "flat"
                filtered_dE_dX_Exp = fn_baseline_correction(IC_dE_dX_ExptArr_iEvent, window_size, threshold)

                print("Filtered array: ")
                print(filtered_dE_dX_Exp)
                print(f"Performing Spline optimization for {i_Event}-th event")

                # Fit the spline to the filtered experimental data
                x_scaling_, y_scaling_, reduced_chi2 = fit_spline(filtered_dE_dX_Exp, x_spline, y_spline)
                sp_Range = calculateRange(x_scaling_, y_scaling_, x_spline, y_spline)
                sp_TotalE = calculateTotalEnergy(filtered_dE_dX_Exp, x_scaling_, y_scaling_, x_spline, y_spline)
                sp_PeakE = calculatePeaKEnergy(y_scaling_, y_spline)

                # ::::: Fill values to the Variables :::::
                splineRange[0] = sp_Range
                splineTotalE[0] = sp_TotalE
                splinePeakE[0] = sp_PeakE
                x_scaling[0] = x_scaling_
                y_scaling[0] = y_scaling_
                chi2[0] = reduced_chi2
                avg3padE[0] = np.sum(filtered_dE_dX_Exp[:3]) / 3.0  # Average of the first 3 pads
                
                # print chi2 and avg3padE values.. later comment this
                print(f"Event {i_Event}: Reduced chi2 = {reduced_chi2}, avg3padE = {avg3padE[0]}")

        
        # ::::: Fill the branches with the calculated values :::::
        splineRange_branch.Fill()
        splineTotalE_branch.Fill()
        splinePeakE_branch.Fill()
        x_scaling_branch.Fill()
        y_scaling_branch.Fill()
        chi2_branch.Fill()
        avg3padE_branch.Fill()

    
    # Write the updated tree to the file    
    rootFile.cd()
    rootFile.Write("", ROOT.TObject.kOverwrite)
    rootFile.Close()
    print(f"Total NUmber of Events processed: {countGoodEvent}")





if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python3 fitSpline.py RUN_Nbr")
        sys.exit(1)

    RUN_Nbr = int(sys.argv[1])


    # RUN_Nbr = 1052
    start_time = datetime.now()
    main(RUN_Nbr)
    end_time = datetime.now()
    print(f"Time taken for processing: {end_time - start_time}")


















# ======================================================================
# # Fit only the non-zero part of experimental data
# mask_exp = dE_dX_Exp > 0
# x_data = ICSegment[mask_exp]
# y_target = dE_dX_Exp[mask_exp]

# start_time = datetime.now()
# x_scaling, y_scaling = fit_spline(dE_dX_Exp, x_spline, y_spline)
# # x_scaling, y_scaling = fit_spline_chi2_min(dE_dX_Exp, x_spline, y_spline)
# end_time = datetime.now()
# print(f"Time taken for fitting: {end_time - start_time}")

# calculateRange(x_scaling, y_scaling, x_spline, y_spline)

# calculateTotalEnergy(dE_dX_Exp, x_scaling, y_scaling, x_spline, y_spline)
# print("Experimental Total Energy:", np.sum(dE_dX_Exp))

# # calculateTotalEnergy2(x_scaling, y_scaling, x_spline, y_spline)

# calculatePeaKEnergy(y_scaling, y_spline)

# indices_x = [np.abs(y_spline - min(y_target)).argmin()]
# indices_y = [np.abs(x_spline*x_scaling - val).argmin() for val in x_data]


# scaled_x_spline = x_spline * x_scaling
# scaled_y_spline = y_spline * y_scaling
# xbin_means, xbin_edges, xbinnumber = binned_statistic(scaled_x_spline, scaled_x_spline, statistic='mean', bins=30)
# ybin_means, ybin_edges, ybinnumber = binned_statistic(scaled_x_spline, scaled_y_spline, statistic='mean', bins=30)


# # print("Indices XXXXX",indices_x)
# # print("Min y_target: ",min(y_target))
# # # print("Matched y_spline value at min y_target:", y_spline[indices_x])
# # print("Matched x_spline value at min y_target:", x_spline[indices_x])

# # print("Matched indices_y in x_spline:", indices_y[:-1])
# # print("Matched x_spline values:", x_spline[indices_y[:-1]])
# # print("Matched y_spline values:", y_spline[indices_y[:-1]])

# plt.figure(figsize=(10, 6))
# plt.plot(x_spline*x_scaling, y_spline*y_scaling, label='Spline Fit', color='orange')
# plt.scatter(ICSegment, dE_dX_Exp, label='Experimental dE/dX', color='blue')
# plt.scatter(x_spline[indices_y]*x_scaling, y_spline[indices_y]*y_scaling, label='Matched Spline Points', color='red')
# # plt.scatter(xbin_means, ybin_means, label='Matched Spline Points', color='purple')
# # plt.scatter(x_spline[indices_x]*x_scaling, y_spline[indices_x]*y_scaling, label='Min Spline Point', color='green')
# plt.grid(True)
# plt.legend()
# plt.show()
# ======================================================================


# ======================================================================
# Convert to interpolator for evaluating at scaled positions
# spline_interp = interp1d(x_spline, y_spline, bounds_error=False, fill_value=0)
# print(spline_interp)
# # Cost function
# def cost(params):
#     a, b = params
#     x_scaled = a * x_data
#     y_model = b * spline_interp(x_scaled)
#     return np.sum((y_model - y_target) ** 2)

# # Brute-force parameter search (grid search)
# rranges = (slice(0.8, 2, 0.005), slice(0.5, 2, 0.005))
# best_params = brute(cost, rranges, full_output=True, finish=None)
# a_opt, b_opt = best_params[0]
# print(f"Best a = {a_opt:.4f}, b = {b_opt:.4f}")

# # Plot fit
# x_fit = np.linspace(0, 29, 300)
# y_fit = b_opt * spline_interp(a_opt * x_fit)

# # y_fit = 1.2 * spline_interp(1.75 * x_fit)

# plt.figure(figsize=(10, 6))
# plt.plot(x_fit, y_fit, label='Brute-force Fit', color='green')
# plt.scatter(ICSegment, dE_dX_Exp, label='Experimental dE/dX', color='blue')
# plt.xlabel("IC Segment")
# plt.ylabel("dE/dX (MeV/unit)")
# plt.title("Fitting Scaled Spline to Experimental Data (Brute-force)")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()
# ======================================================================



# ======================================================================
# Don't Remove
# ICSegment = np.arange(30)
# mask_exp = dE_dX_Exp > 0
# # x_data = np.array(ICSegment[mask_exp])
# # y_target = np.array(dE_dX_Exp[mask_exp])
# x_data = np.array(ICSegment)
# y_target = np.array(dE_dX_Exp)
# print("x_data:", x_data)
# print("y_target:", y_target)

# x = ICSegment
# y = EnergyLoss
# x_spline = np.linspace(x.min(), x.max(), 300)
# spline_fn = make_interp_spline(x, y, k=3)
# y_spline = spline_fn(x_spline)

# # indices_x = [np.abs(y_spline - min(y_target)).argmin()]
# # x_scaling =  x_data[-1] / x_spline[indices_x]
# # indices_y = [np.abs(x_spline*x_scaling - val).argmin() for val in x_data]
# # matched_y_spline = y_spline[indices_y[:-1]]

# from scipy.stats import binned_statistic

# xbin_means, xbin_edges, xbinnumber = binned_statistic(x_spline, x_spline, statistic='mean', bins=30)
# ybin_means, ybin_edges, ybinnumber = binned_statistic(x_spline, y_spline, statistic='mean', bins=30)


# mask_spline_binnned = ybin_means > 1
# # print("Masked spline binned:", len(mask_spline_binnned))
# print("Binned statistic X:", xbin_means[mask_spline_binnned])
# print("Binned statistic Y:", ybin_means[mask_spline_binnned])
# print("Bin edges Y:", ybin_edges)
# # print("Bin number Y:", ybinnumber)
# bin_width = ybin_edges[1] - ybin_edges[0]



# # plt.figure(figsize=(10, 6))
# # plt.plot(x_spline, y_spline, label='Spline Fit', color='orange')
# # plt.scatter(xbin_means[mask_spline_binnned], ybin_means[mask_spline_binnned], label='Binned Statistic', color='red')
# # # plt.hlines(ybin_means, ybin_edges[:-1], ybin_edges[1:], colors='red', label='Binned Statistic', alpha=0.5)
# # # plt.plot((ybinnumber - 0.5) * bin_width, y_spline, 'g.', alpha=0.5)
# # plt.legend()
# # plt.grid(True)
# # plt.show()

# x_scaling_range = np.linspace(0.5, 2.0, 200)
# y_scaling_range = np.linspace(0.5, 2.0, 200)
# # Create a grid of parameters
# param_grid = np.array(list(product(x_scaling_range, y_scaling_range)))
# # To compute chi-square for given parameters
# # x_spline_matched = xbin_means[mask_spline_binnned]
# # y_spline_matched = ybin_means[mask_spline_binnned]
# x_spline_matched = xbin_means
# y_spline_matched = ybin_means

# dict_chi2 = {}

# for params in param_grid:
#     x_scaling, y_scaling = params
#     x_chi2 = np.sum((x_spline_matched * x_scaling - x_data) ** 2)
#     y_chi2 = np.sum((y_spline_matched * y_scaling - y_target) ** 2)
#     total_chi2 = x_chi2 + y_chi2
#     dict_chi2[tuple(params)] = total_chi2
# # Find the optimal parameters with the minimum chi-square
# optimal_params = min(dict_chi2, key=dict_chi2.get)
# x_scaling_opt, y_scaling_opt = optimal_params
# chi2_opt = dict_chi2[optimal_params]
# print("Optimal x_scaling:", x_scaling_opt)
# print("Optimal y_scaling:", y_scaling_opt)
# print("Optimal Chi-square:", chi2_opt)

# # Plot the results
# plt.figure(figsize=(10, 6))
# plt.plot(x_spline, y_spline, label='Spline Fit', color='orange')
# plt.scatter(xbin_means[mask_spline_binnned], ybin_means[mask_spline_binnned], label='Binned Statistic', color='red')
# plt.scatter(x_data, y_target, label='Experimental dE/dX', color='blue')
# plt.scatter(x_spline * x_scaling_opt, y_spline * y_scaling_opt, label='Scaled Spline Fit', color='green')
# plt.legend()
# plt.grid(True)
# plt.show()
# ====================================================================== 