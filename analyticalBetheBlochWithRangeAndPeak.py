import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import product
from datetime import datetime
from multiprocessing import Pool
from functools import partial
import os
import ROOT
import ctypes

# ===================================================================================
available_cores = os.cpu_count()
print(f"Number of available CPU cores: {available_cores}")
use_cores = int(available_cores*0.13)
print(f"Running calculation on {use_cores} CPU cores. To increase number of cores, increase threshold.")
# ===================================================================================


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
    index = np.where(ElossArray<0.1)[0]   # used threshold 0.1 from where we consider the last hit happened
    # index = np.where(ElossArray<0.1)[0] + 1
    if index.size == 0:
        # index = len(ElossArray)-1
        index = len(ElossArray)
    else:
        index = index[0]
    print("index", index)
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

    # # Parameter ranges for fine brute force
    # Xpeak_range = np.arange(5,27,1)
    # Aexp_range = np.linspace(7,30,25)
    # Bexp_range = np.linspace(0.05,0.7,20)
    # sigmaGauss_range = np.linspace(1,1.8,10)
    # ampGauss_range = np.linspace(30,45,15)

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

def BraggPeakFunction(x, par):
    X = x[0]
    Xpeak, Aexp, Bexp, sigmaGauss, ampGauss = par
    if X < Xpeak:
        return Aexp + np.exp(Bexp * X)
    else:
        return ampGauss * np.exp(-((X - Xpeak) ** 2) / (2 * sigmaGauss ** 2))

def minuitBraggFit(IC_dE_dX_ExptArr_iEvent):
    # Find initial estimates
    dEdX_Max_i = np.max(IC_dE_dX_ExptArr_iEvent)
    Xpeak_i = np.argmax(IC_dE_dX_ExptArr_iEvent)
    Aexp_avg = np.mean(IC_dE_dX_ExptArr_iEvent[:4])
    
    # ROOT fit setup
    n = len(IC_dE_dX_ExptArr_iEvent)
    X_vals = np.arange(n, dtype=float)
    Y_vals = np.array(IC_dE_dX_ExptArr_iEvent, dtype=float)

    graph = ROOT.TGraph(n, X_vals, Y_vals)

    # Define Bragg Peak function properly using TF1 with C++ syntax
    braggFit = ROOT.TF1("braggFit", "[1] + exp([2] * x) * (x < [0]) + [4] * exp(-((x - [0]) * (x - [0])) / (2 * [3] * [3])) * (x >= [0])", 0, n, 5)

    # Set initial parameter guesses
    braggFit.SetParameters(Xpeak_i, Aexp_avg, 0.1, 1.0, dEdX_Max_i)
    braggFit.SetParNames("Xpeak", "Aexp", "Bexp", "sigmaGauss", "ampGauss")
    
    # Set parameter ranges
    braggFit.SetParLimits(0, 0, n-1)                        # Xpeak range 
    braggFit.SetParLimits(1, Aexp_avg-6, Aexp_avg+4)        # Aexp range
    braggFit.SetParLimits(2, 0.05, 0.7)                     # Bexp range
    braggFit.SetParLimits(3, 0.2, 1.5)                      # sigmaGauss range
    braggFit.SetParLimits(4, dEdX_Max_i-5, dEdX_Max_i+4)    # ampGauss range

    # Perform Minuit Fit
    fit_result = graph.Fit(braggFit, "RS")  # "RQS": Range, Quiet, Save
    
    print(f"Fit status: {fit_result.Status()}")
    if fit_result.Status() != 0:
        print("Fit did not converge properly.")

    # Get optimized parameters
    best_params = [braggFit.GetParameter(i) for i in range(5)]
    
    return best_params



# Define custom Bragg Peak function in C++
ROOT.gInterpreter.Declare("""
double BraggPeak(double *x, double *par)
{
    double X = x[0];
    double X0 = par[0];         // Start range
    double Xpeak = par[1];      // Peak position
    double Aexp = par[2];       // Exponential coefficient A
    double Bexp = par[3];       // Exponential coefficient B
    double sigmaGauss = par[4]; // Gaussian sigma
    double ampGauss = par[5];   // Gaussian amplitude

    // Exponential part: X < Xpeak
    if (X < Xpeak)
    {
        return Aexp + std::exp(Bexp * X);
    }
    // Gaussian part: X >= Xpeak
    else
    {
        return ampGauss * std::exp(-((X - Xpeak) * (X - Xpeak)) / (2 * sigmaGauss * sigmaGauss));
    }
}
""")

def minuitBraggFit2(IC_dE_dX_ExptArr_iEvent):
    dEdX_Max_i = np.max(IC_dE_dX_ExptArr_iEvent)
    Xpeak_i = np.argmax(IC_dE_dX_ExptArr_iEvent)
    Aexp_avg = np.mean(IC_dE_dX_ExptArr_iEvent[:4])
    
    # Convert NumPy arrays to ROOT-compatible arrays
    n = len(IC_dE_dX_ExptArr_iEvent)
    X_vals = ROOT.std.vector('double')(list(np.arange(n, dtype=float)))
    Y_vals = ROOT.std.vector('double')(list(IC_dE_dX_ExptArr_iEvent))

    graph = ROOT.TGraph(n, X_vals.data(), Y_vals.data())

    # Use the new custom function
    braggFit = ROOT.TF1("braggFit", ROOT.BraggPeak, 0, n, 6)

    # Set initial parameter estimates
    braggFit.SetParameters(0, Xpeak_i, Aexp_avg, 0.1, 1.0, dEdX_Max_i)
    braggFit.SetParNames("X0", "Xpeak", "Aexp", "Bexp", "sigmaGauss", "ampGauss")

    # Set parameter ranges
    braggFit.SetParLimits(0, 0, 10)                         # X0 range should be 0 and 10
    braggFit.SetParLimits(1, 0, n-1)                        # Xpeak range 
    braggFit.SetParLimits(2, Aexp_avg-6, Aexp_avg+4)        # Aexp range
    braggFit.SetParLimits(3, -1, 1)                         # Bexp range
    braggFit.SetParLimits(4, 0.2, 5)                        # sigmaGauss range
    braggFit.SetParLimits(5, dEdX_Max_i-5, dEdX_Max_i+4)    # ampGauss range

    # Perform the Minuit fit
    fit_result = graph.Fit(braggFit, "RS")  # "RQS": Range, Save, Quiet

    print(f"Fit status: {fit_result.Status()}")
    if fit_result.Status() != 0:
        print("Fit did not converge properly.")

    # Extract optimized parameters
    best_params = [braggFit.GetParameter(i) for i in range(6)]

    return best_params


# ==================================================================================================
# x,y = fnExpLorenz(0, 30, 15, 20, 0.2, 20, 1.2)
# x,y = fnExpGaussian(0,30,15,20,0.2,1.5,40)
# plt.plot(x,y)
# plt.show()

# dE_dX_Exp = [
#     19.30883442006729, 20.07757334620357, 18.961192169495646, 20.863828861326464,
#     20.48170201163091, 21.36438847886429, 21.63018097893729, 21.78063574955428,
#     23.27558740839977, 24.935359534156735, 23.598545133143375, 26.10346947974618,
#     25.56415320656543, 28.99915234738814, 32.45695704909536, 36.3346556802946,
#     41.47629481050695, 42.45185952712907, 33.99725735284845, 8.003822782169488,
#     6.399985336097896, 0.0, 0.0, 0.0,
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0
# ]  # Add all points


# dE_dX_Exp = [23.15828857, 26.25686581, 27.12598498, 27.88286197, 29.55377179, 30.43793649,
#  32.06677529, 33.60875744, 36.65823772, 39.9846389,  36.14198983, 33.06684032,
#  10.34674719,  0.,          0.,         0.,          0.,          0.,
#   0.,          0. ,         0.,          0.,          0.    ,      0.,
#   0.    ,      0.  ,        0. ,         0.    ,      0.   ,       0.        ]  # Add all points

# dE_dX_Exp = [11.28088047, 12.95792504, 12.30882978, 14.59142329, 13.69589249, 15.24906241,
#  14.77191947, 15.86285996, 16.05677303, 17.09337617, 16.3128902,  17.728973,
#  17.87914965, 18.631978,   21.41606764, 22.879183,   25.84038813, 26.74905082,
#  30.08926704, 31.2520464,  35.60608161, 37.04810821, 38.32255692, 19.71487461,
#   0.,          0. ,         0.  ,        0.      ,    0.    ,      0.        ]

# dE_dX_Exp =  [ 9.168516,   10.01462886, 10.2900509,  11.87485822, 11.7350413 , 12.04700109,
#  11.8688846,  12.40485406, 13.97172702 ,13.91271879 ,13.48389765 ,14.27730995,
#  13.60621551 ,16.68141319, 17.82849901 ,20.74381663, 21.33160314, 22.05066379,
#  24.77273985, 26.24818362 ,31.27253815 ,31.73278556 ,36.89880655 ,35.82970397,
#  18.34126153 , 0.  ,        0.     ,     0.  ,        0. ,         0.        ]

# dE_dX_Exp = np.array([
#     6.86746039, 6.86174839, 7.22660142, 7.18102129, 6.69747302, 7.40070116,
#     8.28359464, 6.86598553, 8.83833186, 8.46581126, 8.97176758, 9.1037164,
#     8.74981295, 9.48991126, 9.94919109, 11.67198452, 12.0182097, 12.15037258,
#     13.99077146, 13.37034899, 15.32036894, 15.31188849, 16.1815145, 18.00098574,
#     20.59605422, 22.6842125, 20.00806372, 6.1574421, 1.24400165, 1.30416514
# ])

# dE_dX_Exp = np.array([
#     7.99115087, 8.8731933, 8.62346948, 9.3155056, 9.10757759, 9.82433461,
#     9.82761015, 10.21501279, 10.865725, 10.83596941, 11.11509707, 11.70866108,
#     11.44794399, 12.60570605, 13.85539448, 14.68005782, 16.15993464, 15.47958455,
#     17.659119, 18.28436637, 19.52873129, 20.00990523, 22.19400068, 22.51017344,
#     27.32119477, 28.96951372, 31.4497906, 32.6843946, 21.50476169, 1.18248749
# ])    # event Nbr 415653 of 1005 file

# ==================================================================================================
# To extract one event
# ==================================================================================================
ICSegmentList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
ICBranchNameList = generateICBranchNameList(ICSegmentList, "IC_E_Cal_")

# Read the root file
# rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root")
# rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_1053new.root");
# rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_1058new.root");
rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/24sharaq12phys_1063new.root");
tree = rootFile.Get("tree_new")     # this method also works...

# 2024-1063 file events
tree.GetEntry(1996110)    
# tree.GetEntry(1996114)

# 2024-1058 file events
# tree.GetEntry(1989009)  
# tree.GetEntry(1988770)  

# 2024-1053 file events
# tree.GetEntry(863569)    # ~ 15.4932 MeV/u   can see double peak
# tree.GetEntry(95853)       # ~ 13.5 MeV/u
# tree.GetEntry(1797474)    # 13.4825 MeV/u
# tree.GetEntry(1361440)   # ~ 13.5 MeV/u

# 2022-902_All file events
# tree.GetEntry(363409)
# tree.GetEntry(358)
# tree.GetEntry(10)
# tree.GetEntry(660)

dE_dX_Exp = []
for idx in range(len(ICBranchNameList)):
    IC_dE_dX_idx = getattr(tree, ICBranchNameList[idx])
    dE_dX_Exp.append(IC_dE_dX_idx)
dE_dX_Exp = np.array(dE_dX_Exp)
print("RAW array: ")
print(dE_dX_Exp)

# Calculate the last hit index. We take average of last 3 channels, and then in the entire array, 
# the index where we have value less than average+0.5, we consider it as last hit index
# avgOfLast3Channel = np.average(dE_dX_Exp[27:])
# idxLastHit = np.where(dE_dX_Exp<=avgOfLast3Channel+0.2)[0]
# avgOfLast4Channel = np.average(dE_dX_Exp[26:])
# idxLastHit = np.where(dE_dX_Exp<=avgOfLast4Channel+0.2)[0]            
# if len(idxLastHit)>0:
#     idxLastHit = idxLastHit[0]                  # Take the first occurrence
#     dE_dX_Exp[idxLastHit:] = 0    # Set all elements from idxLastHit onwards to 0
# print("Filtered array: ")
# print(dE_dX_Exp)

# Find the baseline
# baseline_indices = find_baseline(dE_dX_Exp, threshold, window_size, tail_length)
# if len(baseline_indices) == 0:
#     raise ValueError("No sufficient flat region found at the end. Try adjusting 'threshold' or 'min_points'.")
# elif len(baseline_indices) >= 1:
#     idxLastHit = baseline_indices[0]              # Take the first occurrence
#     dE_dX_Exp[idxLastHit:] = 0    # Set all elements from idxLastHit onwards to 0

# params for baseline correction
window_size = 2
threshold = 0.3  # You can adjust this to define what you consider "flat"

dE_dX_Exp = fn_baseline_correction(dE_dX_Exp, window_size, threshold)
print("Filtered array: ")
print(dE_dX_Exp)

# ==================================================================================================
# Brute Force Optimization
# ==================================================================================================

# Define some parameters for the brute-force space
dEdX_Max_i = np.max(dE_dX_Exp)
Xpeak_i = np.where(dE_dX_Exp==dEdX_Max_i)[0]
Xpeak_i = Xpeak_i[0]
Aexp_avg = np.average(dE_dX_Exp[:4])
# Parameter ranges
Xpeak_range = np.arange(Xpeak_i-2,Xpeak_i+3, 1)
Aexp_range = np.linspace(Aexp_avg-6, Aexp_avg+4, 11)
Bexp_range = np.linspace(0.05,0.7,20)    # np.linspace(0.1, 0.7, 7)
sigmaGauss_range = np.linspace(0.6, 1.5, 10)  # np.linspace(0.8, 1.5, 8)
ampGauss_range = np.linspace(dEdX_Max_i-5, dEdX_Max_i+4, 10)

start_time = datetime.now()
best_params = bruteForce_minimization_parallel(dE_dX_Exp)
end_time = datetime.now()
print(f"Opt time for the event: {end_time-start_time}")


# print(count)
Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best = best_params
print(f"Best Parameters: Xpeak = {Xpeak_best:.3f}, Aexp = {Aexp_best:.3f}, Bexp = {Bexp_best:.3f}, sigmaGauss = {sigmaGauss_best:.3f}, ampGauss = {ampGauss_best:.3f}")
# print(f"Minimum Chi-Square: {best_chi2:.3f}")

x_sim, fitted_curve = fnExpGaussian(0,30,Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best)

x_simSplined, fitted_curveSplined = fnExpGaussianSplined(0,30,1000,Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best)
# print(fitted_curveSplined)

# Get the Range
print(f"Range (Experimental) of the ion: {fnCalculateRange(np.array(dE_dX_Exp)):0.5f} in mm")
print(f"Range (Coarse) of the ion: {fnCalculateRange(fitted_curve):0.5f} in mm")
print(f"Range (fine tuned) of the ion: {fnCalculateRange(fitted_curveSplined):0.5f} in mm")

# Get the total Eloss
print(f"Total Eloss from Experimental Data: {fnCalculateTotalEloss(np.array(dE_dX_Exp))} MeV")
print(f"Total Eloss from Simulated Data: {fnCalculateTotalEloss(fitted_curve)} MeV")

x_exp = np.arange(0, 30, 1)

plt.figure(figsize=(10, 6))
plt.errorbar(x_exp, dE_dX_Exp, fmt='o', label='Experimental Data', color='red')
plt.plot(x_sim, fitted_curve, label='Fitted Curve', color='blue', linewidth=2)
plt.plot(x_simSplined,fitted_curveSplined,label='Smoothly Fitted Curve', linewidth=2)
plt.legend()
plt.show()
# ==================================================================================================


# # ==================================================================================================
# # Minuit Fit
# # ==================================================================================================

# start_time = datetime.now()
# # best_params = minuitBraggFit(dE_dX_Exp)
# best_params = minuitBraggFit2(dE_dX_Exp)
# end_time = datetime.now()
# print(f"Opt time for the event: {end_time-start_time}")

# X0, Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best = best_params
# print(f"Best Parameters from Minuit Bragg Fit: Xpeak = {Xpeak_best:.3f}, Aexp = {Aexp_best:.3f}, Bexp = {Bexp_best:.3f}, sigmaGauss = {sigmaGauss_best:.3f}, ampGauss = {ampGauss_best:.3f}")

# # braggFit = ROOT.TF1("braggFit", "[1] + exp([2] * x) * (x < [0]) + [4] * exp(-((x - [0]) * (x - [0])) / (2 * [3] * [3])) * (x >= [0])", 0, 29)
# braggFit = ROOT.TF1("braggFit", ROOT.BraggPeak, 0, 30, 6)
# braggFit.SetParameters(X0, Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best)
# braggFit.SetParNames("X0", "Xpeak", "Aexp", "Bexp", "sigmaGauss", "ampGauss")

# canvas = ROOT.TCanvas("canvas", "Bragg Peak Fit", 800, 600)

# graph = ROOT.TGraph(len(dE_dX_Exp), np.arange(len(dE_dX_Exp),dtype=float), np.array(dE_dX_Exp,dtype=float))

# graph.SetMarkerStyle(20)  # Set point style
# graph.SetMarkerSize(1.2)  # Set point size
# graph.SetMarkerColor(ROOT.kRed)  # Set color to red
# graph.Draw("AP")  # Draw the experimental data
# braggFit.Draw("same")
# # canvas.Draw()
# canvas.Update()
# ROOT.gPad.Update()
# ROOT.gPad.WaitPrimitive()
# # ==================================================================================================