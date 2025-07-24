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
# 
# 
# ==================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import product
from datetime import datetime
import ROOT



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
    index = np.where(ElossArray<1)[0]
    index = index[0]
    dx = icStripTotalWidth / len(ElossArray)
    Range = dx * index
    return Range

def fnCalculateTotalEloss(ElossArray):
    """
    Returns total Energy loss in Ion-Chamber in MeV
    """
    index = np.where(ElossArray<0.5)[0]    # can change the threshold. Here I put 0.5
    index = index[0]
    newElossArray = ElossArray[:index]
    totalEloss = np.sum(newElossArray)
    return totalEloss

def chi_square(params, X0, X1, expData):
    Xpeak, Aexp, Bexp, sigmaGauss, ampGauss = params
    x_, model_y = fnExpGaussian(X0, X1, Xpeak, Aexp, Bexp, sigmaGauss, ampGauss)
    chi2 = np.sum(((expData - model_y)) ** 2)
    return chi2


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

# x,y = fnExpLorenz(0, 30, 15, 20, 0.2, 20, 1.2)
# x,y = fnExpGaussian(0,30,15,20,0.2,1.5,40)
# plt.plot(x,y)
# plt.show()



# ===================================================================================
# 
# ===================================================================================
script_start_time = datetime.now()

# Read the root file
rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root")
tree = rootFile.Get("tree_new")     # this method also works...
# Get the number of entries in the tree
num_entries = tree.GetEntries()
# Print the number of entries
print(f"Number of entries in the TTree: {num_entries}")

ICSegmentList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
ICBranchNameList = generateICBranchNameList(ICSegmentList, "IC_E_Cal_")

# Read through all the Events
countGoodEvent = 0
for i_Event in range(num_entries):
    tree.GetEntry(i_Event)
    IC_dE_dX_ExptArr_iEvent = []
    for idx in range(len(ICBranchNameList)):
        IC_dE_dX_idx = getattr(tree, ICBranchNameList[idx])
        IC_dE_dX_ExptArr_iEvent.append(IC_dE_dX_idx)

    # if(i_Event % 1000 == 0):
    #     print(f"{i_Event} entries done...")

    if IC_dE_dX_ExptArr_iEvent[0]>6.5:    # added a threshold of 6.5
        countGoodEvent=countGoodEvent+1
        IC_dE_dX_ExptArr_iEvent = np.array(IC_dE_dX_ExptArr_iEvent)
        idxLastHit = np.where(IC_dE_dX_ExptArr_iEvent<=7)[0]
        if len(idxLastHit)>0:
            idxLastHit = idxLastHit[0]    # Take the first occurrence
            IC_dE_dX_ExptArr_iEvent[idxLastHit:] = 0    # Set all elements from idxLastHit onward to 0
        print(IC_dE_dX_ExptArr_iEvent)
        print(f"Performing Brute-force optimization for {i_Event}-th event")
        start_time = datetime.now()
        best_params_iEvent = bruteForce_minmization(IC_dE_dX_ExptArr_iEvent)
        end_time = datetime.now()
        print(f"Opt time for {i_Event}-th event: {end_time-start_time}")
        Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best = best_params_iEvent
        print(f"Best Parameters: Xpeak = {Xpeak_best:.3f}, Aexp = {Aexp_best:.3f}, Bexp = {Bexp_best:.3f}, sigmaGauss = {sigmaGauss_best:.3f}, ampGauss = {ampGauss_best:.3f}")



script_end_time = datetime.now()
print(f"Script Run time for {countGoodEvent} good events: {script_end_time-script_start_time}")

# ===================================================================================

# dE_dX_Exp = np.array([
#     19.30883442006729, 20.07757334620357, 18.961192169495646, 20.863828861326464,
#     20.48170201163091, 21.36438847886429, 21.63018097893729, 21.78063574955428,
#     23.27558740839977, 24.935359534156735, 23.598545133143375, 26.10346947974618,
#     25.56415320656543, 28.99915234738814, 32.45695704909536, 36.3346556802946,
#     41.47629481050695, 42.45185952712907, 33.99725735284845, 8.003822782169488,
#     6.399985336097896, 0.0, 0.0, 0.0,
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0
# ])  # Add all points

# # Parameter ranges for brute force
# Xpeak_range = np.arange(5,25,1)
# Aexp_range = np.linspace(7,30,25)
# Bexp_range = np.linspace(0.1,0.7,7)
# sigmaGauss_range = np.linspace(1,1.5,5)
# ampGauss_range = np.linspace(30,45,15)

# # Perform Brute-force search
# count = 0
# best_params = None
# best_chi2 = float('inf')
# for params in product(Xpeak_range,Aexp_range,Bexp_range,sigmaGauss_range,ampGauss_range):
#     chi2 = chi_square(params, 0, 30, dE_dX_Exp)
#     count = count+1
#     # print(chi2)
#     if chi2 < best_chi2:
#         best_chi2 = chi2
#         # print(best_chi2)
#         best_params = params



# print(f"Good Events: {countGoodEvent}")
# print(f"Brute Force Optimization Counts: {count}")
# Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best = best_params
# print(f"Best Parameters: Xpeak = {Xpeak_best:.3f}, Aexp = {Aexp_best:.3f}, Bexp = {Bexp_best:.3f}, sigmaGauss = {sigmaGauss_best:.3f}, ampGauss = {ampGauss_best:.3f}")
# print(f"Minimum Chi-Square: {best_chi2:.3f}")

# x_sim, fitted_curve = fnExpGaussian(0,30,Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best)
# # print(fitted_curve)

# x_simSplined, fitted_curveSplined = fnExpGaussianSplined(0,30,1000,Xpeak_best, Aexp_best, Bexp_best, sigmaGauss_best, ampGauss_best)
# # print(fitted_curveSplined)

# print(f"Range (fine tuned) of the ion: {fnCalculateRange(fitted_curveSplined):0.5f} in mm")
# print(f"Range (Coarse) of the ion: {fnCalculateRange(fitted_curve):0.5f} in mm")

# print(f"Total Energy Loss (fine tuned): {fnCalculateTotalEloss(fitted_curveSplined):0.5f} MeV")
# print(f"Total Energy Loss (Coarse): {fnCalculateTotalEloss(fitted_curve):0.5f} MeV")
# print(f"Total Energy Loss (Experimental): {fnCalculateTotalEloss(dE_dX_Exp):0.5f} MeV")

# # index = np.where(fitted_curveSplined<0.5)[0]
# # index = index[0]
# # newElossArray = fitted_curveSplined[:index]
# # print(len(newElossArray))

# x_exp = np.arange(0, 30, 1)
# popt, pcov = curve_fit(fnExpGaussianFineCalculation, x_exp, dE_dX_Exp, p0=best_params)    # p0=best_params is the initial guess
# Xpeak_best_ls, Aexp_best_ls, Bexp_best_ls, sigmaGauss_best_ls, ampGauss_best_ls = popt
# print(f"Best Parameters from Least-Square: Xpeak = {Xpeak_best_ls:.3f}, Aexp = {Aexp_best_ls:.3f}, Bexp = {Bexp_best_ls:.3f}, sigmaGauss = {sigmaGauss_best_ls:.3f}, ampGauss = {ampGauss_best_ls:.3f}")

# plt.figure(figsize=(10, 6))
# plt.errorbar(x_exp, dE_dX_Exp, fmt='o', label='Experimental Data', color='red')
# plt.plot(x_sim, fitted_curve, label='Fitted Curve', color='blue', linewidth=2)
# plt.plot(x_simSplined,fitted_curveSplined,label='Smoothly Fitted Curve', linewidth=2)
# plt.plot(x_exp, fnExpGaussianFineCalculation(x_exp, *popt), color="black", label=" Fine Calculation Fitted Curve", linewidth=2)
# plt.legend()
# plt.show()

# ===================================================================================