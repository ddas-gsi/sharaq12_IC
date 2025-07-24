# ==================================================
# This code is to Fit Bragg Curve simulated from
# LISE++ simulation to the Experimental data 
# for all the events
# ==================================================


import ROOT
import numpy as np
import scipy 
import json
# import uproot     # not there in lxpool python installation
# import pandas     # not there in lxpool python installation

def calculateSlope(X1, Y1, X2, Y2):
    slope = (Y2 - Y1)/(X2 - X1)
    return slope

def calculateDy(Y1, Y2):
    dY = Y2 - Y1
    return dY

def getEnergyValForLowestChiSqr(chiSqrMap):
    EnergyValForLowestChiSqr = min(chiSqrMap, key = chiSqrMap.get)
    return EnergyValForLowestChiSqr

def getSegmentIDfromBraggCurve(ELossList, value):
    index = ELossList.index(value)
    return index

def generateICBranchNameList(ICSegmentList, branchNameStr):
    """Returns the IC Branch Names as a List of strings"""
    ICBranchNameList = []
    for segmentNbr in ICSegmentList:
        ICBranchName = branchNameStr + str(segmentNbr)
        ICBranchNameList.append(ICBranchName)
        # print(ICBranchName)
    return ICBranchNameList

def readEnergyLossData(filename):
    """
    Read Eloss data from txt file of LISE++ Simulation for different energies
    """
    # data = np.loadtxt(filename, delimiter="\t", skiprows=1)
    data = np.genfromtxt(filename, delimiter="\t", skip_header=1, dtype=float)
    
    energyArray = []
    segmentID = []
    for d in data:
        segmentID.append(d[0])
        energyArray.append(d[1])

    # print(energyArray)
    return segmentID,energyArray


# readEnergyLossData("/u/ddas/software/work/artemis-oedo/output/Analysis/12.0MeV.txt")


# define BeamEnergies
BeamEnergies = [12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5]

# Load the ELoss data and store in LISE Eloss Data Map/Dictionary
LISEElossDataMap = {}
LISEElossFilePath = "/u/ddas/software/work/artemis-oedo/output/Analysis/"
for energy in BeamEnergies:
    fileName = LISEElossFilePath + str(energy) + "MeV.txt"
    print(fileName)
    segmentID, Eloss = readEnergyLossData(fileName)
    LISEElossDataMap[energy] = Eloss

# print(LISEElossDataMap[17.5])

# ----------------------------------------------------------------------------------------

rootFile = ROOT.TFile.Open("/u/ddas/software/work/artemis-oedo/output/Analysis/902_All.root")
# tree = rootFile["tree_new"]   # this method does not work
# tree = rootFile.tree_new    # this method works... use the rootFile.treeName method
tree = rootFile.Get("tree_new")     # this method also works...
# print(type(tree))

# for branch in tree.GetListOfBranches():
#     print(branch.GetName())


# Create a TCanvas
# canvas = ROOT.TCanvas("c1", "IC Canvas", 2500, 750)
# Divide the canvas into 2 parts (1 row, 2 columns)
# canvas.Divide(2, 1)  # (columns, rows)

# Create a 2D histogram
# hist2d_IC_E_Cal = ROOT.TH2F("IC_E_Cal", "IC Calibrated Energy Data; Segment; IC Cal Energy", 31, -0.5, 30.5, 1000, 0, 50)

# i = 0
# for entry in tree:
    # print(type(entry))
#     # print(entry.IC_E_Cal_10)  # entry.BranchName
#     val = getattr(entry, "IC_E_Cal_10")    # this also works
#     print(val)  
#     i=i+1
# print(i)

ICSegmentList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
ICBranchNameList = generateICBranchNameList(ICSegmentList, "IC_E_Cal_")

# tree.GetEntry(363409)
# IC_E_Cal_ExptArr_iEvent = []
# for idx in range(len(ICBranchNameList)):
#     # print(getattr(tree, ICBranchNameList[idx]))
#     ic_E_Cal_idx = getattr(tree, ICBranchNameList[idx])
#     IC_E_Cal_ExptArr_iEvent.append(ic_E_Cal_idx)
    # hist2d_IC_E_Cal.Fill(idx, ic_E_Cal_idx)

# canvas.cd(1)
# hist2d_IC_E_Cal.Draw("BOX")

def calculateChiSqr(ElossExptArray, LISEElossDataMap, BeamEnergies):
    chiSqrMap = {}
    for energy in BeamEnergies:
        chiSqr = 0
        LISEElossAt_iEnergy = LISEElossDataMap[energy]
        for idx in range(len(ElossExptArray)-1):
            dYExptData = calculateDy(ElossExptArray[idx], ElossExptArray[idx+1])
            slopeLISESim = calculateSlope(idx, LISEElossAt_iEnergy[idx], idx+1, LISEElossAt_iEnergy[idx+1])

            if (slopeLISESim > 0 or dYExptData > -5.0):
                dY = ElossExptArray[idx] - LISEElossAt_iEnergy[idx]
                dYSquare = dY ** 2
                chiSqr = chiSqr + dYSquare

        # print(chiSqr)
        chiSqrMap[energy] = chiSqr
    
    return chiSqrMap


# chiSqrMap = calculateChiSqr(IC_E_Cal_ExptArr_iEvent, LISEElossDataMap, BeamEnergies)
# # print(chiSqrMap)
# energyValWithLowestChiSqr = getEnergyValForLowestChiSqr(chiSqrMap)
# print(f"Energy of Bragg Curve with Lowest ChiSqr: {energyValWithLowestChiSqr} and Value of ChiSqr: {chiSqrMap[energyValWithLowestChiSqr]}")

# Get the number of entries in the tree
num_entries = tree.GetEntries()
# Print the number of entries
print(f"Number of entries in the TTree: {num_entries}")

# define the container for tagging the Events
eventTaggingContainer = {}

for i_Event in range(num_entries):
    tree.GetEntry(i_Event)
    IC_E_Cal_ExptArr_iEvent = []
    for idx in range(len(ICBranchNameList)):
        # print(getattr(tree, ICBranchNameList[idx]))
        ic_E_Cal_idx = getattr(tree, ICBranchNameList[idx])
        IC_E_Cal_ExptArr_iEvent.append(ic_E_Cal_idx)

    chiSqrMap = calculateChiSqr(IC_E_Cal_ExptArr_iEvent, LISEElossDataMap, BeamEnergies)
    energyValWithLowestChiSqr = getEnergyValForLowestChiSqr(chiSqrMap)
    print(f"Event Number: {i_Event} >> Energy of Bragg Curve: {energyValWithLowestChiSqr} and Value of ChiSqr: {chiSqrMap[energyValWithLowestChiSqr]}")

    eventTaggingContainer[i_Event] = energyValWithLowestChiSqr


with open("/u/ddas/software/work/artemis-oedo/output/Analysis/EventTaggingJson.txt", "w") as eventTaggingFile:
    eventTaggingFile.write(json.dumps(eventTaggingContainer))

with open("/u/ddas/software/work/artemis-oedo/output/Analysis/EventTagging.csv", "w") as f:  
    for key, value in eventTaggingContainer.items():  
        f.write('%s,%s\n' % (key, value))


# Update the canvas
# canvas.Update()

# Close the root file
rootFile.Close()

# Keep the Canvas open in different ways
# Keep the canvas open
# input("Press Enter to exit...")
# Keep the canvas open with the ROOT event loop
# ROOT.gApplication.Run()