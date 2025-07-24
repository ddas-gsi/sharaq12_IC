#!/usr/bin/python3

import atima, math
import numpy as np
from itertools import product
from datetime import datetime
from multiprocessing import Pool
from functools import partial
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import os
import csv
import ROOT
import subprocess
import json

cpu_usage = 64

def worker_func_range(params, corrected_expt_Range):
    Z, A, Beam_E = params
    try:
        result = subprocess.run(
            ['python3', 'run_atima.py', str(A), str(Z), str(Beam_E)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )

        output = None
        for line in result.stdout.splitlines():
            try:
                output = json.loads(line)
                print(output)  # Clean final output
                break  # Stop at first valid JSON line
            except json.JSONDecodeError:
                continue

        if output is None:
            print(f"No valid JSON output for A={A}, Z={Z}")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")
            return (params, float("inf"), None)

        calc_range = output["range"]
      #   calc_eloss = output["eloss"]
        chi2_range = (corrected_expt_Range - calc_range) ** 2
        return (params, chi2_range, output)

    except subprocess.CalledProcessError as e:
        print(f"Subprocess error for A={A}, Z={Z}:\n{e.stderr}")
        return (params, float("inf"), float("inf"), None)


def rangeFit(Beam_E, corrected_expt_Range):
   #  Z_range = np.linspace(18, 21, 31)
   #  A_range = np.linspace(48, 52, 41)
    Z_range = np.linspace(18, 21, 16)
    A_range = np.linspace(48, 52, 21)
    param_combinations = list(product(Z_range, A_range, [Beam_E]))

    print("Starting calculation...")
    print(f"Number of combinations: {len(param_combinations)}")
    startTime = datetime.now()

    results = []
   #  with ThreadPoolExecutor(max_workers=4) as executor:
   #      future_to_params = {executor.submit(chi2, params, corrected_expt_Range): params for params in param_combinations}
   #      for future in as_completed(future_to_params):
   #          result = future.result()
   #          results.append(result)

    with ProcessPoolExecutor(max_workers=cpu_usage) as executor:

      future_to_params = {executor.submit(worker_func_range, params, corrected_expt_Range): params for params in param_combinations}
      for future in as_completed(future_to_params):
         result = future.result()
         results.append(result)

    best_result = min(results, key=lambda x: x[1])  # x = (params, chi2, output)
    best_params, best_chi2, best_output = best_result

    print(f"Total time taken: {datetime.now() - startTime}")
    print(f'Best Params: Z: {best_params[0]}, A: {best_params[1]}, Beam_E: {best_params[2]}')
    print(f"Best chi2: {best_chi2}")
    if best_output:
        print(f"Best range [mm]: {best_output['range']}")
        print(f"Best energy loss [MeV]: {best_output['eloss']}")

    return best_params

def worker_func_rangeEloss(params, corrected_expt_Range, IC_Eloss):
    Z, A, Beam_E = params
    try:
        result = subprocess.run(
            ['python3', 'run_atima.py', str(A), str(Z), str(Beam_E)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )

        output = None
        for line in result.stdout.splitlines():
            try:
                output = json.loads(line)
                # print(output)  # Clean final output
                break  # Stop at first valid JSON line
            except json.JSONDecodeError:
                continue

        if output is None:
            print(f"No valid JSON output for A={A}, Z={Z}")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")
            return (params, float("inf"), None)

        calc_range = output["range"]
        calc_eloss = output["eloss"]
        chi2_range = (corrected_expt_Range - calc_range) ** 2
        chi2_eloss = (IC_Eloss - calc_eloss) ** 2
        chi2_total = chi2_range + chi2_eloss
        return (params, chi2_range, chi2_total, output)

    except subprocess.CalledProcessError as e:
        print(f"Subprocess error for A={A}, Z={Z}:\n{e.stderr}")
        return (params, float("inf"), float("inf"), None)

def rangeElossFit(Z_range, A_range, Beam_E, corrected_expt_Range, IC_Eloss):
    # Z_range = np.linspace(18, 21, 16)
    # A_range = np.linspace(48, 52, 21)

    param_combinations = list(product(Z_range, A_range, [Beam_E]))

    print("Starting calculation...")
    print(f"Number of combinations: {len(param_combinations)}")
    startTime = datetime.now()

    results = []
   #  with ThreadPoolExecutor(max_workers=4) as executor:
   #      future_to_params = {executor.submit(chi2, params, corrected_expt_Range): params for params in param_combinations}
   #      for future in as_completed(future_to_params):
   #          result = future.result()
   #          results.append(result)

    with ProcessPoolExecutor(max_workers=cpu_usage) as executor:

      future_to_params = {executor.submit(worker_func_rangeEloss, params, corrected_expt_Range, IC_Eloss): params for params in param_combinations}
      for future in as_completed(future_to_params):
         result = future.result()
         results.append(result)

    best_result = min(results, key=lambda x: x[2])  # x = (params, chi2_range, chi2_total, output)
    best_params, best_chi2_range, best_chi2_total, best_output = best_result

    print(f"Total time taken: {datetime.now() - startTime}")
    print(f'Best Params: Z: {best_params[0]}, A: {best_params[1]}, Beam_E: {best_params[2]}')
    print(f"Best chi2: {best_chi2_total}")
    if best_output:
        print(f"Best range [mm]: {best_output['range']}")
        print(f"Best energy loss [MeV]: {best_output['eloss']}")

    return best_params

# rangeFit(Beam_E = 14.22, corrected_expt_Range = 410)
# Z_range = np.linspace(19, 21, 11)
# A_range = np.linspace(49, 51, 11)
# best_params = rangeElossFit(Z_range, A_range, Beam_E = 14.22, corrected_expt_Range = 410, IC_Eloss = 427.0)
# print(best_params)



if __name__ == "__main__":
    
    RUN_Nbr = 1052

    inputFileName = "/u/ddas/software/work/artemis-oedo/output/Analysis/physOutput/outputRange/rangeCorrection/AoQ/" + str(RUN_Nbr) + "_Analysis_2024.root"
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
    
    FE9cut50Ca20 = ROOT.gROOT.FindObject("FE9pidCut_50Ca_1053_2024")
    FE9cut51Sc21 = ROOT.gROOT.FindObject("FE9pidCut_51Sc_1053_2024")
    FE9cut49K19 = ROOT.gROOT.FindObject("FE9pidCut_49K_1053_2024")


    # Enable all branches first
    tree.SetBranchStatus("*", 1)

    # Define new variables for the new branches
    Z = np.zeros(1, dtype=np.float64)
    A = np.zeros(1, dtype=np.float64)

    # Add new branches to the tree
    Z_branch = tree.Branch("Z", Z, "Z/D")
    A_branch = tree.Branch("A", A, "A/D")

    for i_Event in range(num_entries):
        tree.GetEntry(i_Event)

        if (i_Event % 100000 == 0):
            print(f"Event Number: {i_Event}")

        # Initialize Z and A to -1000000.0 for each event
        Z[0] = -1000000.0
        A[0] = -1000000.0

        # Get the values of the branches
        PID_T_0 = tree.PID_T_0
        FE9_X_0 = tree.FE9_X_0
        deltaAll_stS1Y = tree.deltaAll_stS1Y
        betaS1_0 = tree.betaS1_0
        Beam_E_2_0 = tree.Beam_E_2_0
        IC_Eloss = tree.totalE 
        corrected_expt_Range = deltaAll_stS1Y * betaS1_0**2 * 100

        if (Beam_E_2_0>1.0 and IC_Eloss> 1.0 and corrected_expt_Range>1.0):
            count += 1

            # ***** Check which cut the event falls into and set Z and A accordingly *****

            if FE9cut50Ca20.IsInside(PID_T_0, FE9_X_0):    # For 50Ca20
                print("=================== Inside 50Ca20 cut =====================")
                Z_range = np.linspace(19, 21, 11)
                A_range = np.linspace(49, 51, 11)
                Beam_E = Beam_E_2_0
                best_params = rangeElossFit(Z_range, A_range, Beam_E, corrected_expt_Range, IC_Eloss)
                Z[0] = best_params[0]
                A[0] = best_params[1]
            elif FE9cut51Sc21.IsInside(PID_T_0, FE9_X_0):  # For 51Sc21
                print("=================== Inside 51Sc21 cut =====================")
                Z_range = np.linspace(20, 22, 11)
                A_range = np.linspace(50, 52, 11)
                Beam_E = Beam_E_2_0
                best_params = rangeElossFit(Z_range, A_range, Beam_E, corrected_expt_Range, IC_Eloss)
                Z[0] = best_params[0]
                A[0] = best_params[1]
            elif FE9cut49K19.IsInside(PID_T_0, FE9_X_0):  # For 49K19
                print("=================== Inside 49K19 cut =====================")
                Z_range = np.linspace(18, 20, 11)
                A_range = np.linspace(48, 50, 11)
                Beam_E = Beam_E_2_0 
                best_params = rangeElossFit(Z_range, A_range, Beam_E, corrected_expt_Range, IC_Eloss)
                Z[0] = best_params[0]
                A[0] = best_params[1]
                

        # ::: Fill the new branches with values :::
        Z_branch.Fill()
        A_branch.Fill()


    # Write the new branches to the tree
    rootFile.cd()
    tree.Write("", ROOT.TObject.kOverwrite)
    # Close the ROOT file
    rootFile.Close()
    print(f"Total number of events processed: {count}")
        
















