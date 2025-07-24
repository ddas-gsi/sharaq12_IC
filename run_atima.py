import sys
import json
import atima


def calcRangeEloss(A, Z, Ebeam = 14.22):
   
   proj = [A, Z]

   # SRPPAC-Mylar: C₁₀H₈O₄ ===========================
   targ_mylar = [
    [12.011, 6, 10],       # Carbon
    [1.00784, 1, 8],       # Hydrogen
    [15.999, 8, 4]         # Oxygen
   ]
   density_mylar = 1.38                      # Mylar density on g/cm3
   isgas_mylar = 0                           # Is the target a gass 0/1 --> No/Yes
   targThicknessMylar = 9.936                # SRPPAC-Mylar thickness in mg/cm2.. This considers two layers of 0.0036 cm thick Mylar
   # targThicknessMylar = 9.936 / 2.0          # SRPPAC-Mylar thickness in mg/cm2.. This considers one layers of 0.0036 cm thick Mylar. Maybe this is the correct one

   # IC-Kapton Window: C₂₂H₁₀N₂O₄ ===========================
   targ_kapton = [
    [12.011, 6, 22],       # Carbon
    [1.00784, 1, 10],      # Hydrogen
    [14.0067, 7, 2],       # Nitrogen
    [15.999, 8, 4]         # Oxygen
   ]
   density_kapton = 1.42                     # Kapton density on g/cm3
   isgas_kapton = 0                          # Is the target a gass 0/1 --> No/Yes
   targThicknessKapton = 7.1                 # IC-Kapton window thickness in mg/cm2

   # CF4: C₁F₄ ==============================
   targ_CF4 = [
    [12.011, 6, 1],          # 1 Carbon atom
    [18.9984032, 9, 4]       # 4 Fluorine atoms
   ]
   density_CF4 = 0.000578                    # Target density on g/cm3
   isgas_CF4 = 1                             # Is the target a gass 0/1 --> No/Yes
   
   targThicknessICDeadLayer = 6.43025        # IC dead layer thickness in mg/cm2
   targThicknessICActiveRegion = 43.75       # IC Active region target thickness in mg/cm2
   # Ebeam_IC = 9.59            # Projectile energy in MeV/u

   # proj = [ 51, 20 ]                         # Projectile [Mass, Charge]
   # Ebeam = 14.22                             # Projectile energy in MeV/u

   # Start main calcualtion
   # Calculate the energy loss in the SRPPAC-Mylar
   (Ein_mylar, Eout_mylar, projrange ,dEdXin, dEdXout, RemainingRange,  \
      RangeStraggling , EnegergyStraggling, AngularStraggling, \
      Tof, InterpolatedTargetThickness) \
         = atima.eloss(proj, Ebeam, targ_mylar, density_mylar, isgas_mylar, targThicknessMylar)

   Ebeam = Eout_mylar

   # Calculate the energy loss in the IC-Kapton window
   (Ein_kapton, Eout_kapton, projrange ,dEdXin, dEdXout, RemainingRange,  \
      RangeStraggling , EnegergyStraggling, AngularStraggling, \
      Tof, InterpolatedTargetThickness) \
         = atima.eloss(proj, Ebeam, targ_kapton, density_kapton, isgas_kapton, targThicknessKapton)

   # Calculate the energy loss in the IC dead layer
   (Ein_ICDead, Eout_ICDead, projrange ,dEdXin, dEdXout, RemainingRange,  \
      RangeStraggling , EnegergyStraggling, AngularStraggling, \
      Tof, InterpolatedTargetThickness) \
         = atima.eloss(proj, Eout_kapton, targ_CF4, density_CF4, isgas_CF4, targThicknessICDeadLayer)

   # Calculate the energy loss in the IC active region
   (Ein_ICActive, Eout_ICActive, projrange ,dEdXin, dEdXout, RemainingRange,  \
      RangeStraggling , EnegergyStraggling, AngularStraggling, \
      Tof, InterpolatedTargetThickness) \
         = atima.eloss(proj, Eout_ICDead, targ_CF4, density_CF4, isgas_CF4, targThicknessICActiveRegion)

   # # Print the results
   # print('Ein_mylar [MeV/u]: ' , Ein_mylar)
   # print('Eout_mylar [MeV/u]: ' ,Eout_mylar)
   # print('Ein_kapton [MeV/u]: ' , Ein_kapton)
   # print('Eout_kapton [MeV/u]: ' ,Eout_kapton)
   # print('Ein_ICDead [MeV/u]: ' , Ein_ICDead)
   # print('Eout_ICDead [MeV/u]: ' ,Eout_ICDead)
   # print('Ein_ICActive [MeV/u]: ' , Ein_ICActive)
   # print('Eout_ICActive [MeV/u]: ' ,Eout_ICActive)

   # print('Range [mg/cm2]: ',projrange)   
   # # print('dEdX_in [MeV/mg/cm2]: ' ,dEdXin*proj[0])
   # # print('dEdX out [MeV/mg/cm2]: ' ,dEdXout*proj[0])
   # print('Remaing Range[mg/cm2]  : ' ,RemainingRange)
   # print('Range straggling [mg/cm2]  : ' ,RangeStraggling)
   # # print('Energy Straggling  [MeV] : ' ,EnegergyStraggling*dEdXout*proj[0])
   # # print('Angular Straggling [mrad]  : ' ,AngularStraggling*1000.)
   # # print('Time of Flight [ns] : ' ,Tof)
   # print('Interpolated Target thickness [mg/cm2]: ', InterpolatedTargetThickness)
   # print('IC Active Region Thickness [mg/cm2]: ', targThicknessICActiveRegion)   

   # print()
   # print('Energy Loss in IC Active Region [MeV]:',(Ein_ICActive-Eout_ICActive)*proj[0])
   # print()

   range_mm = (projrange / (density_CF4 * 1000)) * 10.0    # range in mm
   eloss = (Ein_ICActive-Eout_ICActive)*proj[0]  # Energy loss in MeV

   # print('Range [mm]: ', range_mm)
   # print('Total Energy Loss in IC Active Region [MeV]: ', eloss)

   return {"range": range_mm, "eloss": eloss}




if __name__ == "__main__":   
   A = int(float(sys.argv[1]))
   Z = int(float(sys.argv[2]))
   Beam_E = float(sys.argv[3])

   try:
      result = calcRangeEloss(A,Z, Ebeam=Beam_E)
      # print(result)
      print(json.dumps(result))

   except Exception as e:
      print(f"Error: {str(e)}", file=sys.stderr)
      sys.exit(1)
