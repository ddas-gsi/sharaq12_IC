#!/usr/bin/env python3

from scipy.interpolate import make_interp_spline
import numpy as np
import matplotlib.pyplot as plt

# for 16 MeV/u 50Ca20 ion
ICSegment = list(range(30))
EnergyLoss = [
    18.562, 18.961, 19.404, 19.881, 20.395, 20.945, 21.548, 22.208, 22.933, 23.734,
    24.626, 25.628, 26.761, 28.058, 29.558, 31.303, 33.371, 35.812, 38.078, 32.5,
    6.919, 0, 0, 0, 0, 0, 0, 0, 0, 0
]

x = np.array(ICSegment)
y = np.array(EnergyLoss)

# Only interpolate non-zero region to avoid artifacts
# x_nonzero = x[y > 0]
# y_nonzero = y[y > 0]
x_nonzero = x
y_nonzero = y

# print("x_nonzero:", x_nonzero)
# print("y_nonzero:", y_nonzero)

# Create smooth x values
x_smooth = np.linspace(x_nonzero.min(), x_nonzero.max(), 300)
spline = make_interp_spline(x_nonzero, y_nonzero, k=3)
y_smooth = spline(x_smooth)

# print(x_smooth)
# print(y_smooth)

# plt.figure(figsize=(10, 6))
# plt.plot(x_smooth, y_smooth, label='Smoothed Energy Loss', color='green')
# plt.plot(x_smooth*0.5, y_smooth*1.2, label='Smoothed Energy Loss 2', color='blue')
# plt.plot(x_smooth, y_smooth*1.1, label='Smoothed Energy Loss 2', color='red')
# plt.scatter(ICSegment, EnergyLoss, color='blue', label='Original Data')
# plt.grid(True)
# plt.xlabel("ICSegment")
# plt.ylabel("EnergyLoss (MeV)")
# plt.title("Smoothed Parametric Curve")
# plt.legend()
# plt.tight_layout()
# plt.show()


dE_dX_Exp = np.array([
    23.15828857, 26.25686581, 27.12598498, 27.88286197, 29.55377179, 30.43793649,
    32.06677529, 33.60875744, 36.65823772, 39.9846389, 36.14198983, 33.06684032,
    10.34674719, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
])

plt.figure(figsize=(10, 6))
plt.plot(x_smooth, y_smooth, label='Parametric Curve', color='green')
plt.plot(x_smooth*0.6, y_smooth*1.15, label='Fitted Smoothed Energy Loss Curve', color='red')
plt.scatter(ICSegment, dE_dX_Exp, color='blue', label='Original Data')
plt.grid(True)
plt.xlabel("ICSegment")
plt.ylabel("EnergyLoss (MeV)")
plt.title("Smoothed Parametric Curve")
plt.legend()
plt.tight_layout()
plt.show()





# # part 2
# x_nonzero2 = np.array(list(range(10)))

# # Create smooth x values
# x_smooth2 = np.linspace(x_nonzero2.min(), x_nonzero2.max(), 300)
# # spline = make_interp_spline(x_nonzero, y_nonzero, k=3)
# y_smooth2 = spline(x_smooth2)

# plt.figure(figsize=(10, 6))
# plt.plot(x_smooth2, y_smooth2, label='Smoothed Energy Loss 2', color='green')

# plt.legend()
# plt.tight_layout()
# plt.show()
