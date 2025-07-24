import numpy as np
import matplotlib.pyplot as plt

def fn_delta150(x):
    y = np.full_like(x, 150.0)   # Create an array of 150.0 with the same shape as x
    return y

def deltaS1Y_correction_function(x, *par):
    y = par[0] + par[1]*x
    return y

def deltaS1Y_correction_factor(x, *par):
    dy = fn_delta150(x) - deltaS1Y_correction_function(x, *par)
    return dy

paramsDeltaS1Y = (149.405, 0.0957561)
x = np.linspace(-150, 150, 1000)
y = deltaS1Y_correction_function(x, *paramsDeltaS1Y)
# dy = deltaS1Y_correction_factor(x, *paramsDeltaS1Y)

plt.plot(x, deltaS1Y_correction_factor(x, *paramsDeltaS1Y), label="Correction Factor")
plt.plot(x, fn_delta150(x), label='y = 150.0')
plt.plot(x, y, label='deltaCa_st @S1Y correction function')
plt.title('Delta @S1Y correction')
plt.legend()
plt.grid()
plt.show()