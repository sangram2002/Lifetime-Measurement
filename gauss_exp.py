import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

def custom_function(x, y0, A, t0, w, xc):
    z = ((x - xc) / w) - w / t0
    exp_term = np.exp(0.5 * np.clip((w / t0) ** 2 - (x - xc) / t0, -500, 500))
    y = y0 + A / t0 * exp_term * (erf(z / np.sqrt(2)) + 1) / 2
    return y

# Load data from the file
data = np.loadtxt('133Cs_timing.dat')
x_data = data[:, 0]
y_data = data[:, 1]
errors = data[:, 2]

# Set parameter values
y0 = 95.83
A = 52169.65
t0 = 5.5e-13
w = 22.4
xc = 2019.0

# Generate x values
x_values = x_data

# Calculate corresponding y values
y_values = custom_function(x_values, y0, A, t0, w, xc)

# Plot the function
plt.plot(x_values, y_values, label='y = y0 + A/t0 * exp(0.5*(w/t0)^2-(x-xc)/t0)*(erf(z / (sqrt(2)) +1)/2')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
