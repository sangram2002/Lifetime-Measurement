
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import erf
from math import sqrt

from matplotlib.ticker import AutoMinorLocator

# Define the function to fit
def my_function(params, x, y_data, errors):
    y0, A, t0, w = params
    z = ((x - xc) / w) - w / t0
    arg = 0.5 * (w / t0)**2 - (x - xc) / t0
    erf_arg = z / sqrt(2) + 1
    y_fit = y0 + (A / t0) * np.exp(arg) * (erf(erf_arg) / 2)
    chi_sq = np.sum(((y_fit - y_data) / errors)**2)
    return chi_sq

# Load data from the file
data = np.loadtxt('133Cs_alldata.dat')
x_data = data[1985:2203, 0]
y_data = data[1985:2203, 1]
errors = data[1985:2203, 2]
xc = 2005.99999999999999
# Initial guess for parameters y0 A t0 w
initial_params = [-1.0, 100000, 15, 40]

# Perform chi-square minimization using the minimize function
result = minimize(my_function, initial_params, args=(x_data, y_data, errors))

# Extract the optimized parameters and covariance matrix
optimal_params = result.x
cov_matrix = result.hess_inv

# Calculate errors from the diagonal elements of the covariance matrix
errors_params = np.sqrt(np.diag(cov_matrix))

# Generate fitted data using the optimized parameters
z_opt = ((x_data - xc) / optimal_params[3]) - optimal_params[3] / optimal_params[2]
arg_opt = 0.5 * (optimal_params[3] / optimal_params[2])**2 - (x_data - xc) / optimal_params[2]
erf_arg_opt = z_opt / sqrt(2) + 1
y_fit = optimal_params[0] + (optimal_params[1] / optimal_params[2]) * np.exp(arg_opt) * (erf(erf_arg_opt) / 2)

# Plot the data points, fitted curve, and error bars
plt.errorbar(x_data, y_data, yerr=errors, fmt='o', label='Data with Errors', markersize=3, color='blue', linewidth=0.8)
plt.plot(x_data, y_fit, label=f'Fitted Curve: $y = y_0 + \\frac{{A}}{{t_0}} \\exp\\left(0.5\\left(\\frac{{w}}{{t_0}}\\right)^2 - \\frac{{x - x_c}}{{t_0}}\\right) \\left(\\frac{{\\text{{erf}}\\left(\\frac{{z}}{{\\sqrt{{2}}}} + 1\\right)}}{{2}}\\right)$', color='red')

# Add minor ticks on both x and y axes
plt.minorticks_on()

# Set the minor locator for both x and y axes
plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
plt.gca().yaxis.set_minor_locator(AutoMinorLocator())

# Display the optimized parameters and their errors
plt.text(0.95, 0.75, f'y0: {optimal_params[0]:.3f} ± {errors_params[0]:.3f}', transform=plt.gca().transAxes, ha='right', va='center', color='black', fontsize=10)
plt.text(0.95, 0.70, f'A: {optimal_params[1]:.3e} ± {errors_params[1]:.3e}', transform=plt.gca().transAxes, ha='right', va='center', color='black', fontsize=10)
plt.text(0.95, 0.65, f't0: {optimal_params[2]:.3f} ± {errors_params[2]:.3f}', transform=plt.gca().transAxes, ha='right', va='center', color='black', fontsize=10)
plt.text(0.95, 0.60, f'w: {optimal_params[3]:.3f} ± {errors_params[3]:.3f}', transform=plt.gca().transAxes, ha='right', va='center', color='black', fontsize=10)

# Display the expression for z
plt.text(0.25, 0.75, f'where $z = \\frac{{(x - x_c)}}{{w}} - \\frac{{w}}{{t_0}}$', transform=plt.gca().transAxes, ha='left', va='center', color='black', fontsize=10)

plt.xlabel('Channel Number  (1 Ch=200 ps)')
plt.ylabel('Number of Counts')
plt.title("Fitting the time difference spectrum with the convolution function \n associated with the single exponential decay")
plt.legend()


# Save the plot as a JPG file
plt.savefig('1.pdf', format='pdf')

plt.show()

# Display the optimized parameters
print("Optimized Parameters:")
print(f"y0: {optimal_params[0]} ± {errors_params[0]}")
print(f"A: {optimal_params[1]} ± {errors_params[1]}")
print(f"t0: {optimal_params[2]} ± {errors_params[2]}")
