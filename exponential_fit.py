import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Define the function to fit
def my_function(params, x):
    y0, A, t0 = params
    return y0 + A * np.exp(-x / t0)

# Define the chisquare function
def chi_square(params, x, y_data, errors):
    y_fit = my_function(params, x)
    chi_sq = np.sum(((y_fit - y_data) / errors)**2)
    return chi_sq

# Numerical approximation of the Hessian matrix
def numerical_hessian(func, params, args, epsilon=1e-5):
    hessian_matrix = np.zeros((len(params), len(params)))
    for i in range(len(params)):
        for j in range(len(params)):
            params_plus_plus = np.copy(params)
            params_plus_plus[i] += epsilon
            params_plus_plus[j] += epsilon
            chi_plus_plus = func(params_plus_plus, *args)

            params_minus_plus = np.copy(params)
            params_minus_plus[i] -= epsilon
            params_minus_plus[j] += epsilon
            chi_minus_plus = func(params_minus_plus, *args)

            params_minus_minus = np.copy(params)
            params_minus_minus[i] -= epsilon
            params_minus_minus[j] -= epsilon
            chi_minus_minus = func(params_minus_minus, *args)

            params_minus_minus[j] += 2 * epsilon

            hessian_matrix[i, j] = (chi_plus_plus - chi_minus_plus - chi_minus_minus + func(params_minus_minus, *args)) / (4 * epsilon**2)

    return hessian_matrix

# Load data from the file
data = np.loadtxt('133Cs_only_exp.txt')  # Replace 'your_data_file.txt' with your actual data file
x_data = data[:, 0]
y_data = data[:, 1]
errors = data[:, 2]  # Assuming the 3rd column contains errors

# Initial guess for parameters
initial_params = [-2.69, 1.0e22, 40.0]

# Perform chisquare minimization using the minimize function
result = minimize(chi_square, initial_params, args=(x_data, y_data, errors), method='L-BFGS-B')

# Extract the optimized parameters
y0_opt, A_opt, t0_opt = result.x

# Compute the Hessian matrix numerically
hessian_matrix = numerical_hessian(chi_square, result.x, args=(x_data, y_data, errors))

# Calculate errors from the diagonal elements of the covariance matrix
errors_params = np.sqrt(np.diag(np.linalg.inv(hessian_matrix)))

# Generate fitted data using the optimized parameters
y_fit = my_function([y0_opt, A_opt, t0_opt], x_data)

# Plot the data points, fitted curve, and error bars
plt.errorbar(x_data, y_data, yerr=errors, fmt='o', label='Data with Errors', markersize=3, color='blue', linewidth=0.5)
plt.plot(x_data, y_fit, label=f'Fitted Curve: $y_0 + A \cdot e^{{-x / t_0}}$', color='red')

# Display the optimized parameters and their errors on the plot
plt.text(0.95, 0.8, f'y0: {y0_opt:.3f} ± {errors_params[0]:.3f}', transform=plt.gca().transAxes, ha='right', va='center', color='black', fontsize=10)
plt.text(0.95, 0.75, f'A: {A_opt:.3e} ± {errors_params[1]:.3e}', transform=plt.gca().transAxes, ha='right', va='center', color='black', fontsize=10)
plt.text(0.95, 0.7, f't0: {t0_opt:.3f} ± {errors_params[2]:.3f}', transform=plt.gca().transAxes, ha='right', va='center', color='black', fontsize=10)

# Customize plot settings
plt.xlabel('Channel Number    (1 Ch=200 ps) ')
plt.ylabel('Number of Counts')
plt.title('Fitting the Convolution of a Prompt with an Exponential Function \n to find lifetime of the desired state')
plt.legend()

# Save the plot as a JPG file
plt.savefig('3.pdf', format='pdf')


# Display the optimized parameters and their errors
print("Optimized Parameters:")
print(f"y0: {y0_opt} ± {errors_params[0]}")
print(f"A: {A_opt} ± {errors_params[1]}")
print(f"t0: {t0_opt} ± {errors_params[2]}")

# Show the plot
plt.show()


