import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#fldr='01_run_regular_dns/z_reconstruct'
#fldr='02_run_dns_df'
#fldr='03_run_dns_dfh/z_combined'
#fldr='04_run_1'
#fldr='05_increase_Re'
fldr='06_increase_Re_h'
fname = os.path.join(fldr, 'ubar.dat')
# Step 1: Load the data from the file
data = np.loadtxt(fname)

# Step 2: Extract the two columns (x and y)
x_data = data[:, 0]
y_data = data[:, 1]

print(f'nsteps = {len(x_data):d}')

# Step 3: Define the exponential function
def exponential(x, a, b, c):
    return c - a * np.exp(-b * x)

# Step 4: Perform the exponential fit
params0, covariance0 = curve_fit(exponential, x_data, y_data)

# Step 5: Extract the fitted parameters
a0, b0, c0 = params0
print(f"Fitted parameters: a = {a0}, b = {b0}, c = {c0}")

# Step 6: Generate the fitted curve
y_fit0 = exponential(x_data, *params0)

# Step 4: Perform the exponential fit
step = 10
srange = list(range(10,1001,10)) + list(range(1100, 10001, 100)) + list(range(11000, 40001, 1000))
av, bv, cv = [], [], []
for i in srange:
   params, covariance = curve_fit(exponential, x_data[:i], y_data[:i])

   # Step 5: Extract the fitted parameters
   a, b, c = params
   av.append(a)
   bv.append(b)
   cv.append(c)
   print(f"Fitted parameters: a = {a}, b = {b}, c = {c}")

plt.figure()
plt.title(' y = c - a * e^(-b * x)')
plt.plot(srange, abs(av-a0)/a0, label='a')
plt.plot(srange, abs(bv-b0)/b0, label='b')
plt.plot(srange, abs(cv-c0)/c0, label='c')
ax = plt.gca()
ax.set_yscale('log')
plt.legend()

# Step 6: Generate the fitted curve
y_fit = exponential(x_data, *params0)

# Step 7: Plot the data and the fit
plt.figure()
plt.scatter(x_data, y_data, label='Data', color='black')
plt.plot(x_data, y_fit0, label=f'Exponential Fit: y = {c0:.8f} - {a0:.8f} * e^(-{b0:.8f} * x)', color='red')
plt.axhline(c0)
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()