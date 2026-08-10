import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 14,
    "font.size": 14,
    "legend.fontsize": 12,
})

def load_time_energy(filename):
    """
    Loads data of the form:
    time   <value>   energy   <value>
    and returns numpy arrays time, energy
    """
    times = []
    energies = []
    with open(filename) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            # expected: ["time", t, "energy", e]
            t = float(parts[1])
            e = float(parts[3])
            times.append(t)
            energies.append(e)
    return np.array(times), np.array(energies)

# Load both files
t1, e1 = load_time_energy("energy_3D.txt")
t2, e2 = load_time_energy("energy_2Dh.txt")

# Normalize energies
e1_norm = e1 / e1[0]
e2_norm = e2 / e2[0]

# --------------------------------------
# Plot 1: semilogy of file 1 normalized
# --------------------------------------
plt.figure(figsize=(4,4))
plt.semilogy(t1, e1_norm, label="3D")
plt.xlabel(r"$t$")
plt.ylabel(r"$E / E_0$")
plt.grid(True)
plt.legend()
plt.xlim([0, 15])
plt.ylim([0.015, 2])
plt.tight_layout()

# ---------------------------------------------------------
# Plot 2: same as Plot 1 + scatter of file 2 every 100 steps
# ---------------------------------------------------------
plt.figure(figsize=(4,4))
plt.semilogy(t1, e1_norm, label="3D")

# Scatter file 2 every 100 points
step = 10
plt.scatter(t2[::step], e2_norm[::step], color='red', s=15, label="2.5D")

plt.xlabel(r"$t$")
plt.ylabel(r"$E / E_0$")
plt.grid(True)
plt.legend()
plt.xlim([0, 15])
plt.ylim([0.015, 2])
plt.tight_layout()
plt.show()