import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.stats import binned_statistic
import os

def sort_nums(filename):
    return int(filename.split('_')[-1].split('.')[0])
    
# path = "/Volumes/Simulations/emissivity_warped_5_10"
# directory_path = os.path.join(os.path.expanduser("~"), path)
# unsorted_files = [os.path.join(directory_path, name) for name in os.listdir(directory_path)]
# filtered_files = [file for file in unsorted_files if not os.path.basename(file).startswith("._")]
# organized_files = sorted(filtered_files, key=sort_nums)

image_filename = 'emissivity2d.fits'#organized_files[0]
data = np.loadtxt(image_filename)

radius = data[:, 0]
phi = data[:, 1]
emissivity = data[:,5]

radius = radius.reshape(50,50)
phi = phi.reshape(50,50) #np.array([r1 + r2 + r3 + r4 for r1, r2, r3, r4 in zip(phi1, phi2, phi3, phi4)])
emissivity = emissivity.reshape(50, 50)#np.array([r1 + r2 + r3 + r4 for r1, r2, r3, r4 in zip(emissivity1, emissivity2, emissivity3, emissivity4)]).reshape(50,50)

from matplotlib.animation import FuncAnimation

# Assuming you have already defined 'radius', 'phi', and 'emissivity'

# Create a figure and axis
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlabel('Radius')
ax.set_ylabel('Emissivity')
ax.set_title('Emissivity vs. Radius')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-3, 1e6)
ax.grid(True)

# Initialize the line to be updated in the animation
line, = ax.plot([], [], label='Phi = 0', lw=2)

def animate(i):
    ax.collections.clear()  # Clear the previous plot
    line.set_data(radius[:, i], emissivity[:, i])
    line.set_label(f'Phi = {np.round(phi[0][i], 3)} = {np.round(np.rad2deg(phi[0][i]), 3)} degrees')
    ax.legend()

# Create the animation
ani = FuncAnimation(fig, animate, frames=range(50), repeat=False)

plt.show()