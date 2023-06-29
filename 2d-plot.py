import sys
sys.path.append("/local/jweatson/i2e_plotters")
from plot2d import Image, plot_circles
from i2_plot import Simulation
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import os
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "serif",
  "font.serif": ["Computer Modern Roman"],
})

working_dir = "./"
filename = working_dir+sys.argv[1]
image = Image(filename)

# norm = Normalize(vmin=3000, vmax=6000)

f, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [3, 1, 1]})

f.set_figheight(10)
f.set_figwidth(6)
pcm = a0.imshow(image.heat.data,cmap="inferno")
# Plot circles for mantle and core
a0.set_xlabel("X (km)")
a0.set_ylabel("Y (km)")
# Build the colour bar
cbar = f.colorbar(pcm,ax=a0,label=image.temp.label)

figname = "figures/fig_"+str(000)+".png"
plt.savefig(figname,dpi=600)
print("Finished "+filename+" -> "+figname)