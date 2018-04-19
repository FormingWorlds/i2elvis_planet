# -*- coding: utf-8 -*-
"""
Visualisation routine for I2ELVIS_planet .prn output files

Created on Thu Feb 26 17:05:18 2015 @author: Tim Lichtenberg
Edit 2015-05-15: Change contour to heatmap, various layout changes

based on the Matlab 2D plotting routine by Gregor Golabek
"""

# Import packages
import matplotlib # Plotting routine
matplotlib.use('Agg') # to counteract problems with torque usage
import os #https://stackoverflow.com/questions/3964681/find-all-files-in-directory-with-extension-txt-with-python
import sys
import glob
import struct # Decode binary data
import numpy as np 
from natsort import natsorted #https://pypi.python.org/pypi/natsort

# Packages for plotting
#import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, DrawingArea, HPacker
from matplotlib import rc
from matplotlib import cm
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import colorlover as cl #http://nbviewer.ipython.org/gist/jackparmer/4696226c9b6b931bbdf6
from IPython.display import HTML
import locale # http://stackoverflow.com/questions/1823058/how-to-print-number-with-commas-as-thousands-separators
import datetime

# https://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib
params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
plt.rcParams.update(params)

###########################
###### PLOT SETTINGS ######
###########################

# Define path to run folder
working_dir = os.getcwd()

# Define path to folder with .prn files, examples:
PDIR = "/Volumes/astro/meyer/timli/19_55cnc"
PDIR = working_dir+"/../tests"

# Define whicploth RUNS you want to plot
# RUNS = [ "test1[8,9]*", "test2[0,1,2]*" ]
RUNS = [ "test*" ]

# Define which FILES you want to plot
FILES = [ '*_*.prn' ] # '*_*[8,9]??.prn'

# # Define path where to store output images
# image_dir = "/Users/tim/Dropbox/work/Projects/18_pebble/figures/"

###########################
### Initial definitions ###
###########################

al_canonical = 5.25e-5
fe_canonical = 1.15e-8

# Choose way to plot for runs and files
run_mode = 'glob' # options: 'glob', 'specific'
file_mode = 'glob' # options: 'glob', 'specific'
# print "Run mode: ", run_mode
# print "File mode: ", file_mode

### Define what kind of output to produce
# Specify list of node quantities to plot
# Available node quantities: "eta", "rho", "tmp", "press", "kt"
plot_nodes    = [ "rho", "tmp" ] #"eta , "prs"
# Specifiy list of marker quantities to plot
# Available marker quantities: "comp", "por"
plot_markers  = [ "comp", "tmax", "acc"  ] # "por", "comp", "tmax", "acc" 

# Shall isotherm show the surface?
isotherm      =  0
# Which temperature shall be used for isotherm plot?
iso_temp      =  290
# Which depths [m] shall be used for isolines plot?
depth_out     =  1.1e4
depth_in      =  1.8e4
# Which peridotite solidus temperature? 
# 0: Hirschmann (2000), 1: Herzberg et al. (2000)
silicate_melt =  1
# Which iron melting temperature? 
# 0: Pure iron based on Boehler (1993), 1: Fe-FeS eutectic based on 
# data compilation by Chudinovskikh & Boehler (2007)
iron_melt     =  1

# Set run list to analyse
if run_mode == 'glob':
    os.chdir(PDIR)
    RUN_LIST = []
    for i in range(0,len(RUNS)):
        glob_runs=glob.glob(RUNS[i])
        for item in glob_runs:
            RUN_LIST.append(item)
if run_mode == 'specific':
    RUN_LIST = RUNS

# sort run list
RUN_LIST=natsorted(RUN_LIST)

print "List of runs to analyse:", RUN_LIST

##########################
### Loop over all runs ###
##########################

for RUN in RUN_LIST:

    # Inform about current run
    print "##########", "Analyse", RUN, "(", PDIR+"/"+RUN, ")", "##########"
    
    # Start the main plotting routine

    # print RUN

    # FILES = [ '*.prn' ] 

    # Loop over all chosen datafiles in directory
    if file_mode == 'glob':
        os.chdir(PDIR+"/"+RUN)
        FILE_LIST = []
        for i in range(0,len(FILES)):
            glob_files=glob.glob(FILES[i])
            for item in glob_files:
                FILE_LIST.append(item)
    if file_mode == 'specific':
        FILE_LIST = FILES

    # sort file list
    FILE_LIST=natsorted(FILE_LIST)

    print "List of files to analyse:", FILE_LIST

    for FILE in FILE_LIST:
        
         # Read in file
        # http://www.tutorialspoint.com/python/python_files_io.htm
        # https://docs.python.org/2/library/struct.html#struct-alignment
        # http://www.tutorialspoint.com/python/python_tuples.htm
        # https://stackoverflow.com/questions/8092469/reading-a-binary-file-in-python
        fdata = open(PDIR+"/"+RUN+"/"+FILE, "rb")
        print datetime.datetime.now(), fdata.name        

        # Read sizes of variables
        # A = struct.unpack("BBBB",fdata.read(4))[0]
        A           = fdata.read(4)
        # Read model parameters
        # Grid resolution
        xnumx       = struct.unpack("q",fdata.read(8))[0]
        ynumy       = struct.unpack("q",fdata.read(8))[0]
        # Markers per cell
        mnumx       = struct.unpack("q",fdata.read(8))[0]
        mnumy       = struct.unpack("q",fdata.read(8))[0]
        # Number of markers
        marknum     = struct.unpack("q",fdata.read(8))[0]
        # Model sizes
        xsize       = struct.unpack("d",fdata.read(8))[0]
        ysize       = struct.unpack("d",fdata.read(8))[0]
        # Efficieny of impact heating 
        gamma_eff   = struct.unpack("d",fdata.read(8))[0]
        # Thermal memory of iron material
        memory_fe   = struct.unpack("d",fdata.read(8))[0]
        # Thermal memory of silicate material
        memory_si   = struct.unpack("d",fdata.read(8))[0]
        # Porosity parameter
        por_init    = struct.unpack("d",fdata.read(8))[0]
        # Porosity parameter
        growth_model = struct.unpack("q",fdata.read(8))[0]
        # Grain size parameter
        gr_init     = struct.unpack("d",fdata.read(8))[0]
        # Grid resolution in z direction
        znumz       = struct.unpack("q",fdata.read(8))[0]
        # 2D cylindrical or 3D spherical
        corr2d3d    = struct.unpack("d",fdata.read(8))[0]
        # Pressure value
        pinit       = fdata.read(40)
        # Gravity
        GXKOEF      = struct.unpack("d",fdata.read(8))[0]
        GYKOEF      = struct.unpack("d",fdata.read(8))[0]
        # Ambient temperature for sticky air
        tmp_ambient = struct.unpack("d",fdata.read(8))[0]
        # 26Al/27Al ratio
        al2627_init = struct.unpack("d",fdata.read(8))[0]
        # 60Fe/56Fe ratio
        fe6056_init = struct.unpack("d",fdata.read(8))[0]
        # Number of rocks
        rocknum     = struct.unpack("i",fdata.read(4))[0]
        # Number of Boundary conditions
        bondnum     = struct.unpack("q",fdata.read(8))[0]
        # Stage
        n1          = struct.unpack("i",fdata.read(4))[0]
        # Time
        Time        = struct.unpack("d",fdata.read(8))[0]
        
        # print A, xnumx, ynumy, mnumx, mnumy, marknum
        # print xsize, ysize, gamma_eff, memory_fe, memory_si
        # print por_init, znumz, corr2d3d, pinit, GXKOEF
        # print GYKOEF, rocknum, bondnum, n1, Time

        # Define cursor positions within the binary file in bytes
        # 1*int (A), x*int, y*float/long
        curpos_nodes    = 4 + 2*4 + 27*8 + rocknum*(8*24+4)
        # Skip all nodes
        curpos_grid     = curpos_nodes+(4*20+8*4)*xnumx*ynumy
        # Skip gridlines and boundary conditions
        curpos_markers  = curpos_grid+4*(xnumx+ynumy)+(4*4+8*3)*(bondnum-1)

        print "Read in node properties..."

        # Skip rock properties
        #curpos0=4+2*4+22*8+rocknum*(8*24+4)
        fdata.seek(curpos_nodes)

        prs   = np.zeros((ynumy,xnumx))  # pressure
        vx    = np.zeros((ynumy,xnumx))  # x component of the velocity
        vy    = np.zeros((ynumy,xnumx))  # y component of the velocity
        exx   = np.zeros((ynumy,xnumx))  # xx component of strain rate tensor
        eyy   = np.zeros((ynumy,xnumx))  # yy component of strain rate tensor
        exy   = np.zeros((ynumy,xnumx))  # xy component of strain rate tensor
        sxx   = np.zeros((ynumy,xnumx))  # xx component of stress tensor
        syy   = np.zeros((ynumy,xnumx))  # yy component of stress tensor
        sxy   = np.zeros((ynumy,xnumx))  # xy component of stress tensor
        rho   = np.zeros((ynumy,xnumx))  # density
        eta   = np.zeros((ynumy,xnumx))
        etd   = np.zeros((ynumy,xnumx))
        mu    = np.zeros((ynumy,xnumx))
        gp    = np.zeros((ynumy,xnumx))  # gravity potential
        ep    = np.zeros((ynumy,xnumx))  # viscosity 1
        et    = np.zeros((ynumy,xnumx))  # viscosity 2
        tmp   = np.zeros((ynumy,xnumx))  # temperature
        cp    = np.zeros((ynumy,xnumx))  # heat capacity
        kt    = np.zeros((ynumy,xnumx))  # thermal conductivity
        ht    = np.zeros((ynumy,xnumx))  # radiogenic heating
        eii   = np.zeros((ynumy,xnumx))  # strain rate tensor
        sii   = np.zeros((ynumy,xnumx))  # stress tensor
        dis   = np.zeros((ynumy,xnumx))  # dissipation
        fm    = np.zeros((ynumy,xnumx))  # silicate melt fraction
        fm_fe = np.zeros((ynumy,xnumx))  # iron melt fraction
        radius= np.zeros((ynumy,xnumx))  # distance from planetary center [m]

        # Read nodes information
        for i in range(0,xnumx):
            for j in range(0,ynumy):
                prs[j,i]=struct.unpack("f",fdata.read(4))[0]
                vx[j,i]=struct.unpack("f",fdata.read(4))[0]
                vy[j,i]=struct.unpack("f",fdata.read(4))[0]
                vbuf1=fdata.read(3*8)
                exx[j,i]=struct.unpack("f",fdata.read(4))[0]
                eyy[j,i]=struct.unpack("f",fdata.read(4))[0]
                exy[j,i]=struct.unpack("f",fdata.read(4))[0]
                sxx[j,i]=struct.unpack("f",fdata.read(4))[0]
                syy[j,i]=struct.unpack("f",fdata.read(4))[0]
                sxy[j,i]=struct.unpack("f",fdata.read(4))[0]
                rho[j,i]=struct.unpack("f",fdata.read(4))[0]
                eta[j,i]=struct.unpack("f",fdata.read(4))[0]
                etd[j,i]=struct.unpack("f",fdata.read(4))[0]
                mu[j,i]=struct.unpack("f",fdata.read(4))[0]
                gp[j,i]=struct.unpack("f",fdata.read(4))[0]
                ep[j,i]=struct.unpack("f",fdata.read(4))[0]
                et[j,i]=struct.unpack("f",fdata.read(4))[0]
                tmp[j,i]=struct.unpack("f",fdata.read(4))[0]
                vbuf3=fdata.read(1*8)
                cp[j,i]=struct.unpack("f",fdata.read(4))[0]
                kt[j,i]=struct.unpack("f",fdata.read(4))[0]
                ht[j,i]=struct.unpack("f",fdata.read(4))[0]

        prs = prs*1.000e-5;  # Convert pressure from Pascal to Bar
        
        for i in range(0,xnumx):
            for j in range(0,ynumy):
                
                # Not at the boundaries of the FD grid                    
                if(j <= (ynumy-2) and i <= (xnumx-2)):

                    if(silicate_melt == 0):

                        if (((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) < 10e+4):
                            # Take peridotite solidus temperature [K] formulation from Hirschmann (2000)                        
                            Tsol = 273.15+1120.661+132.899e-4*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)-5.104e-8*np.power(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4),2)
                        else:
                            # For high pressures take peridotite solidus [K] temperature extrapolation based on Hirschmann (2000)                     
                            Tsol = 273.15+1939.251+30.819e-4*(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)-10e+4)

                    # Take peridotite solidus temperature [K] formulation from Herzberg et al. (2000)
                    elif(silicate_melt == 1):
                        if ((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) < 21.50e+4:
                            Tsol = 273.15+1143.04342+58.2946423e-4*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)+52.3439318e-8*np.power(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4),2)-16.3201032e-12*np.power(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4),3)+2.29886314e-16*np.power(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4),4)-0.180865486e-20*np.power(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4),5)+0.00815679773e-24*np.power(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4),6)-0.000197104325e-28*np.power(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4),7)+1.97908526e-38*np.power(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4),8)

                        else:
                            # For high pressures take peridotite solidus temperature [K] extrapolation based on Herzberg et al. (2000)
                            Tsol = 273.15+2157.500+11.7297e-4*(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)-21.50e+4)

    
                    # Liquidus temperature [K] for peridotites taken from Wade & Wood (2005)                     
                    Tliq   = 1973.0+28.57e-4*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)

                    # Iron melting temperature [K] using parametrization by Boehler (1993)           
                    if (iron_melt == 0):
                        if (((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) < 0.100e6):
                            Tmelt = 1761.00+3.100e-3*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)
                        elif ((((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) >= 0.100e6) and (((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) < 0.200e6)):
                            Tmelt = 1863.00+2.080e-3*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)                           

                        elif ((((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) >= 0.200e6) and (((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) < 0.600e6)):
                            Tmelt = 2071.80+1.035e-3*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)                           

                        elif ((((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) >= 0.600e6) and (((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) < 1.000e6)):
                            Tmelt = 2382.80+5.170e-4*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)                            
                    # Fe-FeS eutectic melting temperature [K] using data compilation by Chudinovskikh & Boehler (2007)                
                    elif (iron_melt == 1):               
                        if (((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) <= 0.418e6):
                            Tmelt = 1260.10+3.3171e-4*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)-3.395e-9*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)**2+2.660e-14*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)**3-3.7688e-20*((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)**4
                        elif (((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4) > 0.418e6):
                            Tmelt = 1597.70+4.2635e-4*(((prs[j,i]+prs[j+1,i]+prs[j,i+1]+prs[j+1,i+1])/4)-0.418e6)
                            
                # S I L I C A T E  M E L T I N G
                Tsol = min((Tliq-100),Tsol)
                melt = (tmp[j,i]-Tsol)/(Tliq-Tsol)

                # Make sure melt amount is within range 0 to 100 %
                if melt >= 1.000:
                    melt = 1.000

                if melt <= 0.000:
                    melt = 0.000

                # Discriminate between iron and silicates using their sigificantly
                # different densities
                if rho[j,i] <= 6000.000:             # Material is silicate melt
                    fm[j,i] = min(1,max(0,melt))
                elif rho[j,i] <= 1000.000:           # Material is most probably sticky air 
                    fm[j,i] = -0.001                 # Set to very small negative value (only for plotting purposes)
                elif rho[j,i] > 6000.000:            # Material is most probably iron
                    fm[j,i] = -0.001                 # Set to very small negative value (only for plotting purposes)
                    
                # Sticky air is cold, good way to remove unphysical "melt" in sticky air                   
                if tmp[j,i] <= iso_temp:
                    fm[j,i] = -0.001                 # Set to very small negative value (only for plotting purposes)
           
                # I R O N  M E L T I N G                    
                if tmp[j,i] < Tmelt:
                    melt_fe = 0.000
                if tmp[j,i] >= Tmelt:
                    melt_fe = 1.000

                # Make sure iron melt amount is within range 0 to 100 %
                if melt_fe >= 1.000:
                    melt_fe = 1.000                     
                if melt_fe <= 0.000:
                    melt_fe = 0.000
                    
                # Discriminate between iron and silicates using their sigificantly different densities
                if rho[j,i] <= 7000.000:              # Material is silicate melt or sticky air material
                    fm_fe[j,i] = -0.001             # Set to very small negative value (only for plotting purposes)
                elif rho[j,i] > 7000.000:             # Material is most probably iron
                    fm_fe[j,i] = min(1,max(0,melt_fe))

        # Computing shear heating
        for i in range(0,xnumx-3):
            for j in range(0,ynumy-3):
                dis[j,i]=((sxx[j,i])*exx[j,i]+(syy[j,i])*eyy[j,i]+2.*(sxy[j,i]*exy[j,i]))*1e+06
        # Dissipation is not plotted if it takes place in sticky air or stabilizing material
                if rho[j,i]<= 1000.000:
                    dis[j,i] = 0.000000
                if dis[j,i] <= 0.000000:
                    dis[j,i] = 0.000000

        # Compute distance of each grid point from planetary center [m]
        for i in range(1,xnumx):
            for j in range(1,ynumy):    
                radius[j,i]=(((i-1)*(xsize/(xnumx-1))-(xsize/2.000))**(2.000)+((j-1)*(ysize/(ynumy-1))-(ysize/2.000))**(2.000))**(0.500)

        # Skip all nodes
        fdata.seek(curpos_grid)
        # Read Gridline positions
        gx = np.zeros(xnumx)
        gy = np.zeros(ynumy)
        for i in range(0,xnumx):
            gx[i]=struct.unpack("f",fdata.read(4))[0]
        for i in range(0,ynumy):
            gy[i]=struct.unpack("f",fdata.read(4))[0]
        # end: if plot nodes

        if len(plot_markers) > 0:
            
            print "Read in marker properties..."

            # Skip rock properties, nodes, gridlines and boundary conditions
            fdata.seek(curpos_markers)

            #print 'Visualize marker type and porosity'

            # Marker grid resolution
            xresol=int(xnumx)#*4 #for seemingly higher resolution
            yresol=int(ynumy)#*4
            xstp=xsize/(xresol-1)
            ystp=ysize/(yresol-1)
            mx=np.arange(0,xsize/1000,xstp/1000)
            my=np.arange(0,ysize/1000,ystp/1000)
            # Read markers
            markcom=np.empty((yresol,xresol))
            por_arr=np.empty((yresol,xresol))
            tmax_arr=np.empty((yresol,xresol))
            acc_arr=np.empty((yresol,xresol))
            markcom[:]=np.NAN
            por_arr[:]=np.NAN
            tmax_arr[:]=np.NAN
            acc_arr[:]=np.NAN
            markdis=1e+20*np.ones((yresol,xresol))
            ynplot=0
            nplot=0

            for mm1 in range(1,marknum+1):
                # Read data for current marker
                mbuf=np.zeros(9)
                for i in range(0,9):
                    mbuf[i]=struct.unpack("f",fdata.read(4))[0]
                markmg_old=fdata.read(1)
                markmg_time=fdata.read(1)
                markt=struct.unpack("B",fdata.read(1))[0]

                # Immobile markers
                if (markt>=100):
                    markt=markt-100
                # Rock markers
                if (markt<50):
                    # Read x and z positions of markers
                    markx=mbuf[0]
                    marky=mbuf[1]
                    # Read marker porosity
                    markpor=mbuf[5]
                    # Read marker grain size
                    markgr=mbuf[6]
                    # Read marker maximum temperature
                    marktmax=mbuf[7]
                    # Read marker accretion time
                    markacc=mbuf[8]

                    # Mask sticky air layers
                    if markt == 0:
                        markacc     = np.nan
                        marktmax    = np.nan
                        markpor     = np.nan
                        markgr      = np.nan

                    # Define visualization cell
                    m1=float(int(markx/xstp-0.5))
                    m2=float(int(marky/ystp-0.5))


                    # Define special cases
                    if (m1<1):
                        m1=1
                    if (m1>xresol-1):
                        m1=xresol-1
                    if (m2<1):
                        m2=1
                    if (m2>yresol-1):
                        m2=yresol-1

                    m1 = int(m1)
                    m2 = int(m2)

                    # print m2, m1

                    # Define distance for 4 surrounding nodes
                    # 1 3
                    # 2 4
                    # Node 1
                    dx=markx-(m1-1)*xstp
                    dy=marky-(m2-1)*ystp
                    dd=(dx*dx+dy*dy)**0.5
                    if (dd<markdis[m2,m1]):
                        markcom[m2,m1]      =   markt
                        por_arr[m2,m1]      =   markpor
                        tmax_arr[m2,m1]     =   marktmax
                        acc_arr[m2,m1]      =   markacc
                        markdis[m2,m1]      =   dd
                    # Node 2
                    dy=dy-ystp;
                    dd=(dx*dx+dy*dy)**0.5
                    if (dd<markdis[m2+1,m1]):
                        markcom[m2+1,m1]    =   markt
                        por_arr[m2+1,m1]    =   markpor
                        tmax_arr[m2+1,m1]   =   marktmax
                        acc_arr[m2+1,m1]    =   markacc
                        markdis[m2+1,m1]    =   dd
                    # Node 4
                    dx=dx-xstp;
                    dd=(dx*dx+dy*dy)**0.5
                    if (dd<markdis[m2+1,m1+1]):
                        markcom[m2+1,m1+1]  =   markt
                        por_arr[m2+1,m1+1]  =   markpor
                        tmax_arr[m2+1,m1+1] =   marktmax
                        acc_arr[m2+1,m1+1]  =   markacc
                        markdis[m2+1,m1+1]  =   dd
                    # Node 3
                    dy=dy+ystp;
                    dd=(dx*dx+dy*dy)**0.5;
                    if (dd<markdis[m2,m1+1]):
                        markcom[m2,m1+1]    =   markt
                        por_arr[m2,m1+1]    =   markpor
                        tmax_arr[m2,m1+1]   =   marktmax
                        acc_arr[m2,m1+1]    =   markacc
                        markdis[m2,m1+1]    =   dd

        # PLOTTING ROUTINE

        # Find name of current file
        FILE=FILE.split("_")
        FILE=FILE[-1]
        RUN=RUN.split("/")
        RUN=RUN[-1]
        print "Plotting run", RUN, ", file", FILE

        # Get the settings of the current run
        if RUN.startswith("r"):
            r_start     = str(int(RUN[1:4]))
            t_start     = str(int(RUN[5:8])/100)
            al_start    = str(int(RUN[10:15])*1e-3*1e-5/al_canonical)
            fe_start    = str(int(RUN[17:22])*1e-3*1e-8/fe_canonical)

        # Time formatting
        locale.setlocale(locale.LC_ALL, 'en_US') # http://stackoverflow.com/questions/1823058/how-to-print-number-with-commas-as-thousands-separators
        time_Myr        = float(Time)/1e6
        time_kyr        = float(Time)/1e3
        # t_current       = '%.3f' % time_kyr # str(round(Time/1e6,2))
        t_current       = locale.format("%d", round(Time), grouping=True)
        # t_current       = '%.3f' % Time_Myr # str(round(Time/1e6,2))
        # t_current       = str(Time)

        # por_start   = str(float(RUN[5:-5])/100)
        #print r_start, por_start, t_start, t_current
        #print "Done"

        ## sources
        # http://matplotlib.org/users/annotations_guide.html
        # http://matplotlib.org/users/annotations_intro.html
        # https://nickcharlton.net/posts/drawing-animating-shapes-matplotlib.html
        # http://matplotlib.org/users/tight_layout_guide.html
        # http://matplotlib.org/faq/howto_faq.html#howto-subplots-adjust
        # http://matplotlib.org/users/usetex.html
        # https://stackoverflow.com/questions/24997926/making-a-custom-colormap-using-matplotlib-in-python
        # http://matplotlib.org/users/customizing.html#customizing-matplotlib

        #plt.rc('text', usetex=True)
        #plt.rc('font', family='sans-serif')

        # GENERAL DEFINITIONS

        # rcParams settings, see all with mpl.rcParams.keys()
        label_size = 20
        mpl.rcParams['xtick.labelsize'] = label_size
        mpl.rcParams['ytick.labelsize'] = label_size 
        # mpl.rcParams['font.family'] = 'sans-serif'
        #mpl.rcParams['font.sans-serif'] = 'arial'
        #mpl.rcParams['mathtext.fontset'] = 'custom'
        #mpl.rcParams['mathtext.fallback_to_cm'] = True
        #mpl.rcParams['mathtext.default'] = 'sf'
        #mpl.rcParams['mathtext.sf'] = 'arial'
        #mpl.rcParams['mathtext.rm'] = 'arial'


        #mathtext.cal : cursive
        #mathtext.rm  : serif
        #mathtext.tt  : monospace
        #mathtext.it  : serif:italic
        #mathtext.bf  : serif:bold
        #mathtext.sf  : sans
        #mathtext.fontset : cm # Should be 'cm' (Computer Modern), 'stix',
                               # 'stixsans' or 'custom'
        #mathtext.fallback_to_cm : True  # When True, use symbols from the Computer Modern
                                         # fonts when a symbol can not be found in one of
                                         # the custom math fonts.
        #mathtext.default : it # The default font to use for math.
                               # Can be any of the LaTeX font names, including
                               # the special name "regular" for the same font
                               # used in regular text.


        #image_dir="/Users/tim/Dropbox/work/gfd/poros/matplotlib_scripts/"

        #plot_nodes    = [ "rho", "eta", "tmp", "prs", "kt" ]
        #plot_markers  = [ "comp", "por" ]

        # comp_map = "Spectral"
        comp_map = LinearSegmentedColormap.from_list('mycmap', 
        [(0,"#ffffff"),
        (0.17,"#ffffff"),
        (0.19,"#1f78b4"), 
        (0.21,"#ffffff"),
        (0.23,"#a6cee3"), 
        (0.25,"#ffffff"),
        (0.27,"#000000"), 
        (0.29,"#ffffff"),
        (0.31,"#4c4c4c"),
        (0.33,"#ffffff"), 
        (0.95,"#ffffff"),
        (0.96,"#e31a1c"),
        #[0.97,"#ffffff"],
        (0.98,"#ffffff"),
        (1,"#fb9a99")]
         )

        for var in plot_nodes+plot_markers:
            # for specific_map in cmaps_specific: # cmap tests
            
            # image_dir="/net/astrogate/export/astro/meyer/timli/poros2d_viz/" # when running on bluesky
            image_dir=working_dir+"/../figures/"+RUN+"/" # when running on local machine
            
            # Add dir
            if not os.path.exists(image_dir):
                os.makedirs(image_dir)

            z = np.ones((xnumx,ynumy))

            fig = plt.figure()
            ax = plt.gca()
            fig.set_size_inches(9.5, 9.5)
            ax = fig.add_subplot(111)

            # PARAMETER SPECIFIC SETTINGS
            
            if var == "comp":
                var_name=r"Composition"
                var_cmap=comp_map
                data_array=markcom
                unit='Composition code'
                cmap_range_min=0
                cmap_range_max=26
            if var == "por":
                var_name=r"Porosity"
                var_cmap="terrain_r"
                data_array=por_arr
                unit=r'$\phi \, [\mathrm{non-dim.}]$'
                cmap_range_min=0
                cmap_range_max=0.75
            if var == "tmax":
                var_name=r"Maximum temperature"
                var_cmap="magma"
                data_array=tmax_arr
                unit=r'$T_\mathrm{max} \, [\mathrm{K}]$'
                cmap_range_min=np.nanmin(data_array)
                cmap_range_max=np.nanmax(data_array)
            if var == "acc":
                var_name=r"Accretion time"
                var_cmap="Spectral"
                data_array=acc_arr/1e6
                unit=r'$t_\mathrm{acc} \, [\mathrm{Myr}]$'
                cmap_range_min=0.0
                cmap_range_max=np.nanmax(data_array)
            if var == "rho":
                var_name=r"Density"
                var_cmap="Spectral"# "Blues" "gist_heat_r" #hot_r
                data_array=rho
                unit=r'$\rho \, [\mathrm{kg} \, \mathrm{m}^{-2}]$'
                cmap_range_min=2750#int(np.amin(data_array))#3000
                cmap_range_max=int(np.amax(data_array))
            if var == "eta":
                var_name=r"Viscosity"
                var_cmap="terrain_r"
                data_array=eta
                unit=r'$\eta$ [$\mathrm{Pa} \, \mathrm{s}$]'
                cmap_range_min=np.amin(data_array)
                cmap_range_max=np.amax(data_array)
            if var == "tmp":
                var_name=r"Temperature"
                var_cmap="Reds"#"RdBu_r"
                data_array=tmp
                unit=r'$T \, [\mathrm{K}]$'
                cmap_range_min=int(np.amin(data_array))#290
                cmap_range_max=int(np.amax(data_array))
            if var == "prs":
                var_name=r"Pressure"
                var_cmap="seismic"
                data_array=prs
                unit=r'$p \, [\mathrm{bar}]$'
                cmap_range_min=int(np.amin(data_array))
                cmap_range_max=int(np.amax(data_array))
            if var == "kt":
                var_name=r"Thermal conductivity"
                var_cmap="coolwarm"
                data_array=kt
                unit=r'$k_{\mathrm{eff}} \, [\mathrm{W} \, \mathrm{m}^{-1} \, \mathrm{K}^{-1}]$'
                cmap_range_min=int(np.amin(data_array))
                cmap_range_max=int(np.amax(data_array))
            if var == "fm":
                var_name=r"Silicate melt fraction"
                var_cmap="Greens"
                data_array=fm
                unit=r'$\varphi_{\mathrm{silicate}}$'
                cmap_range_min=0.0
                cmap_range_max=1.0
                
            # Define output figure name
            NUMBER=(FILE[:-4].split("_"))[-1]
            NUMBER="%04d" % (float(NUMBER),)
            # figure_name="2d_"+RUN+"_"+NUMBERING+"_"+var+"_"+specific_map+".pdf" # cmap tests
            figure_name=RUN+"_"+var+"_"+NUMBER+".png"
            
            plt.title(var_name+r'$, \, \mathit{t} = '+t_current+' \, \mathrm{yr}$', fontsize=28)
            # plt.xlabel(r'$\mathrm{x} \, [\mathrm{km}]$', fontsize=26)
            # plt.ylabel(r'$\mathrm{y} \, [\mathrm{km}]$', fontsize=26)
            plt.xlabel('x [km]', fontsize=26)
            plt.ylabel('y [km]', fontsize=26)

            nodex_pr = xnumx-1
            nodey_pr = ynumy-1

            plt.xticks([0, int((1./5.)*nodex_pr), int((2./5.)*nodex_pr), int((3./5.)*nodex_pr), int((4./5.)*nodex_pr), int((5./5.)*nodex_pr)], ("0", str(int((1./5.)*xsize/1e3)), str(int((2./5.)*xsize/1e3)), str(int((3./5.)*xsize/1e3)), str(int((4./5.)*xsize/1e3)), str(int((5./5.)*xsize/1e3))), fontsize = 22)
            plt.yticks([0, int((1./5.)*nodey_pr), int((2./5.)*nodey_pr), int((3./5.)*nodey_pr), int((4./5.)*nodey_pr), int((5./5.)*nodey_pr)], ("0", str(int((1./5.)*ysize/1e3)), str(int((2./5.)*ysize/1e3)), str(int((3./5.)*ysize/1e3)), str(int((4./5.)*ysize/1e3)), str(int((5./5.)*ysize/1e3))), fontsize = 22)


            if var != "comp":
                norm = mpl.colors.Normalize(vmin=cmap_range_min, vmax=cmap_range_max)
                im = ax.imshow(data_array, cmap=var_cmap, norm=norm)
                # create an axes on the right side of ax. The width of cax will be 5%
                # of ax and the padding between cax and ax will be fixed at 0.05 inch.
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="7%", pad=0.1)
                if var == "por":
                    cbar = plt.colorbar(im, cax=cax, norm=norm, label=unit, ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75])
                else:
                    cbar = plt.colorbar(im, cax=cax, norm=norm, label=unit)
                cbar.set_label(unit,size=24)
                cbar.ax.tick_params(labelsize=20) 

                if var == "acc":
                    cbar.cmap.set_under('w')
                    cbar.set_clim(0, cmap_range_max)

            print var, cmap_range_min, cmap_range_max#, data_array

            if var == "comp":
                norm = mpl.colors.Normalize(vmin=cmap_range_min, vmax=cmap_range_max)
                im = ax.imshow(data_array, cmap=var_cmap, norm=norm)
                # im = ax.imshow(data_array, cmap=var_cmap)
                box1 = DrawingArea(50, 20, 0, 0)
                rt1 = plt.Rectangle((0, 5), 20, 13, fc="#1f78b4")
                rt2 = plt.Rectangle((25, 5), 20, 13, fc="#a6cee3")
                box1.add_artist(rt1)
                box1.add_artist(rt2)
                box2 = TextArea("Solid silicates", textprops=dict(color="k", size=20))
                box = HPacker(children=[box1, box2],
                              align="center",
                              pad=0, sep=5)
                anchored_box1 = AnchoredOffsetbox(loc=2,
                                                 child=box, pad=0.,
                                                 frameon=False,
                                                 bbox_to_anchor=(0.025, 0.98),
                                                 bbox_transform=ax.transAxes,
                                                 borderpad=0.,
                                                 )
                ax.add_artist(anchored_box1)

                box1 = DrawingArea(50, 20, 0, 0)
                rt1 = plt.Rectangle((0, 5), 20, 13, fc="#e31a1c")
                rt2 = plt.Rectangle((25, 5), 20, 13, fc="#fb9a99")
                box1.add_artist(rt1)
                box1.add_artist(rt2)
                box2 = TextArea("Molten silicates", textprops=dict(color="k", size=20))
                box = HPacker(children=[box1, box2],
                              align="center",
                              pad=0, sep=5)
                anchored_box1 = AnchoredOffsetbox(loc=2,
                                                 child=box, pad=0.,
                                                 frameon=False,
                                                 bbox_to_anchor=(0.025, 0.94),
                                                 bbox_transform=ax.transAxes,
                                                 borderpad=0.,
                                                 )
                ax.add_artist(anchored_box1)

                box1 = DrawingArea(50, 20, 0, 0)
                rt1 = plt.Rectangle((0, 5), 20, 13, fc="#000000")
                rt2 = plt.Rectangle((25, 5), 20, 13, fc="#4c4c4c")
                box1.add_artist(rt1)
                box1.add_artist(rt2)
                box2 = TextArea("Iron", textprops=dict(color="k", size=20))
                box = HPacker(children=[box1, box2],
                              align="center",
                              pad=0, sep=5)
                anchored_box1 = AnchoredOffsetbox(loc=2,
                                                 child=box, pad=0.,
                                                 frameon=False,
                                                 bbox_to_anchor=(0.025, 0.90),
                                                 bbox_transform=ax.transAxes,
                                                 borderpad=0.,
                                                 )
                ax.add_artist(anchored_box1)

                # Show colorbar to check marker values
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="7%", pad=0.1)
                cbar = plt.colorbar(im, cax=cax, norm=norm, label=unit, ticks=[0, 5, 6, 7, 8, 25, 26])
                cbar.set_label(unit,size=24)
                cbar.ax.tick_params(labelsize=20)


            ttl = ax.title
            ttl.set_position([.5, 1.03])

            plt.savefig(image_dir+figure_name, bbox_inches='tight')
            # plt.savefig(image_dir2+figure_name, bbox_inches='tight')
            
            plt.close()

