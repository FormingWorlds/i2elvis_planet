# Notes For Code Understanding
> J. Eatson

## The Burning Questions

- How to modify initial parameters?
  - Currently set for marslike planet
- What is a marker in this instance?
  - Data read in from loadmart.c, loadconf()
  - Very fortlike read commands, makes sense considering provenance

## What Each Input File Does

### Mode.t3c

- Appears to contain the load file (from initialisation program in2mart.c)
- Includes location of save files and simulation parameters after(?) each time step
  - cyc0max_maxxystep_maxtkstep_maxtmstep	
  - So mainly timestep stuff
- In this case we see that dd_0.prn is the default main file, and the code reads that one in
- Difficult to change save types? Seems extremely redundant if you can loop it but does allow for fine grained control
- Includes a series of system paramters at the end, using key:
  - <value>-<key>
  - General, stokes and heat transfer parameters
- Lots of manual read ins, not deterministic at all, needs better documentation?

### Stop.yn
- Some kind of halt file? Checks every now and then and quits if Y writen to file, seems redundant, though I can imagine that its a good way to close the code without a hard kill


## in2mart.c

- 

## impact.c

- Stores impact growth data, specifically the data for hydrous_silicates.t3c
- hydrous_silicates store guide:
  -  0: T(MYr)
  -  1: Solid fraction
  -  2: Liquid fraction
  -  3: Hydrous fraction
  -  4: Primitive fraction
  -  5: N2/O2 fraction
  -  6: CO/Cl fraction
  -  7: H2O fraction
  -  8: 
  -  9: 
  - 10: 
  - ...
  - 15: maxtk        Mean temperature of entire grid(?)
  - 16: T_max_body   Max temperature of entire body
  - 17: meantk:      Mean temeprature of entire grid(?) 
  - 18: t_mean_body: Mean temperature in entire body
  - 19: Count_toohot: 

- I have modified the code to write out using scientific notation rather than floats, personal preference

## init.t3c
- Initialisation file, ideally the only thing that has to be changed signfiicantly
- I have modified this to contain an Iron core, radius 30km, 10x fe60 distribution and 0 al26, I have also removed the time offset of 1Myr, so no decay occurs

## Main Execution Loop

- Set initial timestep
  - Uses titerate() in heatmart.c
- Perform heat transport
- Move marker(? Check marker in cell method)
- Reset (VX,VY,P) Values (Why?)
- Calcualte impacts, impact()
- Handle pebble accretion, pebbleaccr()
- Calculate ro[], nu[] (Meaning?)
- Check core evolution, grain site evolution
- Save impact and core history
- Update time and exit if final time is reached

