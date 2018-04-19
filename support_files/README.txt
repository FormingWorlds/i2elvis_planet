### Minimum working example to get the code up and running

1) Clone to Euler: 
	>> git clone git@github.com:timlichtenberg/i2elvis_planet.git

2) Compile code:
	>> sh compile.sh

3) Generate initial conditions (*_0.prn file):
	>> ./in2mart

4) Submit (chain-)job to Euler:
	>> sh submitjobs.sh


### Plot the output

1) Copy .prn files to local directory

2) Define correct 'PDIR' and 'image_dir' in plot2d.py

3) Run plotting script:
	>> python plot2d.py
