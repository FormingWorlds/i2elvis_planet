#!/usr/bin/env python3

"""
Set up simulation in separate folder with requisite files
Creates a new init.t3c file, traverses to directory and executes it
"""

t3c_files = ["amir.t3c",
             "core_dynamo.t3c",
             "core_stuff.t3c",
             "file.t3c",
             "hydrous_silicates.t3c",
             "impact_history.t3c",
             "impact_mass.t3c",
             "impact_no.t3c",
             "mode.t3c",
             "pebble_accr.t3c",
             "pebble_history.t3c",
             "start_mass.t3c",
             "mars_crust",
             "mars_mantle"]

from string import Template
from shutil import copyfile
from os import mkdir
from subprocess import Popen
import sys
import argparse
from tqdm import tqdm

def gen_sim(mant_ext,core_ext,al,fe,exit_time,work_dir):
  # Quick checks for input conditions
  if mant_ext <= 0.:
    print("!!! Zero or negative mantle size! Check inputs !!!")
    sys.exit()
  if al < 0 or fe < 0:
    print("!!! Cannot have negative abundance ratios !!!")
    sys.exit()
  if mant_ext <= core_ext:
    print("!!! Core/Mantle ratio > 1, nonsensical !!!")
    sys.exit()
  try:
    float(mant_ext)
    float(core_ext)
    float(al)
    float(fe)
    str(work_dir)
  except TypeError:
    print("!!! Input parameters generally wrong, see help with -h command !!!")
    sys.exit()
  # Begin!
  print("! Initialising environment at {}".format(work_dir))
  repo_path = sys.path[0]
  # Copy files
  print("! Copying files")
  try:
    mkdir(work_dir)
    print("! Created {}!".format(work_dir))
  except FileExistsError:
    print("! Already created {}!".format(work_dir))

  try:
    for file in tqdm(t3c_files):
      src = "{}/{}".format(repo_path,file)
      dst = "{}/{}".format(work_dir,file)
      copyfile(src,dst)
  except FileNotFoundError:
    raise("Cannot find files, make sure they are in the same directory as this script!")
  
  print("! Creating init file... ",end="")
  # Read in template file
  template_filename = "{}/init-template.pytemp".format(repo_path)
  try:
    with open(template_filename) as template_file:
      template_text = template_file.read()
      template = Template(template_text)
  except FileNotFoundError:
    raise("Template file not found! Make sure it is in the same directory as this script")

  # Begin modifying template
  init_params = {}
  # Format string for exit time
  init_params["exit_time"] = "{:.3e}-timeexit(yr)".format(exit_time)
  # Calculate abundance ratio from Solar System values
  al26_ss_abun = 5.250 # Al_26/Al_27 solar system ratio, in units of 1e-5
  fe60_ss_abun = 1.150 # Fe_60/Fe_56 solar system ratio, in units of 1e-8
  al26_abun = al * al26_ss_abun
  fe60_abun = fe * fe60_ss_abun
  # Format strings
  init_params["al_abun"] = "{:.3f}-al2627_init(e-5,ss_0=5.250)".format(al26_abun)
  init_params["fe_abun"] = "{:.3f}-fe6056_init(e-8,ss_0=1.150)".format(fe60_abun)
  # Calculate core and mantle size
  comp = ""
  core_ext_m = int(core_ext * 1000)
  mant_ext_m = int(mant_ext * 1000)
  init_params["core_extent"] = "{}".format(core_ext_m)
  init_params["mantle_extent"] = "{}".format(mant_ext_m)
  # Calculate simulation size, which is twice the diameter of the planetisemal
  sim_ext_m = int(mant_ext_m * 4)
  # Check to see if a core is included
  if core_ext > 0:
    comp += "/Iron_1\n"
    comp += " 3 7 0.5 0.5 0.5 0.5 m{} m0 0 360\n".format(core_ext_m)
  # Build mantle
  comp += "/Wet_Silicates\n"
  comp += " 3 6 0.5 0.5 0.5 0.5 m{} m{} 0 360\n".format(mant_ext_m,core_ext_m)
  # Build surrounding vacuum, "air"
  comp += "/AIR\n"
  comp += " 3 0 0.5 0.5 0.5 0.5 m{} m{} 0 360".format(sim_ext_m,mant_ext_m)
  init_params["composition"] = comp
  # Format string for GRID_PARAMETERS_DESCRIPTIONS
  init_params["x_size"] = "{}-xsize(m)".format(sim_ext_m)
  init_params["y_size"] = "{}-ysize(m)".format(sim_ext_m)
  # Convert template and write
  init_data = template.substitute(init_params)
  init_filename = "{}/init.t3c".format(work_dir)
  with open(init_filename,"w") as init_file:
    init_file.write(init_data)
  print("! Done!\nRunning initialisation programme")

  # Execute in2mart
  in2elvis = "{}/in2mart".format(repo_path)
  initialise = Popen([in2elvis],cwd=work_dir)
  # initialise.wait()

  print("!!! Environment setup in {} finished !!!".format(work_dir))
  return

if __name__ == "__main__":
  parse = argparse.ArgumentParser(description="Initialise the i2elvis code in a separate folder")
  parse.add_argument("planet_size",type=float,help="Planetesimal Size (km)")
  parse.add_argument("core_size",type=float,help="Core Size (km)")
  parse.add_argument("al26_ratio",type=float,help="Al26/Al27 ratio comparative to solar system ratio")
  parse.add_argument("fe60_ratio",type=float,help="Fe60/Fe56 ratio comparative to solar system ratio")
  parse.add_argument("-e","--exit_time",type=float,default=15e6,help="Exit time in years, default is 15Myr")
  parse.add_argument("-f","--folder_name",type=str,help="Folder name, defaults to sim_p_<planet_size>_c_<core_size>_al_<al26_ratio>_fe_<fe60_ratio>")
  args = parse.parse_args()
  # Use argparse as a wrapper
  mant_ext  = args.planet_size
  core_ext  = args.core_size
  al        = args.al26_ratio
  fe        = args.fe60_ratio
  exit_time = args.exit_time
  # Make folder name if none specified
  if args.folder_name == None:
    work_dir = "sim_p_{}_c_{}_al_{}_fe_{}".format(int(mant_ext),int(core_ext),al,fe)
  else:
    work_dir = args.folder_name
  # Start
  gen_sim(mant_ext,core_ext,al,fe,exit_time,work_dir)

