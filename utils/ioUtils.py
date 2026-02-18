
import numpy as np
import math
import json
import sys
import os
import pandas as pd
from utils.constants import constants, updateConstants

################
##PARSE CONFIG##
################
def parseConfig(configFileName):
  with open(configFileName, "r") as f:
    lines = f.readlines()
  lines = [line.split("//", 1)[0] for line in lines]  #drop comments
  text = "".join(lines)
  config = json.loads(text)
  
  #Check all top level keys are there
  required_top_keys = ["detector","source","binning","output"]
  for key in required_top_keys:
    if key not in config:
      print(f"Error! Config missing required key: '{key}'")
      sys.exit()

  #############################
  #Check detector params valid#
  #############################
  #all keys there
  detector_block = config["detector"]
  required_keys = ["isotopes", "atom_fracs", "mass_kg"]
  for k in required_keys:
    if k not in detector_block:
      print(f"Error! 'detector' section missing required key: '{k}'")
      sys.exit()
  #number of isotopes same as num fractions
  if len(detector_block["isotopes"]) != len(detector_block["atom_fracs"]):
    print("Error! 'isotopes' and 'atom_fracs' must be the same length")
    sys.exit()

  ###########################
  #Check source params valid#
  ###########################
  #all keys there
  source_block = config["source"]
  required_keys = ["name", "flavor"]
  for k in required_keys:
    if k not in source_block:
      print(f"Error! 'source' section missing required key: '{k}'")
      sys.exit()
  #Valid sources
  valid_sources = ["reactor","pidar"] #More to be implemented
  if not source_block["name"] in valid_sources:
    print(f"Error! Invalid source specified in config file of {source_block['name']}!\nValid sources are:")
    for source in valid_sources:
      print(f"\t{source}")
    sys.exit()
  source_params_block = source_block["params"]

  ##Reactor specific source params##
  if source_block["name"] == 'reactor':
    #All keys there
    required_keys = ["distance_m","power_GWth","core_isotopes"]
    for k in required_keys:
      if k not in source_params_block:
        print(f"Error! reactor source section missing required key: '{k}'")
        sys.exit()
    #Check flavor is valid
    valid_flavors = ["vebar"]
    if not source_block["flavor"] in valid_flavors:
      print(f"Error! for 'reactor' source, valid flavors are:")
      for flavor in valid_flavors:
        print(flavor)
      print(f"You put {source_block['flavor']}")
      sys.exit()
    #Check isotopes are valid--we need energy release per fission defined for proper normalization
    allowed_isotopes = ["235U","238U","239Pu","241Pu"] 
    for isotope, iso_block in source_params_block["core_isotopes"].items():
      if not isotope in allowed_isotopes:
        print(f"Error! Isotope {isotope} not currently supported!")
        print("If this is not a mistake, add to 'allowed_isotopes' in ioUtils and specify in constants.py reactor_energy_per_fission_MeV dict")
        sys.exit()
    #Check each of the isotopes has fission_frac and filename
    required_isotope_keys = ["fission_frac","filename"]
    for isotope, iso_block in source_params_block["core_isotopes"].items():
      for k in required_isotope_keys:
        if not k in iso_block:
          print(f"Error! Each core isotope must have 'fission_frac' and 'filename' defined")
          sys.exit()
    #check isotope filenames exist
    for isotope, iso_block in source_params_block["core_isotopes"].items():
      fname = iso_block["filename"]
      if not os.path.exists(fname):
        print(f"Error! Missing reactor spectrum file: {fname}")
        sys.exit()

  #pidar source
  elif source_block["name"]=="pidar":
    #Check flavors
    valid_flavors = ["vubar","vu","ve","all"]
    if not source_block["flavor"] in valid_flavors:
      print(f"Error! for 'pidar' source, valid flavors are:")
      for flavor in valid_flavors:
        print(flavor)
      print(f"You put {source_block['flavor']}")
      sys.exit()
    #Check keys
    required_keys = ["beam_energy_MeV","beam_power_MW","nu_per_flavor_per_proton","distance_m"]
    for k in required_keys:
      if k not in source_params_block:
        print(f"Error! pidar source section missing required key: '{k}'")
        sys.exit()

  ############################
  #Check binning params valid#
  ############################
  required_keys = ["Enu_step_size_MeV","Enr_step_size_MeV"]
  for k in required_keys:
    if k not in config["binning"]:
      print(f"Error! Binning section missing required key: '{k}'")
      sys.exit()

  ###########################
  #Physics paramaeters check#
  ###########################
  #If present, check physics params are valid
  if "physics_params" in config:
    phys_block = config["physics_params"]
    for par in phys_block:
      if not par in constants:
        print(f"Error! User tried to specify physics parameter {par} but that does not exist in constants.py")
        sys.exit()
      else:
        constants[par] = phys_block[par]

    #Specific form-factor checks
    allowed_ffs = ["helm","uniform"]
    if "vector_ff_type" in phys_block:
      if not phys_block["vector_ff_type"] in allowed_ffs:
        print("Error! Allowed ff types are:")
        for ff in allowed_ffs:
          print(ff)
        print(f"Your choice of {phys_block['vector_ff_type']} vector_ff_type is invalid! Exiting!")
        print("Exiting!")
        sys.exit()

    if "axial_ff_type" in phys_block:
      if not phys_block["axial_ff_type"] in allowed_ffs:
        print("Error! Allowed ff types are:")
        for ff in allowed_ffs:
          print(ff)
        print(f"Your choice of {phys_block['axial_ff_type']} axial_ff_type is invalid! Exiting!")
        sys.exit()
  
  updateConstants()

  #########################
  #Output parameter checks#
  #########################
  #TODO

  return config


#################################
##Load nuclear data for targets##
#################################
#Loads the masses for calculating thresholds. 
def loadMassDict(massFile="data/mass_1.mas20.txt"):
  widths = [1,3,5,5,5,1,3,4,1,14,12,13,1, #a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,
            10,1,2,13,11,1,3,1,13,12]     #f10.5,1x,a2,f13.5,f11.5,1x,i3,1x,f13.6,f12.6
  names = ["skip1","N-Z","N","Z","A","skip2","El","O","skip3","mass_excess_keV","mass_excess_unc","binding_energy_keV","skip4",
           "binding_energy_unc","skip5","skip6","beta_decay_energy_keV","beta_decay_energy_unc","skip7","atom_mass_1_uamu","skip8","atom_mass_2_uamu","mass_unc"]
  df = pd.read_fwf(massFile, widths=widths, names=names, skiprows=37,dtype=str)

  df = df[["A","Z","El","atom_mass_1_uamu","atom_mass_2_uamu"]].copy()
  df["A"] = df["A"].str.strip()
  df["El"] = df["El"].str.strip()
  df["symbol"] = df["A"]+df["El"]
  df["Z"] = df["Z"].astype(int)
  df["atom_mass_1_uamu"] = df["atom_mass_1_uamu"].str.replace("#", "")
  df["atom_mass_2_uamu"] = df["atom_mass_2_uamu"].str.replace("#", "")
  df["nuc_mass_MeV"] = ((df["atom_mass_1_uamu"].astype(float)*1e6 + df["atom_mass_2_uamu"].astype(float))*1e-6 - df["Z"]*constants["m_e_amu"])*constants["amu_to_MeV"]
  df["atom_mass_MeV"] = ((df["atom_mass_1_uamu"].astype(float)*1e6 + df["atom_mass_2_uamu"].astype(float))*1e-6)*constants["amu_to_MeV"]
  df = df[["symbol","A","Z","nuc_mass_MeV","atom_mass_MeV"]].dropna()

  mass_dict = {}
  for _, row in df.iterrows():
    mass_dict[row["symbol"]] = {"Z": row["Z"], "A":row["A"],"nuc_mass_MeV": row["nuc_mass_MeV"], "atom_mass_MeV": row["atom_mass_MeV"]}

  return mass_dict

def loadSnSpDict(SnSnTableFilename="data/sn_sp.csv"):
  df = pd.read_csv(SnSnTableFilename,delimiter=",",comment="#")
  df = df.set_index("symbol")
  return df.to_dict(orient="index")

#Loads CSV of energy, neutrinos/fission/energy bin, returns np arrays
def loadFluxCSV(fname):
  df = pd.read_csv(fname,comment="#",delimiter=",")
  energy = np.array(df["energy"])
  flux = np.array(df["flux"])
  return energy,flux

