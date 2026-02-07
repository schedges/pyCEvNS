import numpy as np
import sys
from utils.constants import constants

class TargetObject():
  def __init__(self,*,symbol,Z,A,nuc_mass_MeV,atom_mass_MeV):
    self.symbol = symbol
    self.Z = Z
    self.A = A
    self.N = A - Z
    self.nuc_mass_MeV = nuc_mass_MeV
    self.atom_mass_MeV = atom_mass_MeV
    self.atom_frac = None
    self.mass_frac = None
    self.num_atoms = None
    self.mass = None
    self.sn = 0
    self.sp = 0

def makeTargets(*,config,mass_dict,spinDensity_dict):
  targetObjects = {}
  detector_block = config["detector"]

  #Set names, Z, nuclear/atomic masses
  for isotope in detector_block["isotopes"]:
    if not isotope in mass_dict:
      print("Error! Symbol {isotope} not found in mass dictionary! Exiting!")
      sys.exit()
    targetObjects[isotope] = TargetObject(symbol=isotope,
                                          Z=int(mass_dict[isotope]["Z"]),
                                          A=int(mass_dict[isotope]["A"]),
                                          nuc_mass_MeV=mass_dict[isotope]["nuc_mass_MeV"],
                                          atom_mass_MeV = mass_dict[isotope]["atom_mass_MeV"])

  #Normalize & set atomic fractions
  detector_block["atom_fracs"] = np.asarray(detector_block["atom_fracs"])
  frac_sum = np.sum(detector_block["atom_fracs"])
  if frac_sum <= 0:
    print("Error! 'atom_fracs' sum must be > 0")
    sys.exit()
  if abs(frac_sum - 1.0) > 1e-6:
    print("Notice: supplied 'atom_fracs' not normalized. Normalizing to 1")
  for i,isotope in enumerate(detector_block["isotopes"]):
    targetObjects[isotope].atom_frac = detector_block["atom_fracs"][i] / frac_sum
  
  #Calculate and set atomic mass fractions
  avg_atomic_mass = 0
  for target in targetObjects.values():
    avg_atomic_mass += target.atom_frac * target.atom_mass_MeV
  for target in targetObjects.values():
    target.mass_frac = target.atom_frac*target.atom_mass_MeV / avg_atomic_mass

  #Now calculate the number of atoms
  for target in targetObjects.values():
    target.mass = detector_block["mass_kg"]*target.mass_frac
    target.num_atoms = target.mass/(target.atom_mass_MeV*constants["MeV_to_kg"])

  #Sn/Sp
  for target in targetObjects.values():
    if target.symbol in spinDensity_dict:
      target.sn = spinDensity_dict[target.symbol]["Sn"]
      target.sp = spinDensity_dict[target.symbol]["Sp"]
  return targetObjects
