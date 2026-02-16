#Constants--can be overwritten by config file
import math

VERBOSE=1

constants = {
  "m_e_amu":0.000548579905,
  "amu_to_MeV": 931.49410242,
  "MeV_to_kg": 1.78266192*math.pow(10,-30),
  "MeV_per_joule": 6241509343260.2,
  "hbarc_MeVcm": 197.327e-13,
  "GF_MeV2": 1.1663788*math.pow(10,-11), #2024 PDG value 
  "g_vu_LL": 0.3457, #Neutrino-quark couplings, from PDG 2024 edition
  "g_vd_LL": -0.4288, #Neutrino-quark couplings, from PDG 2024 edition
  "g_vu_LR": -0.1553, #Neutrino-quark couplings, from PDG 2024 edition
  "g_vd_LR": 0.0777, #Neutrino-quark couplings, from PDG 2024 edition
  "gA_u": 0.862, #Quark spin charges, from https://arxiv.org/pdf/1909.00485 Table XI
  "gA_d": -0.424,
  "gA_s": -0.0458,
  "axial_ff_type": "uniform",
  "vector_ff_type": "uniform",
  #Helm form factor constants
  "helm_params": {
    "a":0.52, #fm
    "s":0.9 #fm
  },
  
  #Reactor specific
  #From https://arxiv.org/pdf/1212.6625
  "reactor_energy_per_fission_MeV": {
    "235U": 201.92,
    "238U": 205.52,
    "239Pu": 210.99,
    "241Pu": 213.6
  }
}

#Calculated quantities, this fn is called after parsing the input file
def updateConstants():
  global constants
  #Calculate vector couplings to protons and neutrons
  constants["gp_V"] = 2*(constants["g_vu_LL"] + constants["g_vu_LR"]) + (constants["g_vd_LL"]+constants["g_vd_LR"])
  constants["gn_V"] = (constants["g_vu_LL"] + constants["g_vu_LR"]) + 2*(constants["g_vd_LL"]+constants["g_vd_LR"])

  #Calculate axial couplings to protons and neutrons
  constants["gp_A"] = 0.5*(constants["gA_u"] - constants["gA_d"] - constants["gA_s"])
  constants["gn_A"] = 0.5*(-constants["gA_u"] + constants["gA_d"] - constants["gA_s"])

z_dict = {
  "H": 1, "He": 2,
  "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
  "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
  "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
  "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
  "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
  "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
  "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
  "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
  "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
  "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
  "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100
}
def get_ZA_from_symbol(symbol):
  i = 0
  while i < len(symbol) and symbol[i].isdigit():
    i += 1
    
  A = int(symbol[:i])
  el = symbol[i:]

  Z = z_dict[el]
  return Z, A