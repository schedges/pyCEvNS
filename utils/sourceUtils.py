#REQUIRES UNFORM ENERGY BINNING FOR RAW CSV, GIVEN AS FLUX/MeV!!!
import numpy as np
import math
from utils import ioUtils
from utils.constants import constants


class FluxObject():
  def __init__(self,raw_energies,raw_fluxes):
    #These are required of all flux objects
    self.raw_energies = raw_energies
    self.raw_flux = raw_fluxes
    self.interpolated_energies = None
    self.interpolated_flux = None
    self.normalized_flux = None

def makeReactorFluxObjects(config):
  fluxObjects = {}
  source_params_block = config["source"]["params"]

  #Load spectra for individual core components
  for isotope, iso_block in source_params_block["core_isotopes"].items():
    raw_energies, raw_flux = ioUtils.loadFluxCSV(iso_block["filename"])
    fluxObjects[isotope] = FluxObject(raw_energies, raw_flux)

  #Get global max value for defining energy grid
  maxE = 0
  for fluxObject in fluxObjects.values():
    if max(fluxObject.raw_energies) > maxE:
      maxE = max(fluxObject.raw_energies)

  #Define energy grid
  interpolated_bin_size = config["binning"]["Enu_step_size_MeV"]
  interpolated_energies = np.arange(0,maxE+interpolated_bin_size,interpolated_bin_size)
  
  #Interpolate fluxes
  for fluxObject in fluxObjects.values():
    fluxObject.interpolated_energies = interpolated_energies
    fluxObject.interpolated_flux = np.interp(interpolated_energies,fluxObject.raw_energies,fluxObject.raw_flux,left=0.,right=0.)
    fluxObject.interpolated_bin_size = interpolated_energies[1] - interpolated_energies[0]

  #Debugging
  for isotope, fluxObject in fluxObjects.items():
    print(f"Per {isotope} fission, raw values give {np.trapz(fluxObject.raw_flux,x=fluxObject.raw_energies):.2f} neutrinos/fission")
    print(f"Per {isotope} fission, interpolated values give {np.trapz(fluxObject.interpolated_flux,x=fluxObject.interpolated_energies):.2f} neutrinos/fission")

  #Normalize to neutrinos/cm2/sec
  power_GW = source_params_block["power_GWth"]
  joules_per_sec = power_GW*math.pow(10,9)
  MeV_per_sec = joules_per_sec*constants["MeV_per_joule"]
  
  avg_energy_per_fission_MeV = 0
  for isotope, iso_block in source_params_block["core_isotopes"].items():
    avg_energy_per_fission_MeV += (constants["reactor_energy_per_fission_MeV"][isotope] * iso_block["fission_frac"])
  fissions_per_second = MeV_per_sec / avg_energy_per_fission_MeV

  distance_factor = 1./(4*math.pi*math.pow(source_params_block["distance_m"],2))
  for isotope, fluxObject in fluxObjects.items():
    fission_frac = source_params_block["core_isotopes"][isotope]["fission_frac"]
    fluxObject.normalized_flux = (fluxObject.interpolated_flux * fissions_per_second * fission_frac * distance_factor)
    
  totalFlux = FluxObject(np.array([]),np.array([]))
  for i,fluxObject in enumerate(fluxObjects.values()):
    if i==0:
      totalFlux.raw_energies = fluxObject.raw_energies.copy()
      totalFlux.interpolated_energies = fluxObject.interpolated_energies.copy()
      totalFlux.raw_flux = fluxObject.raw_flux.copy()
      totalFlux.interpolated_flux = fluxObject.interpolated_flux.copy()
      totalFlux.normalized_flux = fluxObject.normalized_flux.copy()
    else:
      totalFlux.raw_flux += fluxObject.raw_flux
      totalFlux.interpolated_flux += fluxObject.interpolated_flux
      totalFlux.normalized_flux += fluxObject.normalized_flux
  fluxObjects["total"] = totalFlux

  #Debugging
  print(f"\n{power_GW} GWth reactor produces {fissions_per_second:.3e} fissions/sec")
  print(f"Total neutrino flux {source_params_block['distance_m']}m from a {power_GW} GWth reactor is {np.trapz(fluxObjects['total'].normalized_flux,x = fluxObjects['total'].interpolated_energies):.3e} nu/cm2/sec")
      
  return fluxObjects

def makeFluxObjects(*,config):
  if config["source"]["name"]=="reactor":
    fluxObjects = makeReactorFluxObjects(config)
  return fluxObjects