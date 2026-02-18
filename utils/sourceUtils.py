#REQUIRES UNFORM ENERGY BINNING FOR RAW CSV, GIVEN AS FLUX/MeV!!!
import numpy as np
import math
import sys
from utils import ioUtils
from utils.constants import constants, VERBOSE

class FluxObject():
  def __init__(self,raw_energies,raw_fluxes):
    #These are required of all flux objects
    self.raw_energies_MeV = raw_energies
    self.raw_flux = raw_fluxes
    self.interpolated_energies_MeV = None
    self.interpolated_flux = None
    self.normalized_flux = None

def makeReactorFluxObjects(config):
  fluxObjects = {}
  source_params_block = config["source"]["params"]

  #Load spectra for individual core components
  for isotope, iso_block in source_params_block["core_isotopes"].items():
    raw_energies_MeV, raw_flux = ioUtils.loadFluxCSV(iso_block["filename"])
    fluxObjects[isotope] = FluxObject(raw_energies_MeV, raw_flux)
    #TODO: Check if raw_energies_MeV sorted, if not return error and exist

  #Get global max value for defining energy grid
  maxE = 0
  for fluxObject in fluxObjects.values():
    if max(fluxObject.raw_energies_MeV) > maxE:
      maxE = max(fluxObject.raw_energies_MeV)

  #Define energy grid
  interpolated_bin_size = config["binning"]["Enu_step_size_MeV"]
  interpolated_energies_MeV = np.arange(0+interpolated_bin_size/2.,maxE+interpolated_bin_size/2.,interpolated_bin_size)
  
  #Interpolate fluxes
  for fluxObject in fluxObjects.values():
    fluxObject.interpolated_energies_MeV = interpolated_energies_MeV
    fluxObject.interpolated_flux = np.interp(interpolated_energies_MeV,fluxObject.raw_energies_MeV,fluxObject.raw_flux,left=0.,right=0.)
    fluxObject.interpolated_bin_size = interpolated_energies_MeV[1] - interpolated_energies_MeV[0]

  #Debugging
  if VERBOSE>0:
    for isotope, fluxObject in fluxObjects.items():
      print(f"Per {isotope} fission, raw values give {np.trapz(fluxObject.raw_flux,x=fluxObject.raw_energies_MeV):.2f} neutrinos/fission")
      print(f"Per {isotope} fission, interpolated values give {np.trapz(fluxObject.interpolated_flux,x=fluxObject.interpolated_energies_MeV):.2f} neutrinos/fission")

  #Normalize to neutrinos/cm2/sec
  power_GW = source_params_block["power_GWth"]
  joules_per_sec = power_GW*math.pow(10,9)
  MeV_per_sec = joules_per_sec*constants["MeV_per_joule"]
  
  avg_energy_per_fission_MeV = 0
  for isotope, iso_block in source_params_block["core_isotopes"].items():
    avg_energy_per_fission_MeV += (constants["reactor_energy_per_fission_MeV"][isotope] * iso_block["fission_frac"])
  fissions_per_second = MeV_per_sec / avg_energy_per_fission_MeV

  distance_factor = 1./(4*math.pi*math.pow(source_params_block["distance_m"]*100,2))
  for isotope, fluxObject in fluxObjects.items():
    fission_frac = source_params_block["core_isotopes"][isotope]["fission_frac"]
    fluxObject.normalized_flux = (fluxObject.interpolated_flux * fissions_per_second * fission_frac * distance_factor)
  
  #Calculate total flux
  totalFlux = FluxObject(np.array([]),np.array([]))
  for i,fluxObject in enumerate(fluxObjects.values()):
    if i==0:
      totalFlux.raw_energies_MeV = fluxObject.raw_energies_MeV.copy()
      totalFlux.interpolated_energies_MeV = fluxObject.interpolated_energies_MeV.copy()
      totalFlux.raw_flux = fluxObject.raw_flux.copy()
      totalFlux.interpolated_flux = fluxObject.interpolated_flux.copy()
      totalFlux.normalized_flux = fluxObject.normalized_flux.copy()
    else:
      totalFlux.raw_flux += fluxObject.raw_flux
      totalFlux.interpolated_flux += fluxObject.interpolated_flux
      totalFlux.normalized_flux += fluxObject.normalized_flux
  fluxObjects["total"] = totalFlux

  #Debugging
  if VERBOSE>0:
    print(f"\n{power_GW} GWth reactor produces {fissions_per_second:.3e} fissions/sec")
    print(f"Total neutrino flux {source_params_block['distance_m']}m from a {power_GW} GWth reactor is {np.trapz(fluxObjects['total'].normalized_flux,x = fluxObjects['total'].interpolated_energies_MeV):.3e} nu/cm2/sec")
      
  return fluxObjects

def makePiDARFluxObjects(config):
  fluxObjects = {}
  Enu_vu = (math.pow(constants["m_pi+_MeV"],2) - math.pow(constants["m_mu_MeV"],2)) / (2*constants["m_pi+_MeV"])
  Enu_max = (math.pow(constants["m_mu_MeV"],2) - math.pow(constants["m_e_MeV"],2)) / (2*constants["m_mu_MeV"])

  #Define energy grid
  interpolated_bin_size = config["binning"]["Enu_step_size_MeV"]
  interpolated_energies_MeV = np.arange(0+interpolated_bin_size/2.,Enu_max+interpolated_bin_size/2.,interpolated_bin_size) #Bin centers

  #Calculate spectral shapes
  vu_obj = FluxObject(None,None)
  vu_obj.interpolated_energies_MeV = interpolated_energies_MeV
  vu_obj.interpolated_flux = np.zeros_like(interpolated_energies_MeV)

  #Find bin in interpolated energies with Enu_vu
  vu_bin = np.argmin(np.abs(interpolated_energies_MeV - Enu_vu))
  vu_obj.interpolated_flux[vu_bin] = 1./interpolated_bin_size

  vubar_obj = FluxObject(None,None)
  vubar_obj.interpolated_energies_MeV = interpolated_energies_MeV
  vubar_obj.interpolated_flux = 16*np.power(interpolated_energies_MeV,2)*(3.*constants["m_mu_MeV"]-4*interpolated_energies_MeV)/np.power(constants["m_mu_MeV"],4) 
  #don't allow flux to be negative
  vubar_obj.interpolated_flux = np.maximum(vubar_obj.interpolated_flux, 0.0)

  ve_obj = FluxObject(None,None)
  ve_obj.interpolated_energies_MeV = interpolated_energies_MeV
  ve_obj.interpolated_flux = 96*np.power(interpolated_energies_MeV,2)*(constants["m_mu_MeV"]-2.0*interpolated_energies_MeV)/np.power(constants["m_mu_MeV"],4)
  ve_obj.interpolated_flux = np.maximum(ve_obj.interpolated_flux, 0.0)
  #mask negative

  if config["source"]["flavor"]=="vu":
    fluxObjects = {"vu":vu_obj}
  elif config["source"]["flavor"]=="vubar":
    fluxObjects = {"vubar":vubar_obj}
  elif config["source"]["flavor"]=="ve":
    fluxObjects = {"ve":ve_obj}
  else:
    fluxObjects = {"vu":vu_obj,
                   "ve":ve_obj,
                   "vubar":vubar_obj}

  #Normalize to num neutrinos
  source_params_block = config["source"]["params"]
  beamPower_W = source_params_block["beam_power_MW"]*1e6
  beamEnergy_J = source_params_block["beam_energy_MeV"]/constants["MeV_per_joule"]
  protons_per_sec = beamPower_W/beamEnergy_J
  nus_per_flavor_per_sec = protons_per_sec*source_params_block["nu_per_flavor_per_proton"]
  distance_cm = source_params_block["distance_m"]*100
  nus_per_flavor_per_sec /= (4*math.pi*math.pow(distance_cm,2)) 
  for flavor_name,fluxObject in fluxObjects.items():
    #Normalize area to 1
    area = np.trapz(fluxObject.interpolated_flux, x=fluxObject.interpolated_energies_MeV)
    if area <=0:
      print(f"Error in normalization of pidar source, area of {flavor_name} spectrum <= 0, exiting!")
      sys.exit()
    fluxObject.interpolated_flux /= area
    fluxObject.normalized_flux = fluxObject.interpolated_flux * nus_per_flavor_per_sec

  #Create "all" object if all requested
  if config["source"]["flavor"]=="all":
    vall_obj = FluxObject(None,None)
    vall_obj.interpolated_energies_MeV = interpolated_energies_MeV
    vall_obj.interpolated_flux = ve_obj.interpolated_flux + vu_obj.interpolated_flux + vubar_obj.interpolated_flux
    vall_obj.normalized_flux = ve_obj.normalized_flux + vu_obj.normalized_flux + vubar_obj.normalized_flux
    fluxObjects["total"] = vall_obj

  return fluxObjects

def makeFluxObjects(*,config):
  if config["source"]["name"]=="reactor":
    fluxObjects = makeReactorFluxObjects(config)
  elif config["source"]["name"]=="pidar":
    fluxObjects = makePiDARFluxObjects(config)
  return fluxObjects