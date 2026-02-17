import numpy as np
import math
from utils.constants import constants
from utils import formFactorUtils

class CEvNSObject():
  def __init__(self,*,Enrs_MeV,Enus_MeV,Enu_step_size_MeV,Enr_step_size_MeV,flux):
    self.Enrs_MeV = Enrs_MeV
    self.Enus_MeV = Enus_MeV
    self.Enrs_keV = Enrs_MeV*1000
    self.Enu_step_size_MeV = Enu_step_size_MeV
    self.Enr_step_size_MeV = Enr_step_size_MeV
    self.flux = flux

    self.phys_params = {}
    self.target_params = {}

    self.total_xs_cm2 = {}
    self.recoil_kernel = {}
    self.flux_avg_xs_cm2 = {}
    self.normalized_recoil_spectrum_MeV = {}
    self.normalized_recoil_spectrum_keV = {}
    self.total_counts = {}
    self.axial_ff_exp = {}
    self.vector_ff_exp = {}

#Calculate xs vs recoil energy for the given neutrino energy
def calc_recoil_xs(*,E_nu,Enrs_MeV,target):
  #Load constants
  GF_MeV2 = constants["GF_MeV2"]
  gp_V = constants["gp_V"]
  gn_V = constants["gn_V"]
  gp_A = constants["gp_A"]
  gn_A = constants["gn_A"]
  hbarc_MeVcm = constants["hbarc_MeVcm"]

  M = target.nuc_mass_MeV
  Z = target.Z
  N = target.N
  sn = target.sn
  sp = target.sp

  #Calculate q^2
  q_s_per_cm = np.sqrt(2*M*Enrs_MeV + np.power(Enrs_MeV,2))/constants["hbarc_MeVcm"]
  q_s_per_fm = q_s_per_cm * 1.e-13

  #Energy/momentum transfer specific constants
  vector_form_factors = formFactorUtils.calcFormFactors(target,q_s_per_fm,constants["vector_ff_type"])
  axial_form_factors = formFactorUtils.calcFormFactors(target,q_s_per_fm,constants["axial_ff_type"])

  #Calculate the XS at each possible KE 
  xs = np.zeros(Enrs_MeV.size)

  prefix = math.pow(GF_MeV2,2)/(2.0*math.pi)
  G_V = (gp_V*Z + gn_V*N) * vector_form_factors
  G_A = (gp_A*2*sp + gn_A*2*sn) * axial_form_factors
  Ts = Enrs_MeV
  Tmax_MeV = 2*E_nu**2/(M+2*E_nu)

  #Restrict to valid KEs for this target
  valid = (Ts <= Tmax_MeV)
  Ts = Ts[valid]
  vector_form_factors = vector_form_factors[valid]
  axial_form_factors = axial_form_factors[valid]
  G_V = G_V[valid]
  G_A = G_A[valid]

  #Calculate the xs
  valid_xs = prefix * M * (np.power(G_V+G_A,2) + np.power(G_V-G_A,2) * np.power(1-Ts/E_nu,2) - (G_V**2 - G_A**2)*M*Ts / (E_nu**2))
  valid_xs = valid_xs * hbarc_MeVcm * hbarc_MeVcm #convert to cm2
  
  #Set only the calculated bins of valid_xs to xs, the rest should be zero
  xs[:len(valid_xs)] = valid_xs

  #Calc ff expectation valid, for debugging
  ff_exp = {}
  if len(valid_xs) > 1:
    ff_exp["vector_ff"] = np.trapz(vector_form_factors * valid_xs, Ts) / np.trapz(valid_xs, Ts)
    ff_exp["axial_ff"] = np.trapz(axial_form_factors * valid_xs, Ts) / np.trapz(valid_xs, Ts)
  else:
    ff_exp["vector_ff"] = 1.
    ff_exp["axial_ff"] = 1.

  return xs,ff_exp

  
##################################
##CEvNS Cross Section Calculator##
##################################
#Based on https://arxiv.org/pdf/1803.09183, ignoring radiative corrections
def calc_xs(*,fluxObject,targetObjects,config):

  Enu_step_size_MeV = config["binning"]["Enu_step_size_MeV"]
  Enr_step_size_MeV = config["binning"]["Enr_step_size_MeV"]

  #Calculate max KE of any target to form consistent recoil energy grid
  lightest_target_nuc_mass_MeV = np.inf
  for target in targetObjects.values():  
    if target.nuc_mass_MeV < lightest_target_nuc_mass_MeV:
      lightest_target_nuc_mass_MeV = target.nuc_mass_MeV
  Enu_max = fluxObject.interpolated_energies_MeV[-1]
  T_max_MeV = 2*math.pow(Enu_max,2)/(lightest_target_nuc_mass_MeV+2*Enu_max)

  #Make KE grid
  Enrs_MeV = np.arange(Enr_step_size_MeV/2.,T_max_MeV+Enr_step_size_MeV/2.,Enr_step_size_MeV)
  
  #Make our output object
  result = CEvNSObject(Enrs_MeV=Enrs_MeV,
                       Enus_MeV=fluxObject.interpolated_energies_MeV,
                       Enr_step_size_MeV=Enr_step_size_MeV,
                       Enu_step_size_MeV=Enu_step_size_MeV,
                       flux=fluxObject.normalized_flux
                       )

  #Step through targets, calculating XS and filling recoil distribution grid (kernel)
  for targetName,target in targetObjects.items():    
    
    #Since we're storing this in the target, we should also store the energies grid with the target
    result.total_xs_cm2[targetName] = np.zeros(len(result.Enus_MeV))
    result.recoil_kernel[targetName] = np.zeros((len(result.Enus_MeV),len(result.Enrs_MeV)))
    result.target_params[targetName] = {"nuc_mass_MeV":target.nuc_mass_MeV,"sn":target.sn,"sp":target.sp}
    result.axial_ff_exp[targetName] = []
    result.vector_ff_exp[targetName] = []

    for i,E_nu in enumerate(result.Enus_MeV):

      #Calculate cross section vs. energy
      xs,ff_exp = calc_recoil_xs(E_nu = E_nu,Enrs_MeV=result.Enrs_MeV, target=target)
      result.recoil_kernel[targetName][i] = xs
      result.total_xs_cm2[targetName][i] = np.trapz(xs, result.Enrs_MeV)
      #TODO: add extra contribution from tail to total_xs
      
      result.axial_ff_exp[targetName].append(ff_exp["axial_ff"])
      result.vector_ff_exp[targetName].append(ff_exp["vector_ff"])

  flux_pdf = result.flux/np.trapz(result.flux,result.Enus_MeV)
  for targetName in targetObjects:

    #Calculate flux-averaged cross section
    result.flux_avg_xs_cm2[targetName] = np.trapz(flux_pdf * result.total_xs_cm2[targetName], result.Enus_MeV)

    #Calculate recoil spectrum
    recoil_spectrum = np.trapz(result.flux[:, None] * result.recoil_kernel[targetName], result.Enus_MeV, axis=0)
    result.normalized_recoil_spectrum_MeV[targetName] = recoil_spectrum*target.num_atoms
    result.normalized_recoil_spectrum_keV[targetName] = result.normalized_recoil_spectrum_MeV[targetName] / 1000

    #Total rate
    result.total_counts[targetName] = np.trapz(result.normalized_recoil_spectrum_MeV[targetName], result.Enrs_MeV)

  return result
