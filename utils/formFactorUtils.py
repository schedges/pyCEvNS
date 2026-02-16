import numpy as np
from utils.constants import constants
from scipy.special import spherical_jn

def calcHelmFormFactors(target,qs_per_fm):
  a = constants["helm_params"]["a"]
  s = constants["helm_params"]["s"]
  A = target.A

  c = 1.23*np.power(A,1./3.) - 0.6
  R_sq = c*c + 7./3. * np.pi*np.pi * a*a - 5*s*s
  R = np.sqrt(R_sq)
  qRs = qs_per_fm*R

  t1s = 3*spherical_jn(1,qRs)/qRs 
  t1s = np.where(abs(qRs)<1e-10,1.0,t1s) #Catch small q
  t2s = np.exp(-np.power(qs_per_fm*s,2)/2.)
  ffs = t1s*t2s

  return ffs


def calcFormFactors(target,qs_per_fm,ffType):
  if ffType=="uniform":
    return np.ones_like(qs_per_fm)
  if ffType=="helm":
    return calcHelmFormFactors(target,qs_per_fm)