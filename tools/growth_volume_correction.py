import numpy as np
from scipy.optimize import curve_fit


def allometry(x,xdata):
  return x[0]*xdata**x[2]+x[3]

def growth_volume_correction(cylinder,inputs):
  '''
  ---------------------------------------------------------------------
  GROWTH_VOLUME_CORRECTION.M       Use growth volume allometry approach to 
                                    modify the radius of cylinders.
 
  Version 2.0.0
  Latest update     16 Sep 2021
 
  Copyright (C) 2013-2021 Pasi Raumonen
  ---------------------------------------------------------------------
 
  Use growth volume (= the total volume "supported by the cylinder") 
  allometry approach to modify the radius of too large and too small 
  cylinders. Uses the allometry: 
 
        Radius = a * GrowthVolume^b + c
 
  If cylinder's radius is over fac-times or under 1/fac-times the radius
  predicted from the growth volume allometry, then correct the radius to 
  match the allometry. However, the radius of the cylinders in the branch
  tips are never incresed, only decreased by the correction. More details 
  can be from Jan Hackenberg's "SimpleTree" papers and documents.
  ---------------------------------------------------------------------
  Inputs:
  cylinder    Structure array that needs to contains the following fields: 
    radius (Rad)        Radii of the cylinders, vector
    length (Len)        Lengths of the cylinders, vector
    parent (CPar)       Parents of the cylinders, vector
  inputs.GrowthVolFac   The factor "fac", defines the upper and lower
                          allowed radius from the predicted one:
                          1/fac*predicted_rad <= rad <= fac*predicted_rad
  ---------------------------------------------------------------------
 
  Changes from version 1.0.0 to 2.0.0, 16 Sep 2021:
  1) Changed the roles of RADIUS and GROWTH_VOLUME in the allometry, i.e.
     the radius is now predicted from the growth volume
  2) Do not increase the radius of the branch tip cylinders 
  '''
  print("--------------")
  print("Growth volume base correction of cylinder radii:")

  Rad = cylinder['radius'].astype(float)
  Rad0 = Rad
  Len = cylinder['length'].astype(float)
  CPar = cylinder['parent']
  CExt = cylinder['extension']

  initial_volume = np.round(1000*np.pi*np.sum(Rad**2*Len))
  print(f" Initial_volume (L): {initial_volume}")

  ## Define the child cylinders for each cylinder
  n = len(Rad)
  CChi = {}
  ind = np.transpose(np.array(range(n), dtype=int))

  for i in range(n):
    CChi[i] = ind[CPar==i]
  
  ## Compute the growth volume
  GrowthVol = np.zeros(n) # growth volume
  S = np.array([ len(CChi[k]) for k in CChi.keys()])
  print(f"S: {S}")
  modify = (S == 0)
  GrowthVol[modify] = np.pi*Rad[modify]**2*Len[modify]
  print(f"CChi.keys(): {CChi.keys()}")
  print(f"modify: {modify}")
  parents = np.unique(CPar[modify]).astype(int)
  print(f"1: parents: {parents}")
  if parents[0] ==0:
    parents = parents[1:]
  
  print(f"2: parents: {parents}")
  print(f"2: Rad: {Rad}")
  print(f"2: Len: {Len}")
  print(f"2: type(parents): {type(parents)}")
  print(f"2: type(Rad): {type(Rad)}")
  print(f"2: type(Len): {type(Len)}")
  print(f"2: CPar: {CPar}")
  exit()

  while len(parents)>0:
    V = np.pi*Rad[parents]**2*Len[parents]
    m = len(parents)
    for i in range(m):
      GrowthVol[parents[i]] = V[i]+np.sum(GrowthVol[CChi[parents[i]]])
    parents = np.unique(CPar[parents]).astype(int)
    print(f"2: parents: {parents}")
    if parents[0]==0:
      parents = parents[1:]
  
  ## Fit the allometry: Rad = a*GV^b
  p0 = (0.5, 0.5, 0.0)      # initial guess
  popt, pcov = curve_fit(allometry, GrowthVol, Rad, p0=p0)  # no bounds
  print(f"popt: {popt}")
  print(f"pcov: {pcov}")
  exit()
