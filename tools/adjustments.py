import numpy as np
from tools.distances_to_line import distances_to_line
from tools.surface_coverage2 import surface_coverage2
from tools.orthonormal_vectors import orthonormal_vectors
from least_squares_fitting.least_squares_circle_centre import least_squares_circle_centre

def adjustments(cyl, parcyl, inputs, Regs):
  cyl['radius'] = np.atleast_1d(cyl['radius'])
  cyl['SurfCov'] = np.atleast_1d(cyl['SurfCov'])
  cyl['length'] = np.atleast_1d(cyl['length'])
  cyl['mad'] = np.atleast_1d(cyl['mad'])
  SC = np.atleast_1d(cyl['SurfCov'])
  numC = np.size(cyl['radius'])
  Mod = np.zeros(numC,dtype=bool)

  ## Determine the maximum and the minimum radius
  # The maximum based on parent branch

  if np.size(parcyl['radius']>0):
    MaxR = 0.95*parcyl['radius']
    MaxR =np.maximum(MaxR, inputs['MinCylRad'])
  else:
    # use the maximum from the bottom cylinders
    a = np.min((3,numC))
    MaxR = 1.25*np.max(cyl['radius'][:a])
  if np.sum(SC>0.7) > 0:
    MinR = np.min((cyl['radius'][SC>0.7]))
  else:
    MinR = []
  if np.size(MinR)>0 and np.min(cyl['radius']) < MinR/2:
    MinR = np.min(cyl['radius'][SC>0.4])
  elif np.size(MinR) == 0:
    if np.sum(SC>0.4) > 0:
      MinR = np.min(cyl['radius'][SC>0.4])
    else:
      MinR = []
    if np.size(MinR) == 0:
      MinR = inputs['MinCylRad']
  

  ## Check maximum and minimum radii
  I = cyl['radius'] < MinR
  cyl['radius'][I] = MinR
  Mod[I] = True
  if inputs['ParentCor'] or numC <=3:
    I = (np.atleast_1d(cyl['radius'] > MaxR) & (SC < 0.7)) | np.atleast_1d(cyl['radius']> 1.2*MaxR)

    cyl['radius'][I] = MaxR
    Mod[I] = True

    # For short branches modify with more restrictions
    if numC <=3:
      I = ((cyl['radius'] > 0.75 * MaxR) & (SC <0.7))
      if np.any(I):
        r = np.maximum(SC[I]/0.7*np.atleast_1d(cyl['radius'])[I],MinR)
        cyl['radius'][I] = r
        Mod[I] = True
  
  ## Use taper correction to modify radius of too small and large cylinders
  # Adjust radii if a small SurfCov and high SurfCov in the previus and 
  # following cylinders

  for i in range(1,numC-1):
    if (SC[i] < 0.7) and (SC[i-1] >= 0.7) and (SC[i+1] >= 0.7):
      cyl['radius'][i] = 0.5 *(cyl['radius'][i-1] *cyl['radius'][i+1])
      Mod[i] = True
  
  ## Use taper correction to modify radius of too small and large cylinders
  if inputs['TaperCor']:
    if np.max(cyl['radius']) < 0.001:
      ## Adjust radii of thin branches to be linearly decreasing
      if numC>2:
        r = np.sort(cyl['radius'])
        r =r[1:-1]
        a = 2*np.mean(r)
        if a > np.max(r):
          a = np.min((0.01,np.max(r)))
        b=np.min((0.5*np.min(cyl['radius']),0.001))
        cyl['radius'] = np.transposed(np.linspace(a, b, numC))
      elif numC >1:
        r = np.max(cyl['radius'])
        cyl['radius'] = np.concatenate((r, 0.5*r))
      Mod = np.ones(numC, dtype=bool)
    elif numC>4:
      ## Parabola adjustment of maximum and minimum
      # Define parabola taper shape as maximum (and minimum) radii for
      # the cylinders with low surface coverage

      branchlen = np.sum(cyl['length'][:numC]) # branch length
      L = cyl['length']/2+ np.concatenate(([0], np.cumsum(cyl['length'][:numC-1])))
      Taper = np.concatenate((L, [branchlen]))
      Taper = np.column_stack((Taper,np.concatenate((1.05*cyl['radius'], [MinR]))))
      sc = np.concatenate((SC, [1]))

      # Least square fitting of parabola to "Taper":
      
      A = np.column_stack((np.concatenate(([np.sum(sc*Taper[:,0]**4)], [np.sum(sc*Taper[:,0]**2)])), np.concatenate(([np.sum(sc*Taper[:,0]**2)], [np.sum(sc)]))))

      y = np.concatenate(([np.sum(sc*Taper[:,1]*Taper[:,0]**2)],[np.sum(sc*Taper[:,1])]))

      x = np.linalg.solve(A,y)

      x[0] = np.min((x[0], -0.0001)) # tapering from the base to the tip
      Ru = x[0]*L**2+x[1]
      Ru[Ru < MinR] = MinR
      if np.max(Ru) > MaxR:
        a = np.max(Ru)
        Ru = MaxR/a*Ru
      Rl = 0.75*(Ru) # lower bound parabola
      Rl[Rl < MinR] = MinR

      # Modify radii based on parabola:
      # change values larger than the parabola-values when SC < 70%:
      I = ((cyl['radius'] > Ru) & (SC < 0.7))
      cyl['radius'][I] = Ru[I] + (cyl['radius'][I]-Ru[I])*SC[I]/0.7
      Mod[I] = True
      # change values larger than the parabola-values when SC > 70% and 
      # radius is over 33% larger than the parabola value:
      I =  ((cyl['radius'] > 1.333*Ru) & (SC >= 0.7))
      cyl['radius'][I] = Ru[I]+(cyl['radius'][I] - Ru[I])*SC[I]
      Mod[I] = True
      # change values smaller than the downscaled parabola-values
      I = ((cyl['radius'] < Rl) & ( SC < 0.7) | (cyl['radius']<0.5*Rl))
      cyl['radius'][I] = Rl[I]
      Mod[I] = True
    else:
      ## Adjust radii of short branches to be lenearly decreasing
      R = cyl['radius']
      if np.count_nonzero(SC >= 0.7) > 1:
        a = np.max(R[SC>=0.7])
        b = np.min(R[SC>=0.7])
      elif np.count_nonzero(SC>= 0.7) == 1:
        a = np.max(R[SC>=0.7])
        b = np.min(R)
      else:
        a = np.sum(R*SC/np.sum(SC))
        b = np.min(R)
      Ru = np.transpose(np.linspace(a, b, numC))
      I = ((SC < 0.7) & ~np.atleast_1d(Mod))
      cyl['radius'][I] = Ru[I] + (R[I] - Ru[I])*SC[I]/0.7
      Mod[I] = True
  

  ## Modify starting points by potimising them for given radius and axis
  numR = len(Regs)
  for i in range(numC):
    if Mod[i]:
      if numR == numC:
        Reg = Regs[i]
      elif i > 0:
        Reg = Regs[i-1]
      cyl['start'] = np.atleast_2d(cyl['start'])
      cyl['axis'] = np.atleast_2d(cyl['axis'])
      if np.abs(cyl['radius'][i] - np.atleast_1d(cyl['radius0'])[i]) > 0.005 and ((numR == numC) or ((numR < numC) and (i > 0))):
        P = Reg - cyl['start'][i,:]
        U, V = orthonormal_vectors(cyl['axis'][i,:])
        P = P@ np.column_stack((U, V))
        cir = least_squares_circle_centre(P, np.array([0,0]), cyl['radius'][i])
        if cir['conv'] and cir['rel']:
          cyl['start'][i,:] = cyl['start'][i,:] + cir['point'][0]*np.transpose(U) + cir['point'][1]*np.transpose(V)
          cyl['mad'][i] = cir['mad']

          d_, V, h, B_ = distances_to_line(Reg, cyl['axis'][i,:], cyl['start'][i,:])
          if np.min(h) < -0.001:
            cyl['length'][i] = np.max(h) - np.min(h)
            cyl['start'][i,:] = cyl['start'][i,:]+np.min(h)*cyl['axis'][i,:]
            d_, V, h, B_ = distances_to_line(Reg, cyl['axis'][i,:], cyl['start'][i,:])
          
          a = np.maximum(0.02, 0.2*cyl['radius'])
          numL = int(np.ceil(cyl['length'][i]/a))
          numL = np.max((numL, 4))
          numS = int(np.ceil(2*np.pi*cyl['radius'][i]/a))
          numS = np.max((numS, 10))
          numS = np.min((numS, 36))
          cyl['SurfCov'][i] = surface_coverage2(cyl['axis'][i,:], cyl['length'][i],V,h,numL, numS)
  

  ## Continous branches
  # Make cylinders properly "continuous" by moving the starting points
  # Move the starting point to the plane defined by parent cylinder's top
  if numC > 1:

    for j in range(1,numC):
      U = cyl['start'][j,:]-cyl['start'][j-1,:]-cyl['length'][j-1]*cyl['axis'][j-1,:]
      if np.linalg.norm(U) > 0.0001:
        # First define vector V and W which are orthogonal to the 
        # cylinder axis N
        N = np.transpose(cyl['axis'][j,:])
        if np.linalg.norm(N) > 0:
          V, W = orthonormal_vectors(N)
          # Now define the new starting point
          x = np.linalg.solve(np.column_stack((N,V,W)), np.transpose(U))
          cyl['start'][j,:] = cyl['start'][j,:]-x[0]*np.transpose(N)
          if x[0] > 0:
            cyl['length'][j] = cyl['length'][j]+x[0]
          elif (cyl['length'][j] + x[0]) > 0:
            cyl['length'][j] = cyl['length'][j]+x[0]

  ## Connect far away first cylinder to the parent
  if len(parcyl['radius']) >0:
    d, V, h, B = distances_to_line(cyl['start'][0,:], parcyl['axis'], parcyl['start'])
    d =  d-parcyl['radius']
    if np.all(d > 0.001):
      taper = cyl['start'][0,:]
      E = taper+cyl['length'][0]*cyl['axis'][0,:]
      V = parcyl['radius']*V/np.linalg.norm(V)
      if (h >= 0) and (h<= parcyl['length']):
        cyl['start'][0,:] = parcyl['start']+B+V
      elif h < 0:
        cyl['start'][0,:] = parcyl['start']+V
      else:
        cyl['start'][0,:] = parcyl['start']+parcyl['length']*parcyl['axis']+V
      cyl['axis'][0,:] = E = cyl['start'][1,:]
      cyl['length'][0] = np.linalg.norm(cyl['axis'][0,:])
      cyl['axis'][0,:] = cyl['axis'][0,:] /cyl['length'][0]

  return cyl

