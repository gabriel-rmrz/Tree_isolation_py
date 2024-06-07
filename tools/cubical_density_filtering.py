import numpy as np
from collections import defaultdict

def cubical_density_filtering(P,Points, inputs):
  print('---------')
  print('Filtering based on cubical point density...')
  # Pararmeters for cubical partition

  CubeSize = inputs['CubeSizeDensity']
  numP = len(P[:,0])
  Min = P.min(axis=0).astype(np.double)
  Max = P.max(axis=0).astype(np.double)
  
  N = (np.ceil((Max-Min)/CubeSize)).astype(np.double)

  if np.floor(N[0]) == np.ceil(N[0]):
    N[0] = N[0] + 1
  else:
    N[0] = np.ceil(N[0])
  if np.floor(N[1]) == np.ceil(N[1]):
    N[1] = N[1] + 1
  else:
    N[1] = np.ceil(N[1])
  
  # Compute the cube coordinates
  C = np.floor((P-Min)/CubeSize)+1
  LexOrd = C[:,0] + (C[:,1] - 1)*N[0] + (C[:,2] -1)*N[0]*N[1]

  # Sort the coordinates
  Ind = np.argsort(LexOrd)
  LexOrd = np.sort(LexOrd)

  # Definite inverse sort mapping "ind"
  
  ind = np.zeros(numP, dtype=int)
  for i in range(numP):
    ind[Ind[i]] = i

  # Check the cubes for number of points
  Pass = np.ones(numP, dtype='bool')
  p = 1
  while p <= numP:
    t = 1
    while p+t <= numP and LexOrd[p-1+t] == LexOrd[p-1]:
      t = t+1
    if t < inputs["PointDensity"]: # too few points, filter all
      Pass[p-1:p+t-1] = False
    p = p + t

  Pass = Pass[ind]
  P = P[Pass,:]

  numP2 = len(P[:,0])
  numP0 = len(Points)
  Ind = np.asarray(range(numP0))
  Ind = Ind[Points]
  Points[Ind[~Pass]] = False

  print(f"\t Points before: {numP}")
  print(f"\t Filtered points: {numP-numP2}")
  print(f"\t Points left: {numP2}")

  return P, Points
