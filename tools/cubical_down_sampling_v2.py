import numpy as np
from collections import defaultdict

def cubical_density_filtering(P,Points, inputs):
  print('---------')
  print('Cubical downsampling...')
  # Downsamples the given point cloud by selecting one point from each
  # cube of side length CubeSize

  # The vertice of the big cube containing P

  Min = P.min(axis=0).astype(float)
  Max = P.max(axis=0).astype(float)
  
  # Number of cubes with the edge length "EdgeLength" in the sides
  # of the the big cube
  N = (np.ceil((Max-Min)/inputs["SamplingCubeSize"]) + 1).astype(float)

  numP = len(P[:,0])
  PointInd = np.asarray(range(numP))
  R = defaultdict(list)

  C = np.floor((P-Min)/inputs["SamplingCubeSize"]) + 1
  LexOrd = C[:,0] + (C[:,1] -1)*N[0] + (C[:,2] - 1) * N[1]
  LexOrd, I = np.unique(LexOrd, return_index=True, axis=0)
  PointInd =PointInd[I]

  Pass = np.zeros(numP, dtype='bool')
  Pass[PointInd] = True
  P = P[Pass,:]

  numP2 = len(P[:,0])
  numP0 = len(Points)
  Ind = np.asarray(range(numP))
  Ind = Ind[PointInd]
  Points[Ind[not Pass.all()]] = False

  print(f"\t Points before: {numP}")
  print(f"\t Filtered points: {numP-numP2}")
  print(f"\t Points left: {numP2}")

  return P, Points
