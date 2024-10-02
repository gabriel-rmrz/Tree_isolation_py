import numpy as np
from collections import defaultdict

def remove_lowest_points(P, Points):
  print('---------')
  print('Remove the lowest points...')

  # Remove (some of) the ghost points under the ground level
  ## Basic parameters
  numP = len(P[:,0])
  Min = P.min(axis=0).astype(np.double)
  Max = P.max(axis=0).astype(np.double)
  
  # TODO: Add the SQ parameter to the inputs dictionary in the yaml file.
  SQ = 5. # side lengh of the rectangles for ghost point search
  Nx = (Max[0] - Min[0])/SQ # Number of rectangles in the x-y directions
  Ny = (Max[1] - Min[1])/SQ # Number of rectangles in the x-y directions
  if np.floor(Nx) == np.ceil(Nx):
    Nx = Nx + 1
  else:
    Nx = np.ceil(Nx)
  if np.floor(Ny) == np.ceil(Ny):
    Ny = Ny + 1
  else:
    Ny = np.ceil(Ny)
  
  ## Search for the lowest points

  Par = defaultdict(list)
  R = np.floor((P[:,:2] - Min[:2])/SQ) + 1
  LexOrd = R[:,0] + Nx*(R[:,1]-1)
  I = np.argsort(LexOrd)
  R = (R[I,:]).astype(np.int32)
  LexOrd = np.sort(LexOrd)
  PointInd = np.asarray(range(numP), dtype=np.uint32)
  PointInd = PointInd[I]
  #Points = np.asarray(Points[I])
  n=len(R)
  q = 1
  while q <= n:
    t = 1
    while (q+t <= n) and (LexOrd[q-1] == LexOrd[q-1+t]):
      t = t + 1
    points = PointInd[q-1:q+t-1]
    if t > 1e4:
      J = np.argsort(P[points,2])
      points = points[J[:int(1e4)]]
    Par[tuple([R[q-1,0], R[q-1,1]])] = points
    q = q + t
  
  ## Remove the lowest points
  Pass = np.ones(numP, dtype='bool')
  for i in range(int(Nx)):
    for j in range(int(Ny)):
      points = Par[tuple([i,j])]
      if np.any(points):
        N2, Edges = np.histogramdd(P[points,2], 100)

        CS = np.cumsum(N2)/len(points)*100
        k=0
        while CS[k] < 1:
          k = k+1
        if k > 0:
          k = k - 1
        Pass[points] = P[points,2] > Edges[0][k]

  numP0 = len(Points)
  Ind = np.asarray(range(numP0), dtype=np.uint32)
  Ind = Ind[Points]
  Points[Ind[~Pass]] = False
  P = P[Pass,:]
  numP2 = len(P[:,0])
  print(f"\t Points before: {numP}")
  print(f"\t Filtered points: {numP-numP2}")
  print(f"\t Points left: {numP2}")



  return P, Points
