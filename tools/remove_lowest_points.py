DEBUG=False
import numpy as np
from collections import defaultdict

def remove_lowest_points(P, Points):
  print('---------')
  print('Remove the lowest points...')

  # Remove (some of) the ghost points under the ground level
  ## Basic parameters
  numP = len(P[:,0])
  Min = P.min(axis=0).astype(float)
  Max = P.max(axis=0).astype(float)
  
  # TODO: Add the SQ parameter to the inputs dictionary in the yaml file.
  SQ = 5 # side lengh of the rectanglees for ghost point search
  N = (Max - Min)/SQ # Number of rectangles in the x-y directions
  if np.floor(N[0]) == np.ceil(N[0]):
    N[0] = N[0] + 1
  else:
    N[0] = np.ceil(N[0])
  if np.floor(N[1]) == np.ceil(N[1]):
    N[1] = N[1] + 1
  else:
    N[1] = np.ceil(N[1])
  N = N.astype(int)
  if DEBUG:
    print(f"N: {N}")
  
  ## Search for the lowest points

  Par = defaultdict(list)
  R = (np.floor((P[:,:2] - Min[:2])/SQ) + 1).astype(int)
  LexOrd = R[:,0] + N[0]*(R[:,1]-1)
  I = np.argsort(LexOrd)
  if DEBUG:
    print(f"R: {R}")
    print(f"len(R): {len(R)}")
  R = R[I,:]
  if DEBUG:
    print(f"R: {R}")
    print(f"len(R): {len(R)}")
  LexOrd = np.sort(LexOrd)
  PointInd = np.asarray(range(numP))
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
  if DEBUG:
    print(f"######################")
    print(f"N[0]: {N[0]}")
    print(f"range(1,N[0]): {range(1,N[0])}")
    print(f"Par.keys(): {Par.keys()}")
    '''
    for key in Par.keys():
      print(f"len(Par): {len(Par[tuple(key)])}")
    '''
  for i in range(1, N[0] +1):
    for j in range(1, N[1] +1):
      points = Par[tuple([i,j])]
      if DEBUG:
        print(f"pre i: {i}")
        print(f"pre j: {j}")
      if points.any():
        N2, Edges = np.histogramdd(P[points,2], 100)

        #CS = np.cumsum(N)/len(points)*100
        CS = np.cumsum(N2)/len(points)*100
        k=0
        while CS[k] < 1:
          k = k+1
        if k > 0:
          k = k - 1
        Pass[points] = P[points,2] > Edges[0][k]
        if DEBUG:
          #print(f"Edges: {Edges}")
          #print(f"len(Edges[0]): {len(Edges[0])}")
          print(f"i: {i}")
          print(f"j: {j}")
          #print(f"CS: {CS}")
          print(f"len(CS): {len(CS)}")
          print(f"k: {k}")
          print(f"Edges[0][k]: {Edges[0][k]}")
  if DEBUG:
    print(f"sum(Pass): {sum(Pass)}")

  numP0 = len(Points)
  Ind = np.asarray(range(numP0))
  Ind = Ind[Points]
  Points[Ind[~Pass]] = False
  P = P[Pass,:]
  numP2 = len(P[:,0])
  print(f"\t Points before: {numP}")
  print(f"\t Filtered points: {numP-numP2}")
  print(f"\t Points left: {numP2}")

  if DEBUG:
    exit()


  return P, Points
