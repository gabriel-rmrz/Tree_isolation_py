DEBUG=False
import numpy as np
def cubical_down_sampling(P,Points,inputs):
  print('---------')
  print('Cubical downsampling...')
  # Downsamples the given point cloud by selecting one point for each
  # cube od side length CubeSize

  # The vertices of the big cube containing P
  
  Min = P.min(axis=0).astype(np.double)
  Max = P.max(axis=0).astype(np.double)
  if DEBUG:
    print(f"Min: {Min}")
    print(f"Max: {Max}")
  numP0 = len(P[:,0])

  # Number of cubes with edge length "EdgeLength" in the sides
  # of the big cube

  Passed = np.zeros(numP0, dtype="bool")
  N = (np.ceil((Max-Min)/inputs["SamplingCubeSize"])+1).astype(np.double)
  if DEBUG:
    print(f"N: {N}")

  # Compute the cube coordinates of the points
  C = np.floor((P -Min)/inputs["SamplingCubeSize"])+1.
  # Compute the lexicographical order of the cubes
  LexOrd = C[:,0] + (C[:,1] -1)*N[0] + (C[:,2] - 1) * N[0]*N[1]
  LexOrd, I = np.unique(LexOrd, return_index=True, axis=0)
  Passed[I] = True 
  P= P[Passed]
  Ind = np.asarray(range(numP0)).astype(np.uint32)
  Ind = Ind[Points]
  Points[Ind[~Passed]] = False
  numP = len(P[:,0])


  print(f"\t Points before: {numP0}")
  print(f"\t Filtered points: {numP0 - numP}")
  print(f"\t Points left: {numP}")

  if DEBUG:
    print(f"LexOrd: {LexOrd}")
    exit()
  return P, Points
