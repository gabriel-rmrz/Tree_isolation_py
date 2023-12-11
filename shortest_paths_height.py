import numpy as np
from collections import defaultdict

def shortest_paths_height(P, cover, Hei, Base, BaseDist, inputs, Forb=None):
  # Determines the shortest paths from every point to the base

  nei = cover['neighbor'] # the neighbors of the patches
  NeiDis = cover['NeiDis'] # distances between the neighbors
  numB = len(cover['ball']) # number of patches
  Ce = P[cover['center'],:] # the centers of the patches
  ind = np.asarray(range(numB)).astype(int) # indices from 0 to numB-1
  Hei = np.asarray(Hei[cover['center']]/100).astype(float) # the heights of the patches (in m)
  GD = inputs['GD0'] # Intial gap distance that can be brigded over with a new link
  DHRel = inputs['DHRel0'] # Initial maximum allowed path_lenght/height relation

  # If not forbidden patches given, then define all patches allowable
  if Forb ==None:
    For = np.zeros(numB, dtype='bool')
  
  ## Determine the shoretest paths for as many sets as possible.
  # Paths are not computed for the sets in different components than the base
  # sets and sets whose paths distance compared to height is to large
  # Sort the base according to their distance:

  I = np.argsort(BaseDist)
  BaseDist = BaseDist[I]
  Base = Base[I]
  C = Base[0] # The current set whose patch distance is to be determined
  Unvisited = np.zeros(numB, dtype=np.uint32) # patchees/sets not yet visited by paths
  n = len(Base)
  Unvisited[:n] = Base # Base set have been visited
  UV = np.zeros(numB, dtype='bool') # Unvisited sets in logical vercor also
  UV[Base] = True
  a = 1 # Counts how many unvisited sets have been visited by a shortest path
  b = n # how many unvisited sets there currently are
  J = 1 # index of the set C (a+J-1)
  PathLen = 1000*np.ones(numB,dtype=np.single)
  PathNei = np.zeros(numB, dtype=np.uint32)
  EndSet = np.zeros(numB, dtype=np.uint32)
  PathLen[Base] = BaseDist
  EndSet[Base] = Base
  while a <= b:
    nei[C] # neighbors of the set C
    d = NeiDis[C] # the distances of the the neighbors of C
    L = PathLen[C] + d # Path lenght to the neighbors of C
    # Accept the path via C to N if the path lenghts are shoreter, if the 
    # length/height ratio is small enough, if the sets are not forbidden.
    I = L < PathLen[N] & L / Hei[N] < DHRel & ~Forb[N]
    N = N[I]
    if len(N)>0:
      # Update the length, neighbor and endset for the sets N:
      PathLen[N] = L[I]
      PathNei[N] = C
      EndSet[N] = EndSet[C]

    # Define C as visited and N as unvisited
    if J>1: 
      Unvisited[a-J-1-1] = Unvisited[a-1]

    a += 1
    UV[C] = False
    if len(N) > 0:
      I = UV[N]
      N = N[~I]
      if len(N) > 0:
        Unvisited[b+1-1: b+len(N)] = N
        b = b + len(N)
        UV[N] = True

    # Determine the next set "C" form the unvisited sets that has the minimum
    # path length
    _, J = np.argmin(PathLen[Unvisited[a-1:b]])
    C = Unvisited[a + J -1 -1)
  
  ## Iteratively expand trees base on the shortest paths
  # Partition of the cover sets into search space

  Partition, CC, info = cubical_partition(Ce, GD, nargout=3)
  




  EndSet = 0
  PathLen = 0
  return  EndSet, cover, PathLen

