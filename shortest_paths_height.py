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




  EndSet = 0
  PathLen = 0
  return  EndSet, cover, PathLen

