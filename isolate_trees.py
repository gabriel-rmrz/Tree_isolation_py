DEBUG=False
import numpy as  np
import yaml
from plotting.plot_point_cloud import plot_point_cloud
from plotting.plot_segs import plot_segs
from plotting.plot_tree_structure import plot_tree_structure
from filtering_plot import filtering_plot
from cover_sets_plot import cover_sets_plot
from tools.connected_components import connected_components
from shortest_paths_height import shortest_paths_height
from shortest_paths import shortest_paths
from segments_num_path import segments_num_path
from pydict_to_matstruct import pydict_to_matstruct

import scipy.io

def isolate_trees(P, Hei=None, cover=None):
  '''
  This function implements the tree isolation method in Greco et al. 2023
  "Methodology for measuring dendrometric parameters in a mediterranean
  forest with UAVs flying inside forest" in International Journal of
  Applied Earth Observation and Geoinformation.
 
  Copyright (C) 2020-2023 Pasi Raumonen
 
  Inputs:
  P       The plot level point cloud, (n_points x 3)-matrix
 
  Outputs:
  Trees   Cell array where each cell contains the point cloud of an
            isolated tree
 
  The method has the following steps:
  1. Filter the point cloud and define the height from the ground level
  2. Cover the filtered point cloud with pathces (clustering)
  3. Locate candidate stem sections
  4. Determine the shortest paths to the stem sections
  5. Define the trees as those sets whose paths end up at the same StemSec
  6. Extend the trees downwards with shortest paths
  7. Remove trees that are not connected to the ground or are too short
  8. (Optional) Segment the trees into stem and branches based on shortest paths
  '''

  # Import input dictionary
  with open('configs/inputs.yaml', 'r') as  file:
    inputs = yaml.safe_load(file)

  inputs = inputs['inputs']
  
  ## 1. Filter the point cloud and define the height of the points
  if Hei==None:
    '''
    Pass = Logical vector defining the points passing the filtering
    Hei = Vector giving the heights (in cm) of the points passing the 
    filtering from the ground
    '''
    Pass, Hei, Ground, Tri = filtering_plot(np.copy(P),inputs)
    P = P[Pass,:]

  ## 2. Cover the point cloud with patches
  if cover==None:
    cover  = cover_sets_plot(np.copy(P),inputs)


  ## 3. Locate three candidates by selecting stem sections
  print('\t -------------------------------------')
  print('\t Selecting stem sections...')
  print('\t -------------------------------------')
  # Determine connected components of patches with at least 25 patches
  # between 2 and 3 metters above the ground
  H = np.asarray(Hei[cover['center']]) # Height of the patches
  Sec = (H > 200) *(H < 300) # Logical vector of patches

  StemSec, CompSize = connected_components(cover['neighbor'], Sec, 25)
  
  StemSec_k = {f"k{key+1}":key+1 for key in StemSec.keys()}
  StemSec2 = {f"k{key+1}": value+1 for key, value in StemSec.items()}
  pydict_to_matstruct(StemSec2, 'matlab/StemSec.mat','StemSec')
  pydict_to_matstruct(StemSec_k, 'matlab/StemSec_k.mat','StemSec_k')

  n = len(StemSec)
  print(f"{n} stem sections determined")
  ## NOTICE! Here we have another more sophisticated way to select the stem sections:
  # StemSec = define_stem_sections(P,cover,Hei)
  
  # Plot the stem sections
  Ce = P[cover['center'],:]
  sec = (H > 50) *(H < 400) 
  plot_point_cloud(Ce[sec,:], "centers_50_400")
  plot_segs(P, StemSec,5, cover['ball'])

  ## 4. Determine the shortest path of the patches to the "StemSec"
  print('\t -------------------------------------')
  print('\t Computing the shortest paths...')
  print('\t -------------------------------------')
  # Compute the shortest paths, first with the height restrictions.
  # The proceess also identifies gaps in the grah and creates new links
  H = (Hei[cover['center']]/100).astype(np.double) # The height of the patches in meters
  # Define the stem sections as the base/starting
            # point for the shortest paths
  '''
  Base = [] 
  for key in StemSec.keys():
    Base = np.concatenate((Base, StemSec[key]), axis=0)
  '''
  Base = np.concatenate([StemSec[a] for a in StemSec.keys()])
  

  BaseDist = (Hei[Base]/100).astype(np.double) # Path length at the base is the height

  # Restriction and iteration parameters  for the shortest path computation
  inputs['GD0'] = 0.1
  inputs['DHRel0'] = 1.25
  inputs['GD'] = 0.05
  inputs['DHRel'] = 0.25
  inputs['MaxDHRel'] = 1.75
  #shortest_paths_height return  PathNum, PathDist, EndSet
  EndSet, cover, PathLen = shortest_paths_height(np.copy(P), cover, np.copy(Hei), np.copy(Base), np.copy(BaseDist), inputs)

  ## 5. Define the trees as those sets whose path end up at the same stem section
  print('  -----------------------')
  print('   Defining the trees...')
  print('  -----------------------')
  numT = len(StemSec) # number of stem setions
  numB = len(cover["ball"]) # number of sets/patches
  ind = np.array(range(numB))
  Trees = {} # Initialize the ouput dictionary
  Base = {} # tree bases
  for i in range(numT):
    S = StemSec[i]
    n = len(S)
    C = {}
    for j in range(n):
      C[j] = ind[(EndSet == S[j]) & (EndSet>0)]
    T = np.concatenate([C[key] for key in C.keys()])
    Trees[i] = T
    Base[i] = T[H[T] < (min(H[T]) + 0.5)]
  ## 6. Extend the trees downwards with shortest paths
  # Forbiden sets and the base for the paths:
  Forb = np.zeros(numB, dtype='bool') # Sets already assignd to trees
  base = np.concatenate([Base[key] for key in Base.keys()])
  Forb[np.concatenate([Trees[key] for key in Trees.keys()])] = True
  Forb[base] = False
  # Determine the shortest paths (assumes the trees are connected to the gound level)
  _, PathLen, EndSet = shortest_paths(cover, np.copy(base), np.copy(Forb))

  if DEBUG:
    exit()

  ## Define the bottoms of the trees as thos sets whose path_lenght/height
  #  ratio is less than 1.1
  for i in range(numT):
    B = Base[i] # the base sets of the tree "i"
    n = len(B)
    C = {} # each element contains all the set connected to the corresponding base set
    for j in range(n):
      S = ind[EndSet == B[j]] # sets connected to the base set B[j-1]
      I = PathLen[S] < 1.1*H[B[j]] # paths are straight enough
      C[j] = S[I]
    Bottom = np.concatenate([C[key] for key in C.keys()]) # The bottom sets
    Trees[i] = np.concatenate([Trees[i], Bottom]) # include the bottom to the tree
  ## 7. Remove trees that are not connected to the ground of are too short
  # Keep only the trees that are close enough  to the ground level (minimum 
  # height is less than 50 cm) and are tall enough (minimum height is 5 m)
  numT = len(Trees) # number of trees
  H = Hei[cover["center"]] # height of the sets
  Keep = np.ones(numT, dtype='bool') # keep these threes
  Bases = {} # Tree bases
  Ind = np.array(range(numT)).astype(np.uint32)
  for i in range(numT):
    T = Trees[i]
    if (len(T) ==0) or (np.min(H[T]) > 50) or (np.max(H[T]) - np.min(H[T])) < 500:
      Keep[i] = False
    else:
      I = H[T] < 50
      Bases[i] = T[I]
  Trees = {key: Trees[key] for key in Ind[Keep] }
  numT = len(Trees)
  Bases = {key: Bases[key] for key in Ind[Keep] }


  ## 8. Segment the trees into stem and branches based on shortest paths
  # Determine the shortest paths
  base = np.concatenate([Bases[key] for key in Bases.keys()])
  Forb = np.ones(numB, dtype='bool')
  Forb[np.concatenate([Trees[key] for key in Trees.keys()])] = False
  PathNum, _1, _2  = shortest_paths(cover, np.copy(base), np.copy(Forb))
  # Segment each tree and select only the stem
  Segments = {}
  for i in range(numT):
    if not(i in Trees.keys()):
      continue
    Forb = np.ones(numB, dtype='bool')
    Forb[Trees[i]] = False
    segment =segments_num_path(cover, Bases[i],Forb,PathNum)
    Segments[f"{i}"] = segment
    plot_segs(P,{0:Trees[i]},1,cover["ball"],f"tree_{i}")
    #plot_segs(P,segment['segments'][0],1,cover["ball"],f"segs_{i}")
    plot_tree_structure(P,cover,segment,1,10,0, f"segs_{i}")

  return P, inputs,  cover, Trees, Segments

  

