DEBUG = True
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

def test_isolate_trees(P, Hei=None, cover=None):
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
  print(f"Ground[0]: {Ground[0]}")
  print(f"Tri[0]: {Tri[0]}")
  
  print(f"type(Pass[0]): {type(Pass[0])}")
  print(f"type(Ground[0]): {type(Ground[0])}")
  print(f"type(Hei[0]): {type(Hei[0])}")
  print(f"type(Tri[0]): {type(Tri[0])}")
  scipy.io.savemat('matlab/Pass.mat', {'Pass': Pass})
  scipy.io.savemat('matlab/Hei.mat', {'Hei': Hei})
  scipy.io.savemat('matlab/Ground.mat', {'Ground': Ground})
  scipy.io.savemat('matlab/Tri.mat', {'Tri': Tri})
  scipy.io.savemat('matlab/P.mat', {'P': P})
  pydict_to_matstruct(inputs, 'matlab/inputs.mat','inputs')

  exit()
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
  if DEBUG:
    #print(f"H: {H}")
    #print(f"H: {H}")
    #print(f"sum(Sec): {sum(Sec)}")
    print(f"len(Sec): {len(Sec)}")
    print(f"len(cover.center): {len(cover['center'])}")
    print(f"len(cover.neighbor): {len(cover['neighbor'])}")
    #print(sum(Sec))

  StemSec, CompSize = connected_components(cover['neighbor'], Sec, 25)

  if DEBUG:
    #print(f"StemSec[1]: {StemSec[1]}")
    #print(f"StemSec.keys(): {StemSec.keys()}")
    #plot_point_cloud(P[cover['center'][np.concatenate([StemSec[0], StemSec[1], StemSec[2], StemSec[3], StemSec[4], StemSec[5]])],:])
    plot_point_cloud(P[cover['center'][np.concatenate([StemSec[0], StemSec[1]])],:], "two_stems")
    #plot_point_cloud(P[cover['center'][StemSec[0]],:])
    #plot_point_cloud(P[cover['center'][Sec,:],:])
    #plot_point_cloud(P)
    plot_segs(P,{0:StemSec[0], 1:StemSec[1]},1,cover["ball"], 'two_stems')
    #print(StemSec)
    #print(CompSize)
    print(f"len(CompSize): {len(CompSize)}")
  n = len(StemSec)
  print(f"{n} stem sections determined")
  ## NOTICE! Here we have another more sophisticated way to select the stem sections:
  # StemSec = define_stem_sections(P,cover,Hei)
  
  # Plot the stem sections
  Ce = P[cover['center'],:]
  sec = (H > 50) *(H < 400) 
  plot_point_cloud(Ce[sec,:], "centers_50_400")
  #plot_segs(P, StemSec,5, cover['ball'])

  ## 4. Determine the shortest path of the patched to the "StemSec"
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

  # Determine the shortest paths
  #shortest_paths_height return  PathNum, PathDist, EndSet
  if DEBUG:
    print(f"type(P): {type(P)}")
    print(f"type(Hei): {type(Hei)}")
    print(f"type(Base): {type(Base)}")
    print(f"type(BaseDist): {type(BaseDist)}")
  EndSet, cover, PathLen = shortest_paths_height(np.copy(P), cover, np.copy(Hei), np.copy(Base), np.copy(BaseDist), inputs)

  if DEBUG:
    print(f"len(EndSet): {len(EndSet)}")
    print(f"len(cover): {len(cover)}")
    print(f"len(PathLen): {len(PathLen)}")
    #print(f"cover['neighbor']: {cover['neighbor']}")
    #print(CompSize)
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
      C[j] = ind[EndSet == S[j]]
    T = np.concatenate([C[key] for key in C.keys()])
    Trees[i] = T
    Base[i] = T[H[T] < (min(H[T]) + 0.5)]
  if DEBUG:
    plot_point_cloud(P[cover['center'][Trees[0]],:], "0_tree_p5")
    plot_point_cloud(P[cover['center'][Trees[1]],:], "1_tree_p5")
    plot_point_cloud(P[cover['center'][Trees[2]],:], "2_tree_p5")
    plot_point_cloud(P[cover['center'][Trees[3]],:], "3_tree_p5")
    plot_point_cloud(P[cover['center'][Trees[4]],:], "4_tree_p5")
    plot_point_cloud(P[cover['center'][Trees[5]],:], "5_tree_p5")
    #plot_segs(P,Trees,1,cover["ball"],'segs')
    print(f"5: len(Trees): {len(Trees)}")


  ## 6. Extend the trees downwards with shortest paths
  # Forbiden sets and the base for the paths:
  Forb = np.zeros(numB, dtype='bool') # Sets already assignd to trees
  base = np.concatenate([Base[key] for key in Base.keys()])
  Forb[np.concatenate([Trees[key] for key in Trees.keys()])] = True
  Forb[base] = False

  # Determine the shortest paths (assumes the trees are connected to the gound level)
  _, PathLen, EndSet = shortest_paths(cover, np.copy(base), np.copy(Forb))

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
  if DEBUG:
    plot_point_cloud(P[cover['center'][np.concatenate([Trees[0], Trees[1], Trees[2], Trees[3], Trees[4], Trees[5]])],:], "six_stems")
    plot_point_cloud(P[cover['center'][Trees[1]],:], "1_tree_p6")
    plot_point_cloud(P[cover['center'][Trees[0]],:], "0_tree_p6")
    plot_point_cloud(P[cover['center'][Base[1]],:], "1_base_p6")
    plot_point_cloud(P[cover['center'][Base[0]],:], "0_base_p6")
    print(f"6: len(Trees): {len(Trees)}")
  
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
    if (len(T) ==0) or (np.min(H[T]) > 100) or (np.max(H[T]) - np.min(H[T])) < 500:
      Keep[i] = False
    else:
      I = H[T] < 50
      Bases[i] = T[I]
  Trees = {key: Trees[key] for key in Ind[Keep] }
  numT = len(Trees)
  Bases = {key: Bases[key] for key in Ind[Keep] }

  ## 8. Segment the trees into stem and branches based on shortest paths
  # Determine the shortest paths
  base = np.concatenate([Bases[key] for key in Bases.keys()]) # the bases of the paths
  Forb = np.ones(numB, dtype='bool') # the forbiden sets for the paths
  Forb[np.concatenate([Trees[key] for key in Trees.keys()])] = False
  PathNum, PathDist, EndSet  = shortest_paths(cover, np.copy(base), np.copy(Forb))
  if DEBUG:
    print(f"PathNum: {PathNum}")
    print(f"len(PathNum): {len(PathNum)}")
    print(f"len(cover['ball']): {len(cover['ball'])}")
    print(f"type(cover['ball']): {type(cover['ball'])}")
    print(f"len(cover['center']): {len(cover['center'])}")
    print(f"type(cover['center']): {type(cover['center'])}")
    print(f"len(cover['neighbor']): {len(cover['neighbor'])}")
    print(f"type(cover['neighbor']): {type(cover['neighbor'])}")
    print(f"len(cover['NeiDis']): {len(cover['NeiDis'])}")
    print(f"type(cover['NeiDis']): {type(cover['NeiDis'])}")
    print(f"len(cover['BallOfPoint']): {len(cover['BallOfPoint'])}")
    print(f"type(cover['BallOfPoint']): {type(cover['BallOfPoint'])}")
    print(f"len(cover['inputs']): {len(cover['inputs'])}")
    print(f"type(cover['inputs']): {type(cover['inputs'])}")
    print(f"cover['inputs']: {cover['inputs']}")


  # Segment each tree and select only the stem
  for i in range(numT):
    Forb = np.ones(numB, dtype='bool')
    Forb[Trees[i]] = False
    segment =segments_num_path(cover, Bases[i],Forb,PathNum)
    plot_segs(P,{0:Trees[i]},1,cover["ball"],f"tree_{i}")
    print(f"segment: {segment}")
    #plot_segs(P,segment['segments'][0],1,cover["ball"],f"segs_{i}")
    plot_tree_structure(P,cover,segment,1,10,0, f"segs_{i}")


  

