import numpy as  np
import yaml
from ploting.plot_point_cloud import plot_point_cloud
from filtering_plot import filtering_plot
from cover_sets_plot import cover_sets_plot
from tools.connected_components import connected_components
from shortest_paths_height import shortest_paths_height
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
    Pass, Hei, Ground, Tri = filtering_plot(P,inputs)
    P = P[Pass,:]

  ## 2. Cover the point cloud with patches
  if cover==None:
    cover  = cover_sets_plot(P,inputs)


  ## 3. Locate three candidates by selecting stem sections
  print('\t -------------------------------------')
  print('\t Selecting stem sections...')
  print('\t -------------------------------------')
  # Determine connected components of patches with at least 25 patches
  # between 2 and 3 metters above the ground
  H = np.asarray(Hei[cover['center']]) # Height of the patches
  Sec = (H > 200) *(H < 300) # Logical vector of patches
  #print(f"sum(Sec): {sum(Sec)}")
  #print(f"len(Sec): {len(Sec)}")
  #print(f"len(cover.center): {len(cover['center'])}")
  #print(f"len(cover.neighbor): {len(cover['neighbor'])}")
  #print(sum(Sec))
  exit()

  StemSec, CompSize = connected_components(cover['neighbor'], Sec, 25)
  print(StemSec)
  print(CompSize)
  n = len(StemSec)
  print(f"{n} stem sections determined")
  ## NOTICE! Here we have another more sophisticated way to select the stem sections:
  # StemSec = define_stem_sections(P,cover,Hei)
  
  # Plot the stem sections
  Ce = P[cover['center'],:]
  sec = (H > 50) *(H < 400) 
  plot_point_cloud(Ce[sec,:])
  #plot_segs(P, StemSec,1, 5, cover['ball']))

  ## 4. Determine the shortest path of the patched to the "StemSec"
  print('\t -------------------------------------')
  print('\t Computing the shortest paths...')
  print('\t -------------------------------------')
  # Compute the shortest paths, first with the height restrictions.
  # The proceess also identifies gaps in the grah and creates new links
  H = (Hei[cover['center']]/100).astype(float) # The height of the patches in meters
  Base = [] # Define the stem sections as the base/starting
            # point for the shortest paths
  for key in StemSec.keys():
    Base = np.concatenate((Base, StemSec[key]), axis=0)
  
  BaseDist = (Hei[Base]/100).astype(float) # Path length at the base is the height

  # Restriction and iteration parameters  for the shortest path computation
  inputs['GD0'] = 0.1
  inputs['DHRel0'] = 1.25
  inputs['GD'] = 0.05
  inputs['DHRel'] = 0.25
  inputs['MaxDHRel'] = 1.75

  # Determine the shortest paths
  #shortest_paths_height return  PathNum, PathDist, EndSet
  EndSet, cover, PathLen = shortest_paths_height(P, cover, Hei, Base, BaseDist, inputs)


