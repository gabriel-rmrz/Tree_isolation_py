DEBUG = False
import json
import time 
import numpy as np
from plotting.plot_point_cloud import plot_point_cloud as plot_pc
def cylinders(P, cover, Segments, inputs):
  '''
  Computes the following quantities for every segment:
  - radius
  - length
  - axis
  - start
  '''

  for key in Segments.keys():
    trunk = []
    for ind, seg in Segments[key]['segments'][0].items():
      trunk = np.concatenate((trunk, seg), axis=0)
    trunk_points = []

    for t in trunk:
      print(f"cover['center'][t.astype(int)]: {cover['center'][t.astype(int)]}")
      print(f"t: {t}")
      trunk_points.append(cover['center'][t.astype(int)])
    #trunk_points=np.concatenate(trunk_points)
    print(trunk_points)
    plot_pc(P[trunk_points], "segment_0_tree_0")
  print(f"cover.keys(): {cover.keys()}")
  return 0

def analyse_trees(P, inputs, cover, Trees, Segments):
  '''
  analyse_tree.py       Computes volumes and diameters of the trees.
  
  Version 2.0.0
  Latest update       18 Dec 2019

  Copyright (C) 2015-2019 Pasi Raumonen
  ---------------------------------------------------------------------
  Analysing in the following steps:
  1.    (Optional?) Something

  Details of each step:
  1.

  Inputs:
  P             Plot level point cloud.
  inputs Dictionary: Needs the following fiels:
    1.
  cover,
  Trees
  '''

  print('  -----------------------')
  print('   Analysing trees...')
  print('  -----------------------')

  # This is an input in the original script.
  numT = len(Trees)
  # This is the index of the tree that are selected to be analysed. All of them are selected with this configuration.
  Model = np.ones(numT, dtype=bool)
  # TODO: Here are defined some variable for the field measurements. We will only work with the PointCloud measurement for the moment.
  if DEBUG:
    print("-------------------------------------------------------------")
    print(f"Segments.keys(): {Segments.keys()}") # TEST: 0 and 1 (we have two trees)
    print(f"Segments['0'].keys(): {Segments['0'].keys()}")
    print("-------------------------------------------------------------")
    #Segments['0'].keys(): dict_keys(['segments', 'parent', 'children'])
    print(f"type(Segments['0']['segments']): {type(Segments['0']['segments'])}")
    print(f"type(Segments['0']['parent']): {type(Segments['0']['parent'])}")
    print(f"type(Segments['0']['children']): {type(Segments['0']['children'])}")
    print("-------------------------------------------------------------")
    print(f"Segments['0']['segments'].keys(): {Segments['0']['segments'].keys()}")
    print(f"Segments['0']['parent']: {Segments['0']['parent']}")
    print(f"Segments['0']['children'].keys(): {Segments['0']['children'].keys()}")
    print("-------------------------------------------------------------")
    for ind, seg in Segments['0']['segments'].items():
      print(f"{ind}, len(seg): {len(seg)}")
    print("-------------------------------------------------------------")
    trunk = []
    for ind, seg in Segments['0']['segments'][0].items():
      print(f"{ind}, len(seg): {len(seg)}")
      trunk = np.concatenate((trunk, seg), axis=0)

    trunk_points = []
    for t in trunk:
      trunk_points.append(cover['ball'][t])
    trunk_points=np.concatenate(trunk_points)
    print(f"P[trunk_points]: {P[trunk_points]}")
    plot_pc(P[trunk_points], "segment_0_tree_0")
    exit()
    print(f"trunk: {trunk}")

    exit()
  col = [
    [0.00, 0.00, 1.00],
    [0.00, 0.50, 0.00],
    [1.00, 0.00, 0.00],
    [0.00, 0.75, 0.75],
    [0.75, 0.00, 0.75],
    [0.75, 0.75, 0.00],
    [0.25, 0.25, 0.25],
    [0.75, 0.25, 0.25],
    [0.95, 0.95, 0.00],
    [0.25, 0.25, 0.75],
    [0.75, 0.75, 0.75],
    [0.00, 1.00, 0.00],
    [0.76, 0.57, 0.17],
    [0.54, 0.63, 0.22],
    [0.34, 0.57, 0.92],
    [1.00, 0.10, 0.60],
    [0.88, 0.75, 0.73],
    [0.10, 0.49, 0.47],
    [0.66, 0.34, 0.65],
    [0.99, 0.41, 0.23]
  ]

  i_time = time.time() # Initial running time
  Diam = np.zeros([numT,16], dtype=np.int32)
  print(f"We have {len(Segments)} trees.")
  cyls = cylinders(P, cover, Segments, inputs)


  f_time = time.time()

