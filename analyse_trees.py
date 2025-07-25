DEBUG=True
import json
import time 
import numpy as np
def cylinders(P, cover, segment, inputs):
  '''
  Computes the following quantities for every segment:
  - radius
  - length
  - axis
  - start
  '''
  print(segment)
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

  if DEBUG:
    print(len(Trees[0]))
    print(f"cover.keys(): {cover.keys()}")
    print(f"cover['inputs'].keys(): {cover['inputs'].keys()}")
    #print(f"cover): {cover}")
    exit()
  i_time = time.time() # Initial running time
  Diam = np.zeros([numT,16], dtype=np.int32)
  for i in range(numT):
    if Model[i] and len(Trees[i]) > 0:
      segment = {}

      segment['segments'] = Trees[i] # TODO: From this we will only take into account stem segments.
      segment['ParentSegment'] = 0
      segment['ChildSegment'] = {}
      cyls = cylinders(P, cover, segment, inputs)


  f_time = time.time()

  print("TEST message")
