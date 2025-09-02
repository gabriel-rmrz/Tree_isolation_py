DEBUG = True 
if DEBUG:
  import matplotlib.pyplot as plt
import json
import time 
import numpy as np

from sklearn.covariance import MinCovDet

from plotting.plot_point_cloud import plot_point_cloud as plot_pc
from plotting.plot_slices import plot_slices 

"""
def cylinders(P, cover, Segments, inputs):
  '''
  Computes the following quantities for every segment:
  - radius
  - length
  - axis
  - start
  '''
  cyls = {}
  min_np = 100 # Minimal number of points to run the fitting algorithm
  slice_dh = 0.2 # Size of the cut in z (m) to define the cylinder's height.

  for key in Segments.keys():
    trunk = []
    for ind, seg in Segments[key]['segments'][0].items():
      trunk = np.concatenate((trunk, seg), axis=0)
    trunk_points = []

    '''
    #This would be to take into account only he centers of the balls
    for t in trunk:
      trunk_points.append(cover['center'][t.astype(int)])
    '''
    #This would be to take into account all points in the balls
    for t in trunk:
      trunk_points.append(cover['ball'][t.astype(int)])
    trunk_points=np.concatenate(trunk_points)
    split_edges = np.arange(min(P[trunk_points][:,2]), max(P[trunk_points][:,2]), slice_dh)

    for i in range(len(split_edges)-1):
      bottom_edge = split_edges[i]
      top_edge = split_edges[i+1]
      slice_points = trunk_points[P[trunk_points][:,2] > bottom_edge]
      slice_points = slice_points[P[slice_points][:,2] < top_edge]
      print(f"{bottom_edge} {top_edge}")
      print(f"P[slice_points].shape: {P[slice_points].shape}")
      U, S, Vh = np.linalg.svd(P[slice_points])
      print(f"U.shape: {U.shape}")
      print(f"S.shape: {S.shape}")
      print(f"S: {S}")
      print(f"Vh.shape: {Vh.shape}")
      print(f"Vh: {Vh}")
      #print(f"P: {P}")
      xv = [] 
      '''
      xv_mcd = []
      for x in P[slice_points]:
        xdotv = []
        for v in Vh:
          xdotv.append(np.dot(x,v))
        xv.append(xdotv)
        xdotv = np.array(xdotv)
        xv_mcd.append(MinCovDet(random_state=0).fit(xdotv.reshape(1, -1)))
        print(f"MinCovDet(random_state=0).fit(xdotv.reshape(-1, 1)).location_: {MinCovDet(random_state=0).fit(xdotv.reshape(-1, 1)).location_}")
        print(f"xdotv: {xdotv}")
        exit()
      print(f"xv[0]: {xv[0]}")
      print(f"xv_mcd[0].location_: {xv_mcd[0].location_}")
      print(f"xv_mcd[0].covariance_: {xv_mcd[0].covariance_}")
      '''
      for x in P[slice_points]:
        xdotv = []
        for v in Vh:
          xdotv.append(np.dot(x,v))
        xv.append(xdotv)
        
      xv_med = np.median(xv, axis=0)
      num = xv - xv_med
      print(f"X: {P[slice_points]}")
      print(f"xv_med: {xv_med}")
      print(f"num: {num}")
      xv_mad = np.median(np.abs(xv - xv_med), axis=0)
      print(f"xv_mad: {xv_mad}")
      plot_slices(P, slice_points, Vh, f"slice_{i}")
      outl = []
      for num_i in num:
        outl.append(np.max(np.abs(np.divide(num_i, xv_mad))))
      if DEBUG:
        fig, ax = plt.subplots(tight_layout=True)
        ax.hist(outl, bins= 20)
        plt.savefig(f"tests/outl_hist_slice_{i}.png")
      #print(f"outl: {outl}")



      
    #print(trunk_points)
    #print(type(trunk_points))
    #U, S, Vh = np.linalg.svd(P[trunk_points])
    plot_pc(P[trunk_points], "segment_0_tree_0")
  print(f"cover.keys(): {cover.keys()}")
  return 0
"""

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
    print(f"trunk: {trunk}")

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

