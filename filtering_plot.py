import time
import numpy as np

from tools.display_time import display_time
from tools.restrict_size import restrict_size
from tools.cubical_down_sampling import cubical_down_sampling
from tools.remove_lowest_points import remove_lowest_points
from tools.cubical_density_filtering import cubical_density_filtering
from tools.compute_height import compute_height
from tools.remove_bottom import remove_bottom
from tools.remove_small_separate_clusters import remove_small_separate_clusters

def filtering_plot(P, inputs):
  '''
  FILTERING_PLOT.M      Filters plot-level point clouds with multiple steps.
 
  Version 2.0.0
  Latest update       18 Dec 2019
 
  Copyright (C) 2015-2019 Pasi Raumonen
  ---------------------------------------------------------------------
 
  Filtering based on 5 different steps:
  1.    (Optional) Restricts the size of the plot by removing far away points
  2.    (Optional) Downsampling based on cubical partitioning
  3.    (Optional) Removes likely ground points
  4.    (Optional) Filters based on cubical point density
  5.    (Optional) Filters using cover sets, based on point density and clustering
 
  Details for each step:
  1.    Input parameter "Bound" defines the outer boundary. Either it is a
        number giving the radius of the boundary or it is set of planar
        points (x,y-coordinates) describing (convex?) polygon, which is the
        boundary. All points outside the boundary are removed.
  2.    Partitions the point cloud into cubes with side length of
        "CubeSizeSampling" and selects one point from each cube.
  3.    Estime the ground level by dividing the points into rectangles and
        remove the bottom layer.
  4.    Partition the points into cubes with edge length "CubeSizeDensity" to
        remove low point density regions. Remove all the points in the cube
        if there is less than "PointDensity" number of points.
  5.    Use cover sets to remove small separate clusters. Remove all connected
        components with less than "FiltClusterSize" cover sets.
 
  Inputs:
  P                 Plot-level point cloud
  inputs            Structure array. Needs the following fields:
    1. Point cloud restriction:
      Bound             Plot radius or boundary, outside points will be removed
    2. Downsampling:
      CubeSizeSampling  Side length of the cubes used for the downsampling
    3. Ground filtering:
      BottomSquare      Side length of rectangles used to define ground level
      BottomHeight      Height above the ground, below every point is removed
    4. Low point density filtering:
      CubeSizeDensity   Side length of the cubes used for low density filtering
      PointDensity      Threshold for the low point density, the minimum
                            number of points in the cubes.
    5. Cluster filtering:
      filt.PatchDiam    Patch diameter for cover generation with cover_sets_plot
      filt.BallRad      Ball radius for cover generation with cover_sets_plot
      filt.nmin         nmin cover generation with cover_sets_plot
      filt.NCubes       Number of cubes for the partition used in
                            cover generation with cover_sets_plot
      FiltClusterSize   Threshold for numnber of cover sets in an acceptable cluster
 
  Outputs:
  Points    Logical vector indicating the points passing the filtering
  Hei       Height of every point passing the filterings
  Ground    Ground level nodes, (x,y,z)-coordinates
  Tri       Triangulation of the ground level, list of nodes
  pass
  '''

  
  print('  -----------------------')
  print('   Filtering point cloud...')
  print('  -----------------------')



  i_time = time.time() # Initial running time
  #print(P)
  numP0 = len(P[:,0]) # Number of points before any filtering
  bound = np.asarray(inputs["bound"])
  Points = np.ones(numP0, dtype='bool')
  # Restrict the plot size
  if bound[0][0] != 0:
    ti_time = time.time()
    Points = restrict_size(P, inputs)
    P = P[Points]

    display_time(i_time, ti_time,'\t Size restriction:', 1)
  
  if inputs['SamplingCubeSize'] > 0:
    ## Cubical downsampling
    ti_time = time.time()
    P, Points = cubical_down_sampling(P, Points, inputs)

    display_time(i_time, ti_time,'\t Downsampling:', 1)
  ## Remove the lowest points
  # (possible "ghost" points clearly below the ground level) TODO: Check carefully this
  ti_time = time.time()
  P, Points = remove_lowest_points(P, Points)
  display_time(i_time, ti_time,'\t Removing lowest points:', 1)
  
  if inputs['CubeSizeDensity'] > 0:
    ti_time = time.time()
    P, Points = cubical_density_filtering(P, Points, inputs)
    display_time(i_time, ti_time,'\t Cube filtering: ', 1)

  ti_time = time.time()
  Hei, Ground, Tri = compute_height(P, inputs)
  display_time(i_time, ti_time,'\t Height computation: ', 1)

  if inputs['BottomHeight'] > 0:
    ti_time = time.time()
    P, Points, Hei = remove_bottom(P,Points,Hei,inputs)
    display_time(i_time, ti_time,'\t Removing bottom: ', 1)

  if inputs['filt']['PatchDiam'] > 0:
    ti_time = time.time()
    Points = remove_small_separate_clusters(P,Points,inputs)
    display_time(i_time, ti_time,'\t Removing small clusters: ', 1)

  Hei = Hei[Points]
  numP = len(Points[Points])
  print('---------')
  print('Filtering summary:')
  print(f"\t Points filtered: {numP0 - numP}")
  print(f"\t Points left: {numP}")
  print(f"\t Total filtering: {round((numP0-numP)/numP0*1000)/10}")
  print('---------')

  f_time = time.time()

  return Points, Hei, Ground, Tri

