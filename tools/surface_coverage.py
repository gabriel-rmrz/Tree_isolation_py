import numpy as np
from tools.distances_to_line import distances_to_line
from tools.orthonormal_vectors import orthonormal_vectors
from tools.rotation_matrix import rotation_matrix

def surface_coverage(P,Axis,Point,numL,numS,Dmin=None,Dmax=None):

  '''
  ---------------------------------------------------------------------
  SURFACE_COVERAGE.M   Computes point surface coverage measure
 
  Version       1.1.0
  Last update   7 Oct 2021
 
  Copyright (C) 2017-2021 Pasi Raumonen
  ---------------------------------------------------------------------
  Inputs:    
  Axis      Axis direction (1 x 3) 
  Point     Starting point of the cylinder (1 x 3)
  nl        Number of layers in the axis direction used for to partition
                the cylinder surface into layer/sectors
  ns        Number of sectors used to partition the cylinder surface into 
                layer/sectors
  Dmin      (Optional) Minimum point distance from the axis to be included
                into SurfCov calculations
  Dmax      (Optional) Maximum point distance from the axis to be included
                into SurfCov calculations
  
  Output:  
  SurfCov   Number between 0 and 1 descring how big portion of the cylinder
                surface is covered with points
  Dis       (Optional) Mean distances of the distances of the layer/sectors
  CylVol    (Optional) Volume of the cylinder estimated by the mean
                distances of the layer/sectors as cylindrical sections
  dis       (Optional) Same as "Dis" but empty cells are interpolated
  ---------------------------------------------------------------------
  Computes surface coverage (number between 0 and 1) of points on cylinder 
  surface defined by "Axis" and "Point".
 
  Changes from version 1.0.0 to 1.1.0, 7 Oct 2021:
  1) Added two possible inputs, minimum and maximum distance, 
     Dmin and Dmax, which can be used to filter out points for the surface
     coverage calculations
  2) Computes the SurfCov estimate with four baseline directions used in
     the sector determination and selects the largest value
  3) Smalle changes to speed up computations
  '''
  ## Compute the distances and heights of the points
  d, V, h, D_ = distances_to_line(P, Axis, Point)
  h = h-np.min(h)
  Len = np.max(h)

  ## (Optional)
  if Dmin != None:
    Keep = (d > Dmin)
    if Dmax != None:
      Keep = Keep & (d < Dmax)
    V = V[Keep, :]
    h = h[Keep]
  
  ## Compute SurfCov
  # from 4 different baseline directions to determine the angles and select
  # the maximum value
  V0 = V
  U, W = orthonormal_vectors(Axis) # First planar axes
  R = rotation_matrix(Axis, 2*np.pi/(numS*4)) # Rotation matrix to rotate the axes
  SurfCov = np.zeros((1,4))
  for i in range(4):
    ## Rotate the axes
    if i > 0:
      U = R@U
      W = R@W

    ## compute the angles (sectors) of the points
    V = V0@np.column_stack((U,W))
    ang= np.arctan2(V[:,1], V[:,0]) + np.pi

    ## Compute lexicographic order (sector, layer) of every point
    Layer = np.ceil(h/Len*numL).astype(int)
    Layer[Layer <= 0]=1
    Layer[Layer > numL] = numL
    Sector = np.ceil(ang/(2*np.pi)*numS).astype(int)
    Sector[Sector <= 0] = 1
    LexOrd = np.column_stack([Layer, Sector-1])@np.transpose([1, numL]).astype(int)

    ## Compute SurfCov
    Cov = np.zeros((numL,numS))
    flat_Cov = Cov.ravel(order='F')
    print(f"flat_Cov: {flat_Cov}")
    print(f"Cov: {Cov}")
    print(f"numL: {numL}")
    print(f"numS: {numS}")
    flat_Cov[LexOrd -1] = 1
    print(f"flat_Cov: {flat_Cov}")
    SurfCov[0,i] = np.count_nonzero(flat_Cov)/(1.*numL*numS)
  SurfCov =np.max(SurfCov)
  print(SurfCov)

    ## Compute volume estimate


  exit()






