import numpy as np
from tools.distances_to_line import distances_to_line
from tools.orthonormal_vectors import orthonormal_vectors
def surface_coverage_filtering(P,c,lh,numS):
  '''
   ---------------------------------------------------------------------
   SURFACE_COVERAGE_FILTERING.M    Filters a point cloud based on the 
                                     assumption that it samples a cylinder
  
   Version 1.1.0
   Latest update     6 Oct 2021
  
   Copyright (C) 2017-2021 Pasi Raumonen
   ---------------------------------------------------------------------
  
   Filter a 3d-point cloud based on given cylinder (axis and radius) by 
   dividing the point cloud into "ns" equal-angle sectors and "lh"-height
   layers along the axis. For each sector-layer intersection (a region in 
   the cylinder surface) keep only the points closest to the axis. 
  
   Inputs:
   P             Point cloud, (n_points x 3)-matrix
   c             Cylinder, stucture array with fields "axis", "start",
                   "length"
   lh            Height of the layers
   ns            Number of sectors
  
   Outputs:              
   Pass          Logical vector indicating which points pass the filtering
   c             Cylinder, stucture array with additional fields "radius",
                   "SurfCov", "mad", "conv", "rel", estimated from the
                   filtering
   ---------------------------------------------------------------------
  
   Changes from version 1.0.0 to 1.1.0, 6 Oct 2021:  
   1) Small changes to make the code little faster
   2) Change the radius estimation to make it much faster
  
  
   Compute the distances, heights and angles of the points
  '''
  d, V, h, _ = distances_to_line(P,c['axis'],c['start'])
  h = h- np.min(h)
  U, W = orthonormal_vectors(c['axis'])
  V = V @ np.column_stack((U, W))
  ang = np.arctan2(V[:,0], V[:,1]) + np.pi

  # Sort based on lexicographic order of (sector, layer)
  numL = int(np.ceil(c['length']/lh))
  Layer = np.ceil(h/c['length']*numL).astype(int)
  Layer[Layer==0] = 1
  Layer[Layer> numL] = numL
  Sector = np.ceil(ang/(2*np.pi)*numS).astype(int)
  Sector[Sector == 0] = 1
  LexOrd = np.dot(np.column_stack((Layer, Sector-1)),np.transpose([1, numL]))
  SortOrd = np.argsort(LexOrd)
  LexOrd = LexOrd[SortOrd]
  ds = d[SortOrd]

  # Estimate the distances for each sector-layer intersection
  Dis = np.zeros((numL, numS))
  numP = len(P)
  p = 0
  while p < numP:
    t = 1
    while (p+t < numP) and (LexOrd[p] == LexOrd[p+t]):
      t+=1
    D = np.min(ds[p:p+t-1+1])
    flat_Dis = Dis.ravel(order='F')
    flat_Dis[LexOrd[p]-1] = np.min((1.05*D, D+0.02))
    Dis = flat_Dis.reshape(Dis.shape, order='F')
    p+=t

  # Compute the number of sectors based on the estimated radius
  R = np.median(Dis[Dis>0.])
  a = np.max((0.02, 0.02*R))
  numS = int(np.ceil(2*np.pi*R/a))
  numS = int(np.min([36,np.max([numS,8])]))
  numL = int(np.ceil(c['length']/a))
  numL = int(np.max((numL, 3)))


  # Sort based on lexicographic order of (sector, layer)
  numL = int(np.ceil(c['length']/lh))
  Layer = np.ceil(h/c['length']*numL).astype(int)
  Layer[Layer==0] = 1
  Layer[Layer> numL] = numL
  Sector = np.ceil(ang/(2*np.pi)*numS).astype(int)
  Sector[Sector == 0] = 1
  LexOrd = np.dot(np.column_stack((Layer, Sector-1)),np.transpose([1, numL]))
  SortOrd = np.argsort(LexOrd)
  LexOrd = LexOrd[SortOrd]
  ds = d[SortOrd]
  # Filtering for sector-layer intersection
  Dis = np.zeros((numL, numS), dtype=int)
  Pass = np.zeros(numP, dtype=bool)
  p = 0 # index of point under processing
  k = 0 # number of nonempty cells
  r = np.max((0.01, 0.05*R)) # cell diameter from the closest point
  while p < numP:
    t = 1
    while (p+t < numP) and (LexOrd[p] == LexOrd[p+t]):
      t+=1
    ind = np.array(list(range(p,p+t-1+1)))
    D = d[ind]
    Dmin = np.min(D)
    I = np.array((D <= Dmin+r), dtype=bool)
    Pass[ind[I]] = True
    flat_Dis = Dis.ravel(order='F')
    flat_Dis[LexOrd[p]-1] = np.min((1.05*Dmin, Dmin+0.02))
    Dis = flat_Dis.reshape(Dis.shape, order='F')
    p+=t
    k+=1
  d = d[Pass]

  # Sort the "Pass"-vector back to original point cloud order
  n = len(SortOrd)
  InvSortOrd = np.zeros(n, dtype=int)
  InvSortOrd[SortOrd] = np.array(list(range(n)), dtype=int)
  Pass = Pass[InvSortOrd]

  # Compute raius, SurfCov and mad
  R = np.median(Dis[Dis>0])
  mad = np.sum(np.abs(d-R))/len(d)

  c['radius'] = R
  c['SurfCov'] = k/(numL*numS)
  c['mad'] = mad
  c['conv'] = 1
  c['rel'] = 1
  return Pass, c




