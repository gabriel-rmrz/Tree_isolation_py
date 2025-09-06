import numpy as np
from least_squares_fitting.rotate_to_z_axis import rotate_to_z_axis
from least_squares_fitting.func_grad_cylinder import func_grad_cylinder
from least_squares_fitting.form_rotation_matrices import form_rotation_matrices
from tools.surface_coverage import surface_coverage
'''
 This file is part of TREEQSM.
 
 TREEQSM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TREEQSM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.
'''

def least_squares_cylinder(P,cyl0,weight=None,Q=None, I=None):
  '''
  ---------------------------------------------------------------------
  LEAST_SQUARES_CYLINDER.M   Least-squares cylinder using Gauss-Newton.
 
  Version 2.0.0
  Latest update     5 Oct 2021
 
  Copyright (C) 2013-2021 Pasi Raumonen
  ---------------------------------------------------------------------
  Input    
  P         Point cloud
  cyl0      Initial estimates of the cylinder parameters
  weight    (Optional) Weights of the points for fitting
  Q         (Optional) Subset of "P" where the cylinder is intended
  
  Output  
  cyl       Structure array containing the following fields:
    radius      Radius of the cylinder
    length      Length of the cylinder
    start       Point on the axis at the bottom of the cylinder (1 x 3)
    axis        Axis direction of the cylinder (1 x 3) 
    mad         Mean absolute distance between points and cylinder surface
    SurfCov     Relative cover of the cylinder's surface by the points 
    dist        Radial distances from the points to the cylinder (m x 1) 
    conv        If conv = 1, the algorithm has converged 
    rel         If rel = 1, the algorithm has reliable answer in terms of
                    matrix inversion with a good enough condition number
  ---------------------------------------------------------------------
 
  Changes from version 1.3.0 to 2.0.0, 5 Oct 2021:  
  1) Included the Gauss-Newton iterations into this function (removed the 
       call to nlssolver function)
  2) Changed how the updata step is solved from the Jacobian
  3) Simplified some expressions and added comments
  4) mad is computed only from the points along the cylinder length in the
      case of the optional input "Q" is given.  
  5) Changed the surface coverage estimation by filtering out points whose 
      distance to the axis is less than 80% of the radius 
 
  Changes from version 1.2.0 to 1.3.0, 14 July 2020:  
  1) Changed the input parameters of the cylinder to the struct format.
  2) Added optional input for weights
  3) Added optional input "Q", a subset of "P", the cylinder is intended
     to be fitted in this subset but it is fitted to "P" to get better
     estimate of the axis direction and radius
 
  Changes from version 1.1.0 to 1.2.0, 14 Jan 2020:  
  1) Changed the outputs and optionally the inputs to the struct format.
  2) Added new output, "mad", which is the mean absolute distance of the
     points from the surface of the cylinder.
  3) Added new output, "SurfCov", that measures how well the surface of the
     cylinder is covered by the points.
  4) Added new output, "SurfCovDis", which is a matrix of mean point distances 
     from layer/sector-intersections to the axis.
  5) Added new output, "SurfCovVol", which is an estimate of the cylinder's 
     volume based on the radii in "SurfCovDis" and "cylindrical sectors".
  6) Added new optional input "res" which gives the point resolution level
     for computing SurfCov: the width and length of sectors/layers.
 
  Changes from version 1.0.0 to 1.1.0, 3 Oct 2019:  
  1) Bug fix: --> "Point = Rot0'*([par(1) par(2) 0]')..."
  '''
  ## Initialize data and values
  print('Here')
  res = 0.03 # "Resolution level" for computing surface coverage
  maxiter = 50 # maximum number of Gauss-Newton iterations
  iter_ = 0
  conv = False # Did the iteration converge
  rel = True # Are the results reliable (condition number was not very bad)
  print(weight)
  if np.any(weight == None):
    NoWeights = True
  else:
    NoWeights = False
  
  # Transform the data to close to standard position via a translation
  # followed by a rotation
  Rot0, D_, a_ = rotate_to_z_axis(cyl0['axis'])
  Pt = (P-cyl0['start'])@np.transpose(Rot0)

  # Initial estimates
  par = np.transpose([0, 0, 0, 0, cyl0['radius']])


  ## Gauss-Newton algorithm
  # find estimate of rotation-translation-radius parameters that transform
  # the data so that the best-fit cylinder is one in standard position

  while iter_ < maxiter and ~conv and rel:
    ## Calculate the distances and Jacobian
    if NoWeights:
      d0, J = func_grad_cylinder(par, Pt)
    else:
      d0, J = func_grad_cylinder(par, Pt, weight)

    ## Calculate update step
    SS0 = np.linalg.norm(d0) # Squared sum of the distance
    # solve for the system of equations:
    # par[i+1] = par[i] - (J'J)^(-1)*J'd0[par[i]]
    A = np.transpose(J)@J
    b = np.transpose(J)@d0
    p = np.linalg.solve(-A, b)
    par = par+p # update the parameters
    ## Check reliability
    eps_float = np.finfo(np.float32).eps
    rcond_est = 1.0/ np.linalg.cond(-A, p=1)

    if rcond_est < 10000*eps_float:
      rel = False

    ## Check convergence:
    # The distances with the new parameter values:
    if NoWeights:
      dist, _ = func_grad_cylinder(par, Pt)
    else:
      dist, _ = func_grad_cylinder(par, Pt, weight)

    SS1 = np.linalg.norm(dist) # Squared sum of the distances
    if np.abs(SS0-SS1) < 1e-4:
      conv = True

    iter_ +=1
  ## Compute the cylinder parameters and other outputs
  cyl = {}
  cyl['radius'] = float(par[4]) # radius

  # Inverse transformation fo find axis and point on axis
  # correspoding to orignal data
  Rot, DR1_, DR2_ = form_rotation_matrices(par[2:4])
  Axis = np.transpose(Rot0)@np.transpose(Rot)@np.transpose(np.array([0, 0, 1])) # axis direction
  Point = np.transpose(Rot0)@np.transpose(np.array([par[0], par[1], 0])) + np.transpose(cyl0['start']) # axis point

  if np.any(Q != None):
    if len(Q) > 5:
      P = Q
  H = P@Axis # Heights along the axis
  hmin = np.min(H)
  cyl['length'] = float(np.abs(np.max(H) - hmin))
  hpoint = np.transpose(Axis)@Point
  Point = Point - (hpoint- hmin)*Axis # axis point at the cylinder's bottom
  cyl['start'] = np.transpose(Point).astype(float)
  cyl['axis'] =  np.transpose(Axis).astype(float)
  print(f"cyl['axis']: {cyl['axis']}")
  print(f"cyl['start']: {cyl['start']}")
  if np.any(weight != None) and I != None:
    I = (weight == np.max(weight))
    cyl['mad'] = float(np.average(np.abs(dist(I)))) # mean ab
  else:
    cyl['mad'] = float(np.average(np.abs(dist))) # mean ab
  cyl['conv'] = conv
  cyl['rel'] = rel

  # Compute SurfCov, minimum 3*8 grid
  if ~np.any(np.isnan(Axis)) and ~np.any(np.isnan(Point)) and rel and conv:
    numL = int(np.max((3, np.ceil(cyl['length']/res))))
    numS = np.ceil(2*np.pi*cyl['radius']/res)
    numS = int(np.min((36., np.max((numS,8)))))
    SurfCov, Dis_, CylVol_, dis_ = surface_coverage(P,np.transpose(Axis),np.transpose(Point), numL, numS, 0.8*cyl['radius'])

    cyl['SurfCov'] = float(SurfCov)
  else:
    cyl['SurfCov'] = float(0)

  return cyl
  
  exit() 

  return 0, 1
