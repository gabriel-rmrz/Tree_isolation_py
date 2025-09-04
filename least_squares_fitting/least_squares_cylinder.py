import numpy as np
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

def least_squares_cylinder(P,cyl0,weight,Q):
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
  return 0, 1
