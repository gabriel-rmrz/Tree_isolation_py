import numpy as np
from least_squares_fitting.form_rotation_matrices import form_rotation_matrices

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

def func_grad_cylinder(par,P,weight=None):
  '''
  ---------------------------------------------------------------------
  FUNC_GRAD_CYLINDER.M   Function and gradient calculation for 
                 least-squares cylinder fit.
 
  Version 2.2.0
  Latest update     5 Oct 2021
 
  Copyright (C) 2013-2021 Pasi Raumonen
  ---------------------------------------------------------------------
 
  Input 
  par       Cylinder parameters [x0 y0 alpha beta r]'
  P         Point cloud
  weight    (Optional) Weights for the points
  
  Output
  dist      Signed distances of points to the cylinder surface:
                dist(i) = sqrt(xh(i)^2 + yh(i)^2) - r, where 
                [xh yh zh]' = Ry(beta) * Rx(alpha) * ([x y z]' - [x0 y0 0]')
  J         Jacobian matrix d dist(i)/d par(j).
 
  Five cylinder parameters: 
  Location = axis point intersects xy-plane: x0 and y0 (z0 == 0)
  Rotation angles around x and y axis = alpha and beta
  Radius = r
 
  Transformed points:
  Pt = [xh yx zh] = Ry(beta) * Rx(alpha) * (P - [x0 y0 0])
 
  "Plane points":
  Qt = Pt * I2 = [xh yh];
 
  Distance:
  D(x0,y0,alpha,beta,r) = sqrt( dot(Qt,Qt) )-r = sqrt( Qt*Qt' )-r
 
  Least squares = minimize Sum( D^2 ) over x0, y0, alpha, beta and r
 
  rt = sqrt( dot(Qt,Qt) )
  N = Qt/rt
 
  Jacobian for D with respect to x0, y0, alpha, beta:
  dD/di = dot( N,dQt/di ) = dot( Qt/rt,dQt/di )
 
  R = Ry*Rx
  dQt/dx0 = R*[-1 0 0]'
  dQt/dy0 = R*[0 -1 0]'
  dQt/dalpha = (P-[x0 y0 0])*DRx';
  dQt/dalpha = (P-[x0 y0 0])*DRy';
 
  Changes from version 2.1.0 to 2.2.0, 5 Oct 20201:
  1) Minor changes in syntax
 
  Changes from version 2.0.0 to 2.1.0, 14 July 2020:
  1) Added optional input for weights of the points
  '''
  x0 = par[0]
  y0 = par[1]
  alpha = par[2]
  beta = par[3]
  r = par[4]

  # Determine the rotation matrices and their derivatives

  R, DR1, DR2 = form_rotation_matrices([alpha, beta])

  # Calculate the distances
  Pt = (P-[x0, y0, 0])@np.transpose(R)
  xt = Pt[:,0]
  yt = Pt[:,1]
  rt = np.sqrt(xt*xt +yt*yt)
  dist = rt - r # Distances to the cylinder surface
  if np.any(weight != None):
    dist = weight*dist
  
  # form the Jacobian matrix
  N = np.column_stack([xt/rt, yt/rt])
  m = len(P)
  J = np.zeros((m, 5))

  A1 = np.transpose(R@np.transpose([-1,0,0]))
  A = np.identity(2)
  A[0,0] = A1[0]
  A[1,1] = A1[1]
  J[:,0] = np.sum(N[:, 0:2]@A, axis=1)

  A2 = np.transpose(R@np.transpose([0,-1,0]))
  A[0,0] = A2[0]
  A[1,1] = A2[1]
  J[:,1] = np.sum(N[:, 0:2]@A, axis=1)

  A3 = (P- np.array([x0, y0, 0]))@np.transpose(DR1)
  J[:,2] = np.sum(N[:,0:2]*A3[:,0:2], axis=1)
  
  A4 = (P- np.array([x0, y0, 0]))@np.transpose(DR2)
  J[:,3] = np.sum(N[:,0:2]*A4[:,0:2], axis=1)

  J[:,4] = -1*np.ones(m)

  if np.any(weight != None):
    J = np.column_stack([weight*J[:,0], weight*J[:,1], weight*J[:,2], weight*J[:,3], weight*J[:,4]])

  return dist, J

