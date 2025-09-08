import numpy as np

def func_grad_circle_centre(P,par,weight=None):
  '''
  ---------------------------------------------------------------------
  FUNC_GRAD_CIRCLE.M   Function and gradient calculation for 
                            least-squares circle fit.
 
  Version 1.0
  Latest update     20 Oct 2017
 
  Copyright (C) 2017 Pasi Raumonen
  ---------------------------------------------------------------------
 
  Input 
  P         Point cloud
  par       Circle parameters [x0 y0 r]'
  weight    Weights for the points. Weight the distances.
  
  Output
  dist      Signed distances of points to the circle:
                dist(i) = sqrt((xi-x0)^2 + (yi-y0)^2) - r, where 
                
  J         Jacobian matrix d dist(i)/d par(j).
  '''

  # Calculate distances
  Vx = P[:,0] - par[0]
  Vy = P[:,1] - par[1]
  rt = np.sqrt(Vx*Vx + Vy*Vy)
  if weight != None:
    dist = weight*(rt-par[2]) # Weighted distances to the circle
  else:
    dist = rt - par[2] # Distances to the circle
  
  # form the Jacobian matrixx

  m = len(P)
  J = np.zeros((m,2))
  J[:,0] = -Vx/rt
  J[:,1] = -Vy/rt

  # apply the weights
  if weight != None:
    J = np.concatenate((weight*J[:,0], weight*J[:,1]))

  return dist, J 
