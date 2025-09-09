import numpy as np
from least_squares_fitting.func_grad_circle_centre import func_grad_circle_centre

def least_squares_circle_centre(P,Point0,Rad0):
  '''
  ---------------------------------------------------------------------
  LEAST_SQUARES_CIRCLE_CENTRE.M   Least-squares circle fitting such that
                                    radius is given (fits the centre)
 
  Version 1.0.0
  Latest update     6 Oct 2021
 
  Copyright (C) 2017-2021 Pasi Raumonen
  ---------------------------------------------------------------------
  Input    
  P         2d point cloud
  Point0    Initial estimate of centre (1 x 2)
  Rad0      The circle radius
  weight    Optional, weights for each point
  
  Output  
  cir     Structure array with the following fields
    Rad       Radius of the cylinder
    Point     Centre point (1 x 2)
    ArcCov    Arc point coverage (%), how much of the circle arc is covered 
                with points
    conv      If conv = 1, the algorithm has converged 
    rel       If rel = 1, the algorithm has reliable answer in terms of
                matrix inversion with a good enough condition number
  ---------------------------------------------------------------------
 
  Changes from version 1.0.0 to 1.1.0, 6 Oct 2021:  
  1) Streamlining code and some computations
  '''

  ## Initial estimates and other settings
  par = np.transpose(np.concatenate([np.atleast_1d(Point0), np.atleast_1d(Rad0)]))
  maxiter = 200 # maximum number of Gauss-Newton iteration
  iter_ = 0 # number of iteration so far
  conv = False # converge of Gauss-Newton algorithm
  rel = True # the results reliable (system matrix was not badly conditione)

  ## Gauss-Newton iterations
  while (iter_ < maxiter) and (~conv and rel):
    # Calculate the distancees and Jaconian
    dist, J = func_grad_circle_centre(P, par)

    # Calculate update step and gradient
    SS0 = np.linalg.norm(dist) # Squared sum of the distances
    # solve for the system of equations: par[i+1] - (J'J)^-1*J'd[par[i]]
    A = np.transpose(J)@J
    b = np.transpose(J)@dist
    p = np.linalg.solve(-A, b) # solve for the system of equations

    # Update
    par[:2] =par[:2] + p

    # Check if the updated parameters lower the squared sum value
    dist, J_ = func_grad_circle_centre(P, par)
    SS1 = np.linalg.norm(dist)
    if SS1 > SS0:
      # Update did not decreased the squared sum, use update with much
      # shorter update step
      par[:2] = par[:2] - 0.95*p
      dist, J_ = func_grad_circle_centre(P,par)
      SS1 = np.linalg.norm(dist)
    
    # Check reliability
    rcond_est = 1.0/ np.linalg.cond(A, p=1)
    if rcond_est < 10000*np.finfo(float).eps:
      rel = False

    # Check convergence
    if np.abs(SS0-SS1) < 1e-5:
      conv = True

    iter_+= 1
  
  ## Output
  Point = par[:2]
  if conv and rel:
    # Calculate ArcCov, how much of the circle arc is covered with points
    U = P[:,0] - par[0]
    V = P[:,1] - par[1]
    ang = np.arctan(V,U) + np.pi
    I = np.zeros(100, dtype=bool)
    ang = np.ceil(ang/(2*np.pi)*100).astype(int)
    I[ang] = True
    ArcCov = np.count_nonzero(I)/100
    # mean absoluto distance to the circle
    d = np.sqrt(U*U + V*V) - Rad0
    mad = np.mean(np.abs(d))
  else:
    mad = 0
    Arco = 0
  
  cir = {}
  cir['radius'] = Rad0
  cir['point'] =np.transpose(Point)
  cir['mad'] = np.atleast_1d(mad)
  cir['ArcCov'] = ArcCov
  cir['conv'] = conv
  cir['rel'] = rel

  return cir
