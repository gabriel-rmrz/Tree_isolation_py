import numpy as np
import cmath, math

from tools.mat_vec_subtraction import mat_vec_subtraction
from tools.distances_between_lines import distances_between_lines
def parent_cylinder(SPar, SChi, CiS, cylinder, cyl, si):
  # Finds the parent cylinder from the possible parent segment.
  # Does this by checking it the axis of the cylinder, if continued, will
  # cross the nearby cylinders in the parent segment.
  # Adjust the cylinder so that it starts from the surface of its parent.

  rad = cyl['radius']
  len_ = cyl['length']
  sta = cyl['start']
  axe = cyl['axis']

  # PC:   Parent Cylinder
  numC = rad.size
  added = False
  print(f"SPar[si]: {SPar[si]}")
  if np.any(SPar[si] > 0): # parent segment exist, find parent cylinder
    s = SPar[si]
    PC = np.concatenate([CiS[s_] for s_ in s if s_ in CiS]).astype(int).tolist() # the cylinders in the parent segment
    # select the closest cylinders for closer examination
    print(f"HERE: PC: {PC}")
    if len(PC) > 1:
      print(f"HERE: len(PC): {len(PC)}")
      print(f"HERE: cylinder['start'].shape: {cylinder['start'].shape}")
      print(f"HERE: sta.shape: {sta.shape}")
      D = mat_vec_subtraction(-cylinder['start'][PC,:], -np.atleast_2d(sta)[0,:])
      d = np.sum(D*D,1)
      I = np.argsort(d).astype(int).tolist()
      if len(PC) > 3:
        I = I[:4]
      print(f"HERE: I: {I}")
      print(f"HERE: PC: {PC}")
      pc = PC[I]
      ParentFound = False
    elif len(PC) ==1:
      ParentFound = True
    else:
      PC = np.zeros((0,1))
      ParentFound = True

    ## Check possible crossing points
    if not ParentFound:
      pc0 = pc
      n = len(pc)
      # Calculate the possible crossing points of the cylinder axis, when
      # extended, on the surfaces of the parent candidate cylinders
      x = np.zeros((n,2)) # how much the starting point has to move to cross
      h = np.zeros((n,2)) # the crossing point height in the parent
      Axe = cylinder['axis'][pc,:]
      Sta = cylinder['start'][pc,:]
      for j in range(n):
        # Crossing point solver from a quadratic equation
        A = axe[0,:]-(axe[0,:]@np.transpose(Axe[j,:]))@Axe[j,:]
        B = sta[0,:]-Sta[j,:]-(sta[0,:]*np.transpose(Axe[j,:]))@Axe[j,:]+(Sta[j,:]@np.transpose(Axe[j,:]))@Axe[j,:]

        e = A@np.transpose(A)
        f = 2*A@np.transpose(B)
        g = B@np.transpose(B) - cylinder['radius'][pc[j]]**2
        di = cmath.sqrt(f**2 - 4*e*g) # the discriminant
        s1 = (-f + di)/(2*e)
        # how much the starting point must be moved to cross
        s2 = (-f - di)/(2*e)

        has_imag = isinstance(s1, complex) and not math.isclose(s1.imag, 0.0, rel_tol=0, abs_tol=1e-12)
        if not has_imag: #cylinders can cross
          # the height of the crossing points
          x[j,:] = np.array([s1,s2])
          h[j,0] = sta[0,:]@np.transpose(Axe[j,:])+x[j,0]*axe[0,:]@np.transpose(Axe[j,:])-Sta[j,:]@np.transpose(Axe[j,:])
          h[j,1] = sta[0,:]@np.transpose(Axe[j,:])+x[j,1]*axe[0,:]@np.transpose(Axe[j,:])-Sta[j,:]@np.transpose(Axe[j,:])
      
      ## Extend to crossing point in the (extended) parent
      I = (x[:,0] != 0) # Select only candidates with crossing points
      pc = pc0[I]
      x = x[I,:]
      h = h[I,:]
      j = 0
      n = np.count_nonzero(I)
      X = np.zeros((n,3))
      Len = cylinder['length'][pc]
      while j < n and ~ParentFound:
        if x[j,0] > 0 and x[j,1] < 0:
          # sp inside the parent and crosses its surface
          if h[j,0] >= 0 and h[j,0] <= Len[j] and len_[0]-x[j,0] > 0:
            PC = pc[j]
            sta[0,:] = sta[0,:]+x[j,0]*axe[0,:]
            len_[0] = len_[0] - x[j,0]
            ParentFound = True
          elif len_[0] - x[j,0] > 0:
            if h[j,0] < 0:
              X[j,:] = np.array([x[j,0], np.abs(h[j,0]), 0])
            else:
              X[j,:] = np.array([x[j,0], h[j,0]-Len[j], 0])
          else:
            X[j,:] = np.array([x[j,0], h[j,0], 1])
        elif x[j,0] < 0 and x[j,1] > 0 and len_[0]-x[j,1] > 0:
          # sp inside the parent and crosses its surface
          if h[j,1] >= 0 and h[j,1] <= Len[j] and len_[0]-x[j,1] > 0:
            PC = pc[j]
            sta[0,:] = sta[0,:]+x[j,1]*axe[0,:]
            len_[0] = len_[0] - x[j,1]
            ParentFound = True
          elif len_[0] - x[j,1] > 0:
            if h[j,1] < 0:
              X[j,:] = np.array([x[j,1], np.abs(h[j,1]), 0])
            else:
              X[j,:] = np.array([x[j,1], h[j,1]-Len[j], 0])
          else:
            X[j,:] = np.array([x[j,1], h[j,1], 1])

        elif x[j,0] < 0 and x[j,1] < 0 and x[j,1] < x[j,0] and len_[0]-x[j,1] > 0:
          # sp outside the parent and crosses its surface when extended 
          # backwards
          if h[j,0] >= 0 and h[j,0] <= Len[j] and len_[0]-x[j,1] > 0:
            PC = pc[j]
            sta[0,:] = sta[0,:]+x[j,0]*axe[0,:]
            len_[0] = len_[0] - x[j,0]
            ParentFound = True
          elif len_[0] - x[j,0] > 0:
            if h[j,0] < 0:
              X[j,:] = np.array([x[j,0], np.abs(h[j,0]), 0])
            else:
              X[j,:] = np.array([x[j,0], h[j,0]-Len[j], 0])
          else:
            X[j,:] = np.array([x[j,0], h[j,0], 1])

        elif x[j,0] < 0 and x[j,1] < 0 and x[j,1] > x[j,0] and len_[0]-x[j,1] > 0:
          # sp outside the parent and crosses its surface when extended
          # backwards
          if h[j,1] >= 0 and h[j,1] <= Len[j] and len_[0]-x[j,1] > 0:
            PC = pc[j]
            sta[0,:] = sta[0,:]+x[j,1]*axe[0,:]
            len_[0] = len_[0] - x[j,1]
            ParentFound = True
          elif len_[0] - x[j,1] > 0:
            if h[j,1] < 0:
              X[j,:] = np.array([x[j,1], np.abs(h[j,1]), 0])
            else:
              X[j,:] = np.array([x[j,1], h[j,1]-Len[j], 0])
          else:
            X[j,:] = np.array([x[j,1], h[j,1], 1])
            
        elif x[j,0] > 0 and x[j,1] > 0 and x[j,1] < x[j,0] and len_[0]-x[j,1] > 0:
          # sp outside the parent but crosses its surface when extended forward
          if h[j,0] >= 0 and h[j,0] <= Len[j] and len_[0]-x[j,1] > 0:
            PC = pc[j]
            sta[0,:] = sta[0,:]+x[j,0]*axe[0,:]
            len_[0] = len_[0] - x[j,0]
            ParentFound = True
          elif len_[0] - x[j,0] > 0:
            if h[j,0] < 0:
              X[j,:] = np.array([x[j,0], np.abs(h[j,0]), 0])
            else:
              X[j,:] = np.array([x[j,0], h[j,0]-Len[j], 0])
          else:
            X[j,:] = np.array([x[j,0], h[j,0], 1])

        elif x[j,0] > 0 and x[j,1] > 0 and x[j,1] > x[j,0] and len_[0]-x[j,1] > 0:
          # sp outside the parent and crosses it surface when extended forward
          if h[j,1] >= 0 and h[j,1] <= Len[j] and len_[0]-x[j,1] > 0:
            PC = pc[j]
            sta[0,:] = sta[0,:]+x[j,1]*axe[0,:]
            len_[0] = len_[0] - x[j,1]
            ParentFound = True
          elif len_[0] - x[j,1] > 0:
            if h[j,1] < 0:
              X[j,:] = np.array([x[j,1], np.abs(h[j,1]), 0])
            else:
              X[j,:] = np.array([x[j,1], h[j,1]-Len[j], 0])
          else:
            X[j,:] = np.array([x[j,1], h[j,1], 1])
      j+=1

    if not ParentFound and n > 0:
      I = np.argmin(X[:,1])
      H = X[I,0]
      X = X[I,:]
      if X[2] == 0 and H < 0.1*Len[I]:
        PC = pc[I]
        sta[0,:] = sta[0,:]+X[0]*axe[0,:]
        len_[0] = len_[0]-X[0]
        ParentFound = True
      else:
        PC = pc[I]
        if numC > 1 and X[0] <= rad[0] and np.abs(X[2]) <= 1.25 *cylinder['length'][PC]:
          # removes the first cylinder and adjust the second
          S = sta[0,:] + X[0]*axe[1,:]
          V = sta[1,:] +len_[1]*axe[1,:]-S
          len_[1] = np.linalg.norm(V)
          len_ = len_[1:numC]
          axe[1,:] = V/np.linalg.norm(V)
          axe = axe[1:numC,:]
          sta[1,:] = S
          sta = sta[1:numC,:]
          rad = rad[1:numC]
          cyl['mad'] = cyl['mad'][1:numC]
          cyl['SurfCov'] = cyl['SurfCov'][1:numC]
          numC = numC-1
          ParentFound = True
        elif numC > 1:
          # Remove first cylinder
          sta = sta[1:numC,:]
          len_ = len_[1:numC]
          axe = axe[1:numC,:]
          rad = rad[1:numC]
          cyl['mad'] = cyl['mad'][1:numC]
          cyl['SurfCov'] = cyl['SurfCov'][1:numC]
          numC = numC-1
        elif len(SChi[si]) == 0 :
          # Remove the cylinder
          numC = 0
          PC = np.zeros((0,1))
          ParentFound = True
          rad = np.zeros((0,1))
        elif X[0] <= rad[0] and np.abs(X[0]) <= 1.5*cylinder['length'][PC]:
          # Adjust the cylinder
          sta[0,:] = sta[0,:]+X[0]*axe[0,:]
          len_[0] = np.abs(X[0])
          ParentFound = True
    if not ParentFound:
      # The parent is the cylinder in the parent segment whose axis
      # line is the closest to the axis line of the first cylinder
      # Or the parent cylinder is the one whose base, when connected
      # to the first cylinder is the most parallel
      # Add new cylinder
      pc = pc0
      Dist, _, DistOnLines = distances_between_lines(sta[0,:], axe[0,:], cylinder['start'][pc,:], cylinder['axis'][pc,:])

      I = (DistOnLines >=0)
      J = (DistOnLines <= cylinder['length'][pc])
      I = I & J
      if not np.any(I):
        I = DistOnLines >= -0.2*cylinder['length'][pc]
        J = DistOnLines <= 1.2*cylinder['length'][pc]
        I = I&J
      if np.any(I):
        pc = pc[I]
        Dist = Dist[I]
        DistOnLines = DistOnLines[I]
        I = np.argmin(Dist)
        DistOnLines = DistOnLines[I]
        PC = pc[I]
        Q = cylinder['start'][PC,:] + DistOnLines*cylinder['axis'][PC,:]
        V = sta[0,:] - Q
        L = np.linalg.norm(V) 
        V = V/L
        a = np.arccos(V@np.transpose(cylinder['axis'][PC,:]))
        h = np.sin(a) * L
        S = Q + cylinder['radius'][PC]/h*L*V
        L = (h-cylinder['radius'][PC])/h*L
        if L > 0.01 and L/len_[0] > 0.2:
          numC = numC + 1
          sta = np.concatenate((S, sta))
          rad = np.concatenate(([rad[0]], rad))
          axe = np.concatenate((V, axe))
          len_ = np.concatenate(([L], len_))
          cyl['mad'] = np.concatenate(([cyl['mad'][0]], cyl['mad']))
          cyl['SurvCov'] = np.concatenate(([cyl['SurvCov'][0]], cyl['SurvCov']))
          cyl['rel'] = np.concatenate(([cyl['rel'][0]], cyl['rel']))
          cyl['conv'] = np.concatenate(([cyl['conv'][0]], cyl['conv']))
          added = True
      else:
        V = - mat_vec_subtraction(cylinder['start'][pc,:], sta[0,:])
        L0 = np.sqrt(np.sum(V*V, 1))
        V = np.column_stack((V[:,0]/L0, V[:,1]/L0, V[:,2]/L0))
        A = V@np.transpose(axe[0,:])
        I = np.argmax(A)
        A = A[I]
        L1 = L0[I]
        PC = pc[I]
        V = V[I,:]
        a = np.arccos(V@np.transpose(cylinder['axis'][PC,:]))
        h = np.sin(a)*L1
        S = cylinder['start'][PC,:]*cylinder['radius'][PC]/h*L1*V
        L = (h-cylinder[PC])/h*L1
        if L > 0.01 and L/len_[0] > 0.2:
          numC = numC+1
          sta = np.concatenate((S, sta))
          rad = np.concatenate(([rad[0]], rad))
          axe = np.concatenate((V, axe))
          len_ = np.concatenate(([L], len_))
          cyl['mad'] = np.concatenate(([cyl['mad'][0]], cyl['mad']))
          cyl['SurvCov'] = np.concatenate(([cyl['SurvCov'][0]], cyl['SurvCov']))
          cyl['rel'] = np.concatenate(([cyl['rel'][0]], cyl['rel']))
          cyl['conv'] = np.concatenate(([cyl['conv'][0]], cyl['conv']))
          added = True
  else:
    # no parent segment exist
    PC = np.zeros((0,1))
  
  # define the output
  cyl['radius'] = np.atleast_1d(np.atleast_1d(rad)[:numC])
  cyl['length'] = np.atleast_1d(np.atleast_1d(len_)[:numC])
  cyl['start'] = np.atleast_2d(np.atleast_2d(sta)[:numC,:])
  cyl['axis'] = np.atleast_2d(np.atleast_2d(axe)[:numC,:])
  cyl['mad'] = np.atleast_1d(np.atleast_1d(cyl['mad'])[:numC])
  cyl['SurfCov'] = np.atleast_1d(np.atleast_1d(cyl['SurfCov'])[:numC])
  cyl['conv'] = np.atleast_1d(np.atleast_1d(cyl['conv'])[:numC])
  cyl['rel'] = np.atleast_1d(np.atleast_1d(cyl['rel'])[:numC])

  return np.ndarray.flatten(np.atleast_1d(PC)).astype(int).tolist(), cyl, added
