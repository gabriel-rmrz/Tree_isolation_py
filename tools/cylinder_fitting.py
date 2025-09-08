import numpy as np
from tools.surface_coverage_filtering import surface_coverage_filtering
from least_squares_fitting.least_squares_cylinder import least_squares_cylinder


def cylinder_fitting(P, Points, Ind, numL, si):
  #c0={}
  if numL > 6:
    i0 = 0
    i = 3 # indices of the first and las layers of the region
    t = 0
    Reg = {}
    cyls = {}
    regs = {}
    data = np.zeros((11,4), dtype=int)
    while i0 < numL - 3:
      ## Fit at least three cylinders of different lengths
      bot = Points[Ind[i0,0]:Ind[i0+1,1]+1]
      Bot = np.average(P[bot,:], axis=0)
      again = True
      j = 0
      c0={}
      #CylTop = c['start']+c['length']+c['axis']
      while (i+j <=numL-1) and j<= 10 and (j<=2 or again):
        ## Select points and estimate axis
        RegC = Points[Ind[i0,0]:Ind[i+j,1]+1]  #candidate region
        # Top axis point of the region:
        top = Points[Ind[i+j-1,0]:Ind[i+j,1]+1]
        Top = np.average(P[top,:], axis=0)
        Axis = Top-Bot
        if np.any(np.isnan(Axis).any() ):
          exit()
        c0['axis'] = Axis/np.linalg.norm(Axis)
        #Compute the height along the axis
        h = (P[RegC,:]-Bot)@np.transpose(c0['axis'])
        minh = np.min(h)
        # Correct Bot to correspond to the real bottom
        if j==0:
          Bot = Bot + minh*c0['axis']
          c0['start'] = Bot
          h = (P[RegC,:] -Bot)@np.transpose(c0['axis'])
          minh = np.min(h)
        if (i+j)>= numL-1:
          ht = (Top-c0['start'])@np.transpose(c0['axis'])
          Top = Top + (np.max(h)-ht)*c0['axis']
        # Compute the height of the Top
        ht = (Top-c0['start'])@np.transpose(c0['axis'])
        Sec = (h <= ht) & (h>=minh) # only points below the Top
        c0['length'] = ht-minh # length of the region/cylinder
        # The region for the cylinder fitting:
        reg = RegC[Sec]
        Q0 = P[reg,:]
        

        ## Filter points and estimate radius
        if len(Q0) > 20:
          Keep, c0 = surface_coverage_filtering(Q0,c0, 0.02,20)
          reg = reg[Keep]
          Q0 = Q0[Keep,:]
        else:
          c0['radius'] = 0.01
          c0['SurfCov'] = 0.05
          c0['mad'] = 0.01
          c0['conv'] = 1
          c0['rel'] = 1

        ## Fit cylinder
        if len(Q0) > 9:
          if i >= numL-1 and t ==0:
            c = least_squares_cylinder(Q0, c0)
          elif i >=numL-1 and t>0:
            h = (Q0 - CylTop)@np.transpose(c0['axis'])
            I = h >= 0
            Q = Q0[I,:] # the section
            reg = reg[I]
            n2 = len(Q)
            n1 = np.count_nonzero(~I)
            if (n2 > 9) and (n2 >5):
              Q0 = np.concatenate((Q[~I,:],Q)) # the point cloud for cylinder fitting
              W = np.concatenate((1/3*np.ones(n2), 2/3*np.ones(n1))) # the weights
              c = least_squares_cylinder(Q0,c0,W,Q)
            else:
              c = least_squares_cylinder(Q0,c0)
          elif t == 0:
            top = Points[Ind[i+j-3,0]:Ind[i+j-2,1]+1]
            Top = np.average(P[top,:])
            ht = (Top-Bot)@np.transpose(c0['axis'])
            h = (Q0-Bot)@np.transpose(c0['axis'])
            I = (h<=ht)
            Q = Q0[I,:] # the section
            reg = reg[I]
            n2 = len(Q)
            n3 = np.count_nonzero(~I)
            if (n2 > 9) and (n3 > 5):
              Q0 = np.concatenate((Q, Q0[I,:])) # the point cloud for cylinder fitting
              W = np.concatenate((2/3*np.ones(n2), 1/3*np.ones(n3))) # the weights
              c = least_squares_cylinder(Q0, c0, W, Q)
            else:
              c = least_squares_cylinder(Q0,c0)
          else:
            top = Points[Ind[i+j-3,0]:Ind[i+j-2,1]]
            Top = np.average(P[top,:]) # Top axis point of the region
            ht = (Top - CylTop)@np.transpose(c0['axis'])
            h = (Q0-CylTop)@np.transpose(c0['axis'])
            I1 = (h < 0) # the bottom
            I2 = (h >= 0) & (h<=ht) # the section
            I3 = (h > ht) # the top
            Q = Q0[I2, :]
            reg = reg[I2]
            n1 = np.count_nonzero(I1)
            n2 = len(Q)
            n3= np.count_nonzero(I3)
            if n2 > 9:
              Q0 = np.concatenate((Q0[I1,:], Q, Q0[I3,:]))
              W = np.concatenate((1/4*np.ones(n1), 2/4+np.ones(n2), 1/4*np.ones(n3)))
              c = least_squares_cylinder(Q0, c0, W,Q)
            else:
              c = c0
              c['rel'] = 0
          if c['conv'] == 0:
            c = c0
            c['rel'] =0
          if c['SurfCov'] < 0.2:
            c['rel'] = 0
        else:
          c = c0
          c['rel'] = 0

        # Collect fit data
        data[j,:] = np.array([c['rel'], c['conv'], c['SurfCov'], c['length']/c['radius']])
        cyls[j] = c
        regs[j] = reg
        j+=1
        # If reasonable cylinder fitted, then stop fitting new ones
        # (but always fit at least cylinders)
        RL = c['length']/c['radius'] # relative length of the cylinder

        if again and c['rel'] and c['conv'] and RL >2:
          if si == 1 and c['SurfCov'] > 0.7:
            again= False
          elif si > 1 and c['SurfCov'] > 0.5:
            again = False
      ## Select the best of the fitted cylinders
      # based on maximum surface coverage

      OKfit = (data[:j+1,0] & data[:j+1,1] & data[:j+1, 3] > 1.5)

      J = np.transpose(np.array(range(j)))
      if np.any(OKfit):
        J = J[OKfit]
      I = np.argmax(data[J,2] - 0.01*data[J,3])
      J = J[I]
      c = cyls[J]

      ## Update the indices of the layers for the next region:
      CylTop = c['start'] + c['length']*c['axis']
      i0+=1
      bot = Points[Ind[i0,0]:Ind[i0,1]+1]
      Bot = np.average(P[bot,:]) # Bottom axis point of the region
      h = (Bot - CylTop)@np.transpose(c['axis'])
      i00 = i0
      while i0+1 < numL-1 and i0 < i00+5 and h < -c['radius']/3:
        i0+=1
        bot = Points[Ind[i0, 0]:Ind[i0+1,1]+1]
        Bot = np.average(P[bot,:]) # Bottom axis point of the region
        h = (Bot-CylTop)@np.transpose(c['axis'])
      i = i0 + 5
      i = np.min((i,numL-1))

      ## If the next section is very short part of the end of the branch
      # then simply increase the length of the current cylinder
      if numL -1 - i0 +2 < 4:
        reg = Points[Ind[numL-1-5, 0]: Ind[numL-1,1]+1]
        Q0 = P[reg,:]
        ht = (c['start'] + c['length']*c['axis'])@np.transpose(c['axis'])
        h = Q0@np.transpose(c['axis'])
        maxh = np.max(h)
        if maxh > ht:
          c['length'] = c['length'] + (maxh -ht)
        i0 = numL -1
      Reg[t] = regs[j-1]
      t+=1

      if t == 1:
        cyl = c
        names = list(cyl.keys())
        n = len(names)
      else:
        for k in range(n):
          cyl[names[k]] = np.column_stack((np.atleast_1d(cyl[names[k]]), np.atleast_1d(c[names[k]])))

      
      ## compute cylinder top for the definition of the next section
      CylTop = c['start'] + c['length']*c['axis']
    Reg = {t_:Reg[t_] for t_ in range(t)}
    for k in cyl.keys():
      cyl[k] = np.transpose(cyl[k])
      if np.min(np.shape(cyl[k])) == 1:
        cyl[k] = np.ndarray.flatten(cyl[k])
  else:
    c0 = {}
    ## Define the region for small segments
    Q0 = P[Points,:]#
    if len(Q0) > 10:
      ## Define the direction
      bot = Points[Ind[0,0]:Ind[0,1]]
      Bot = np.average(P[bot,:],axis=0)
      top = Points[Ind[numL-1, 0]:Ind[numL-1,1]]
      Top = np.average(P[top,:],axis=0)
      Axis = Top-Bot
      if np.any(np.isnan(Axis).any() ):
        print(f"Axis: {Axis}")
        exit()
      c0['axis'] = Axis/np.linalg.norm(Axis)
      h = Q0@np.transpose(c0['axis'])
      c0['length'] = np.max(h) - np.min(h)
      hpoint = Bot@np.transpose(c0['axis'])
      c0['start'] = Bot - (hpoint - np.min(h))*c0['axis']

      ## Define other ouputs
      Keep, c0 = surface_coverage_filtering(Q0,c0, 0.02, 20)
      Reg = {}
      Reg[0] = Points[Keep]
      Q0 = Q0[Keep, :]
      cyl = least_squares_cylinder(Q0, c0)
      if ~cyl['conv'] or ~cyl['rel']:
        cyl = c0
      t = 1
    else:
      cyl = 0
      t = 0
  # Define Reg as coordinates
  Reg = {i: P[Reg[i],:] for i in range(t)}

  return cyl, Reg
