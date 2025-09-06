DEBUG = True
import numpy as np
from tools.verticalcat import verticalcat
from tools.surface_coverage_filtering import surface_coverage_filtering
from least_squares_fitting.least_squares_cylinder import least_squares_cylinder
# This file is part of TREEQSM.
#
# TREEQSM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TREEQSM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.
#
#
# ---------------------------------------------------------------------
# CYLINDERS.M       Fits cylinders to the branch-segments of the point cloud
#
# Version 3.0.0
# Latest update     1 Now 2018
#
# Copyright (C) 2013-2018 Pasi Raumonen
# ---------------------------------------------------------------------
#
# Reconstructs the surface and volume of branches of input tree with
# cylinders. Subdivides each segment to smaller regions to which cylinders
# are fitted in least squares sense. Returns the cylinder information and
# in addition the child-relation of the cylinders plus the cylinders in
# each segment.
# ---------------------------------------------------------------------
# Inputs:
# P         Point cloud, matrix
# cover     Cover sets
# segment   Segments
# input     Input parameters of the reconstruction:
#   MinCylRad   Minimum cylinder radius, used in the taper corrections
#   ParentCor   Radius correction based on radius of the parent: radii in
#                 a branch are usually smaller than the radius of the parent
#                 cylinder in the parent branch
#   TaperCor    Parabola taper correction of radii inside branches.
#   GrowthVolCor  If 1, use growth volume correction
#   GrowthVolFac  Growth volume correction factor
#
# Outputs:
# cylinder  Structure array containing the following cylinder info:
#   radius        Radii of the cylinders, vector
#   length        Lengths of the cylinders, vector
#   axis          Axes of the cylinders, matrix
#   start         Starting points of the cylinders, matrix
#   parent        Parents of the cylinders, vector
#   extension     Extensions of the cylinders, vector
#   branch        Branch of the cylinder
#   BranchOrder   Branching order of the cylinder
#   PositionInBranch    Position of the cylinder inside the branch
#   mad           Mean absolute distances of points from the cylinder
#                           surface, vector
#   SurfCov       Surface coverage measure, vector
#   added         Added cylinders, logical vector
#   UnModRadius   Unmodified radii
# ---------------------------------------------------------------------
#
# Changes from version 3.0.0 to 3.1.0, 6 Oct 2021:
# 1) Added the growth volume correction option ("growth_volume_correction")
#    back, which was removed from the previous version by a mistake. The
#    "growth_volume_correction" function was also corrected.
# 2) Added the fields "branch", "BranchOrder", "PositionInBranch" to the
#    output structure "cylinder"
# 3) Removed the fields "CylsInSegment" and "ChildCyls" from the output
#    structure "cylinder"
#
# Changes from version 2.0.0 to 3.0.0, 13 Aug 2020:
# Many comprehensive and small changes:
# 1) "regions" and "cylinder_fitting" are combined into "cylinder_fitting"
#   and the process is more adaptive as it now fits at least 3 (up to 10)
#   cylinders of different lengths for each region.
# 2) "lcyl" and "FilRad" parameters are not used anymore
# 3) Surface coverage ("SurfCov") and mean absolute distance ("mad") are
#   added to the cylinder structure as fields.
# 4) Surface coverage filtering is used in the definition of the regions
#   and removing outliers
# 5) "adjustments" has many changes, particularly in the taper corrections
#   where the parabola-taper curve is fitted to all the data with surface
#   coverage as a weight. Adjustment of radii based on the parabola is
#   closer the parabola the smaller the surface coverage. For the stem the
#   taper correction is the same as for the branches. The minimum and
#   maximum radii corrections are also modified.
# 6) Syntax has changed, particularly for the "cyl"-structure
#
# Changes from version 2.1.0 to 2.1.1, 26 Nov 2019:
# 1) Increased the minimum number "n" of estimated cylinders for
#    initialization of vectors at the beginning of the code. This is done
#    to make sure that trees without branches will not cause errors.
#
# Changes from version 2.0.0 to 2.1.0, 3 Oct 2019:
# 1) Bug fix: UnmodRadius is now defined as it should, as the radius after
#    least squares fitting but without parent, taper or growth vol. corrections
# 2) Bug fix: Correction in "least_squares_cylinder.m", calculates the
#    starting point of the cylinder now correctly.
# 3) Bug fix: Correct errors related to combining data when a fitted
#    cylinder is replaced with two shorter ones, in "cylinder_fitting"
# 4) Removed some unnecessary command lines for computing radius estimates
#    in "regions"
#

def cylinder_fitting(P, Points, Ind, numL, si):
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
      print(f"bot: {bot}")
      print(f"len(Points): {len(Points)}")
      print(f"Ind[i0,0]: {Ind[i0,0]}")
      print(f"Ind[i0+1,1]: {Ind[i0+1,1]}")
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
        if DEBUG:
          print(f"Axis: {Axis}")
          print(f"Top: {Top}")
          print(f"Bot: {Bot}")
        if np.any(np.isnan(Axis).any() ):
          exit()
        c0['axis'] = Axis/np.linalg.norm(Axis)
        #Compute the height along the axis
        if DEBUG:
          print(f"P.shape: {P.shape}")
          print(f"c0['axis'].shape: {c0['axis'].shape}")
          print(f"c0['axis']: {c0['axis']}")
        h = (P[RegC,:]-Bot)@np.transpose(c0['axis'])
        if DEBUG:
          print(f"h: {h}")
        minh = np.min(h)
        # Correct Bot to correspond to the real bottom
        if j==0:
          Bot = Bot + minh*c0['axis']
          c0['start'] = Bot
          h = (P[RegC,:] -Bot)@np.transpose(c0['axis'])
          if DEBUG:
            print(f"h: {h}")
          minh = np.min(h)
        if (i+j)>= numL-1:
          ht = (Top-c0['start'])@np.transpose(c0['axis'])
          if DEBUG:
            print(f"ht: {ht}")
          Top = Top + (np.max(h)-ht)*c0['axis']
        # Compute the height of the Top
        ht = (Top-c0['start'])@np.transpose(c0['axis'])
        if DEBUG:
          print(f"ht: {ht}")
          print(f"h: {h}")
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
            print(f"c0['radius']: {c0['radius']}")
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
            print(I2)
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
        print(RL)

        if again and c['rel'] and c['conv'] and RL >2:
          if si == 1 and c['SurfCov'] > 0.7:
            again= False
          elif si > 1 and c['SurfCov'] > 0.5:
            again = False
      ## Select the best of the fitted cylinders
      # based on maximum surface coverage

      OKfit = (data[:j+1,0] & data[:j+1,1] & data[:j+1, 3] > 1.5)

      J = np.transpose(np.array(range(j)))
      t+=1
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

      if t == 1:
        cyl = c
        names = list(cyl.keys())
        n = len(names)
      else:
        print(f"names: {names}")
        print(f"c.keys(): {c.keys()}")
        print(f"cyl.keys(): {cyl.keys()}")
        print(f"c: {c}")
        print(f"cyl: {cyl}")

        for k in range(n):
          print(f"k: {k}")
          print(f"names[k]: {names[k]}")
          cyl[names[k]] = np.concatenate((np.atleast_1d(cyl[names[k]]), np.atleast_1d(c[names[k]])))

      
      ## compute cylinder top for the definition of the next section
      CylTop = c['start'] + c['length']*c['axis']
    Reg = Reg[list(range(t))]
    print(Reg)
    exit()



  return 1,2 

def cylinders(P, cover, segment, inputs):
  # Initialization of variables
  Segs = segment['segments']
  SPar = segment['parent']
  SChi = segment['children']
  #NumOfSeg = np.max([len(s) for s in Segs])
  NumOfSeg = len(Segs) # number of segments
  if DEBUG:
    print(NumOfSeg)
  n = max(2000, min(NumOfSeg, 2e5))
  c = 1 # number of cylinders determined
  CChi = {}
  CiS = {}
  cylinder={}
  cylinder['radius'] = np.zeros(n, dtype=np.float32)
  cylinder['length'] = np.zeros(n, dtype=np.float32)
  cylinder['start'] = np.zeros(n, dtype=np.float32)
  cylinder['axis'] = np.zeros(n, dtype=np.float32)
  cylinder['parent'] = np.zeros(n, dtype=np.uint32)
  cylinder['extension'] = np.zeros(n, dtype=np.uint32)
  cylinder['added'] = np.zeros(n, dtype=bool)
  cylinder['UnmodRadius'] = np.zeros(n, dtype=np.float32)
  cylinder['branch'] = np.zeros(n, dtype=np.uint16)
  cylinder['SurfCov'] = np.zeros(n, dtype=np.float32)
  cylinder['mad'] = np.zeros(n, dtype=np.float32)


  ## Determine suitable order of segments (from trunk to the "youngest" child)
  bases = np.arange(NumOfSeg, dtype=int)
  if DEBUG:
    print(f"SChi: {SChi}")

    print(f"SPar[:,0] == 0: {SPar[:,0] == 0}")
  bases = bases[SPar[:,0] == 0]
  numB = len(bases)
  SegmentIndex= np.zeros(NumOfSeg, dtype=int)
  if DEBUG:
    print(f"bases: {bases}")
    print(f"len(bases): {len(bases)}")
    print(f"SegmentIndex: {SegmentIndex}")
    print(f"len(SegmentIndex): {len(SegmentIndex)}")
  numC = 0
  for i in range(numB):
    if numC == NumOfSeg:
      break
    numC = numC+1
    SegmentIndex[numC -1] = bases[i]
    if DEBUG:
      print(f"SegmentIndex(in for): {SegmentIndex}")
      print(f"bases(in for): {bases}")
    S = [schi for schi in SChi[bases[i]] if schi != 0] #TODO: Look for a way to not include 0 as child.
    
    while S:
      n = len(S)
      SegmentIndex[numC-1+1:numC-1+n+1] = S
      if DEBUG:
        print(f"SegmentIndex(in while): {SegmentIndex}")
      numC+= n
      S = [schi for s in S if s != 0 for schi in SChi[s] ] #TODO: Look for a way to not include 0 as child.
  if DEBUG:
    print(f"SegmentIndex: {SegmentIndex}")

  ## Fit cylinders individually for each segment
  for k in range(NumOfSeg):
    si = SegmentIndex[k]
    if si > -1:
      ## Some initialization about the segment
      Seg = Segs[si] # the current segment under analysis
      numL = len(Seg) # number of cover set layers in the segment
      Sets, IndSets = verticalcat(Seg) # the cover sets in the segment

      numS = len(Sets) # number of cover sets in the current segment
      Points = np.concatenate([ cover['ball'][s] for s in Sets]) # the points of the segments
      numP = len(Points) # number of points in the segment

      # Determine indices of points for faster definition of regions
      #BallSize = [ np.size(cover['ball'][s]) for s in Sets ]
      BallSize = [ np.size(cover['ball'][s]) for s in Sets ]
      IndPoints = np.zeros((numL,2), dtype=int) # indices for points in each layer of the segment
      for j in range(numL):
        IndPoints[j,1] = np.sum(BallSize[IndSets[j,0]:IndSets[j,1]])
      IndPoints[:,1] = np.cumsum(IndPoints[:,1])
      IndPoints[1:,0] = IndPoints[1:,0] + IndPoints[:-1, 1]
      Base = Seg[0] # the base of the segment
      numB = IndPoints[0,1] # number of points in the base

      # Reconstruct only large enough segments
      if (numL > 1) and (numP > numB) and (numS > 2) and (numP > 20) and len(Base) > 0:
        cyl, Reg = cylinder_fitting(P, Points, IndPoints, numL, si)
        exit()

      








  

