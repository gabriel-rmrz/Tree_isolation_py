DEBUG = True
import numpy as np
import cmath, math
from tools.verticalcat import verticalcat
from tools.adjustments import adjustments
from tools.cylinder_fitting import cylinder_fitting
from tools.mat_vec_subtraction import mat_vec_subtraction

#from tools.surface_coverage_filtering import surface_coverage_filtering
#from least_squares_fitting.least_squares_cylinder import least_squares_cylinder
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
  if len(SPar[si]) > 0: # parent segment exist, find parent cylinder
    s = SPar[si]
    PC = [CiS[s_] for s_ in s if s_ in CiS]  # the cylinders in the parent segment
    # select the closest cylinders for closer examination
    if len(PC) > 1:
      D = mat_vec_sumstraction(-cylinder['start'][Pc,:], -sta[1,:])
      print(D)
      d = np.sum(D*D,1)
      if len(PC) > 3:
        I = I[:4]
      pc = PC[i]
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
        if numC > 1


    exit()
  return

def cylinders(P, cover, segment, inputs):
  # Initialization of variables
  Segs = segment['segments']
  SPar = segment['parent']
  SChi = segment['children']
  #NumOfSeg = np.max([len(s) for s in Segs])
  NumOfSeg = len(Segs) # number of segments
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
  bases = bases[SPar[:,0] == 0]
  numB = len(bases)
  SegmentIndex= np.zeros(NumOfSeg, dtype=int)
  numC = 0
  for i in range(numB):
    if numC == NumOfSeg:
      break
    numC = numC+1
    SegmentIndex[numC -1] = bases[i]
    S = [schi for schi in SChi[bases[i]] if schi != 0] #TODO: Look for a way to not include 0 as child.
    
    while S:
      n = len(S)
      SegmentIndex[numC-1+1:numC-1+n+1] = S
      numC+= n
      S = [schi for s in S if s != 0 for schi in SChi[s] ] #TODO: Look for a way to not include 0 as child.

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

        numC = int(np.size(np.atleast_1d(cyl['radius'])))

        ## Search possible parent cylinder
        print(f"numC: {numC}")
        print(f"si: {si}")
        if numC>0 and si > 0:
          parent_cylinder(SPar,SChi,CiS, cylinder,cyl,si)

        elif si == 0:
          PC = []
          added = False
        else:
          added = False
        
        cyl['radius0'] = cyl['radius']

        ## Modify cylinders
        if numC>0:
          # Define parent cylinder
          parcyl = {}
          parcyl['radius'] = cylinder['radius'][PC]
          parcyl['length'] = cylinder['length'][PC]
          parcyl['start'] = cylinder['start'][PC]
          parcyl['axis'] = cylinder['axis'][PC]
          # Modify the cylinders
          cyl = adjustments(cyl,parcyl,inputs,Reg)




      








  

