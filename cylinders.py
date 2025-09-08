DEBUG = True
import numpy as np
from tools.verticalcat import verticalcat
from tools.cylinder_fitting import cylinder_fitting

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

def adjustments(cyl, parcyl, inputs, Regs):
  numC = len(cyl['radius'])
  Mod = np.zeros(numC)
  SC = cyl['SurfCov']

  ## Determine the maximum and the minimum radius
  # The maximum based on parent branch

  if np.size(parcyl['radius']>0):
    MaxR = 0.95*parcyl['radius']
    MaxR =np.max((MaxR, inputs['MinCylRad']))
  else:
    # use the maximum from the bottom cylinders
    a = np.min((3,numC))
    MaxR = 1.25*np.max(cyl['radius'][:a])
  MinR = np.min((cyl['radius'][SC>0.7]))
  if np.size(MinR)>0 and np.min(cyl['radius']) < MinR/2:
    MinR = np.min(cyl['radius'][SC>0.4])
  elif np.size(MinR) == 0:
    MinR = np.min(cyl['radius'][SC>0.4])
    if np.size(minR) == 0:
      MinR = inputs['MinCylRad']
  print(f"MinR: {MinR}")
  

  ## Check maximum and minimum radii
  I = cyl['radius'] < MinR
  cyl['radius'][I] = MinR
  Mod[I] = True
  print(inputs.keys())
  if inputs['ParentCor'] or numC <=3:
    I = (cyl['radius'] > MaxR & SC < 0.7) | (cyl['radius']> 1.2*MaxR)
    cyl['radius'][I] = MarR
    Mod[I] = True

    # For short branches modify with more restrictions
    if numC <=3:
      I = ((cyl['radius'] > 0.75 * MaxR) & (SC <0.7))
      if np.any(I):
        r = np.max((SC[I]/0.7*cyl['radius'][I],MinR))
        cyl['radius'][I] = r
        Mod[I] = True
  
  ## Use taper correction to modify radius of too small and large cylinders
  # Adjust radii if a small SurfCov and high SurfCov in the previus and 
  # following cylinders

  for i in range(1,numC-1):
    if (SC[i] < 0.7) and (SC[i-1] >= 0.7) and (SC[i+1] >= 0.7):
      cyl['radius'][i] = 0.5 *(cyl['radius'][i-1] *cyl['radius'][i+1])
      Mod[i] = True
  
  ## Use taper correction to modify radius of too small and large cylinders
  if inputs['TaperCor']:
    if np.max(cyl['radius']) < 0.001:
      ## Adjust radii of thin branches to be linearly decreasing
      if numC>2:
        r = np.sort(cyl['radius'])
        r =r[1:-1]
        a = 2*np.mean(r)
        if a > np.max(r):
          a = np.min((0.01,np.max(r)))
        b=np.min((0.5*np.min(cyl['radius']),0.001))
        cyl['radius'] = np.linspace(a, b, numC)
      elif numC >1:
        r = np.max(cyl['radius'])
        cyl['radius'] = np.column_stack((r, 0.5*r))
      Mod = np.ones(numC, dtype=bool)
    elif numC>4:
      ## Parabola adjustment of maximum and minimum
      # Define parabola taper shape as maximum (and minimum) radii for
      # the cylinders with low surface coverage

      branchlen = np.sum(cyl['length'][:numC]) # branch length
      L = cyl['length']/2+ np.concatenate(([0], np.cumsum(cyl['length'][:numC-1])))
      print(f"L: {L}")
      print(f"branchlen: {branchlen}")
      Taper = np.concatenate((L, [branchlen]))
      Taper = np.column_stack((Taper,np.concatenate((1.05*cyl['radius'], [MinR]))))
      sc = np.concatenate((SC, [1]))
      print('here')
      exit()

      # Least square fitting of parabola to "Taper":
      #A = np.concatenate(([np.sum(sc*Taper[:,0]**4)], [np.sum(sc*Taper[:,1]**2)], [np.sum(sc*Taper[:,1]*(



  
  exit()



  return cyl
def parent_cylinder(SPar, SChi, CiS, cylinder, cyl, si):
  # Finds the parent cylinder from the possible parent segment.
  # Does this by checking it the axis of the cylinder, if continued, will
  # cross the nearby cylinders in the parent segment.
  # Adjust the cylinder so that it starts from the surface of its parent.

  rad = cyl['radius']
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
        print(numC)
        print(si)

        ## Search possible parent cylinder
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



      








  

