DEBUG = True
import numpy as np
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
    print(bases)

    print(f"bases: {bases}")
    print(f"SPar[:,0] == 0: {SPar[:,0] == 0}")
  bases = bases[SPar[:,0] == 0]
  if DEBUG:
    print(bases)
  numB = len(bases)
  SegmentIndex= np.zeros(NumOfSeg, dtype=int)
  numC = 0
  for i in range(numB):
    numC = numC+1
    SegmentIndex[numC -1] == bases[i]
    if DEBUG:
      print(f"SChi: {SChi}")
      print(f"bases[i]: {bases[i]}")
      for schi in SChi[bases[i]]:
        print(f"schi: {schi}")
    S = [schi for schi in SChi[bases[i]]]
    while len(S):
      n = len(S)
      SegmentIndex[numC-1+1:numC-1+n+1] = S
      numC+= n
      S = [schi for s in S for schi in SChi[s]]
  print(S)









  

