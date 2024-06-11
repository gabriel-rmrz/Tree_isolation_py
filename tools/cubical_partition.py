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
def cubical_partition(P, EL, NE=3, method=1):
  '''
  Partitions the input point cloud into cubes.
  Inputs:
  P             Point cloud, (n_points x 3)-matrix
  EL            Length of the cube edges
  NE            Number of empty edge layers
  
  Outputs:
  Partition     Point cloud partitioned into cubical cells,
                  (nx x ny x nz)-cell, where nx,ny,nz are the number
                  of cubes in x,y,z-directions, respectively
  CC            (n_points x 3)-matrix whose rows are the cube coordinates
                  of each point: x,y,z-coordinates
  Info          Minimum coordinate values, number of cubes in each
                  coordinate direction, EL and NE.
  Cubes         Partition contains only the non-empty cells as a 1-d array
                and Cubes gives each cell its array-index
  '''


  # The vertice of the big cube containing P
  Min = np.min(P, axis=0)
  Max = np.max(P, axis=0)

  # Number of cubes with the edge length "EdgeLength" in the sides
  # of the big cube
  N = np.ceil((Max-Min)/EL)+2.*NE+1.

  while 8*N[0]*N[1]*N[2] > 4e9:
    EL = 1.1*EL
    N = np.ceil((Max-Min)/EL)+2*NE+1

  Info = np.concatenate([Min, N, [EL, NE]])
  #Info = [Min, N, EL, NE]

  # Calculated the cube-coordinates of the points
  CubeCoord = np.floor((P-Min)/EL)+NE+1

  # Sorts the points according a lexicographical order
  LexOrd = CubeCoord[:,0] + N[0]* (CubeCoord[:,1] -1) + N[0]*N[1]*(CubeCoord[:,2] -1)
  CubeCoord = CubeCoord.astype(np.uint16)
  SortOrd = (np.argsort(LexOrd)).astype(np.uint32)
  LexOrd = (LexOrd[SortOrd]).astype(np.uint32)

  # Define "Partition"
  Partition = {}
  numP = len(P[:,0]) # Number of points
  if method== 1:
    p = 1 # The index of the point under comparison
    while p <= numP:
      t = 1
      while (p+t <= numP) and (LexOrd[p-1] == LexOrd[p-1+t]):
        t += 1
      q = SortOrd[p-1]
      Partition[tuple([CubeCoord[q,0],CubeCoord[q,1], CubeCoord[q,2]])] = SortOrd[p-1:p+t-1]
      p += t
    return Partition, CubeCoord, Info
  else:
    numC = len(np.unique(LexOrd))

    # Define "Partition"
    Cubes = np.zeros((N[0], N[1], N[2]), dtype=np.uint32)

    p = 1
    c = 0
    while p <= numP:
      t = 1
      while (p+t <= numP) and (LexOrd[p-1] == LexOrd[p-1+t]):
        t +=1
      q = SortOrd[p-1]
      c += 1
      Partition[c] = SortOrd[p-1:p+t-1]
      Cubes[CubeCoord[q,0],CubeCoord[q,1], CubeCoord[q,2]] = c
      p += t
    return Partition, CubeCoord, Info, Cubes


