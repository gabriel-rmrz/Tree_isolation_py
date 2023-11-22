DEBUG=False
import numpy as np
from collections import defaultdict

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

def connected_components(Nei,Sub,MinSize,Fal=-1):
  if DEBUG:
    print(f"len(Nei): {len(Nei)}")
    print(f"len(Sub): {len(Sub)}")
    exit()
  '''
  ---------------------------------------------------------------------
  CONNECTED_COMPONENTS.M      Determines the connected components of cover
                                    sets using their neighbour-relation
 
  Version 1.2.0
  Latest update     4 Jun 2018
 
  Copyright (C) 2013-2017 Pasi Raumonen
  ---------------------------------------------------------------------
 
  Determines connected components of the subset of cover sets defined
  by "Sub" such that each component has at least "MinSize"
  number of cover sets.
 
  Inputs:
  Nei       Neighboring cover sets of each cover set, (n_sets x 1)-cell
  Sub       Subset whose components are determined,
                length(Sub) < 2 means no subset and thus the whole point cloud
                "Sub" may be also a vector of cover set indexes in the subset
                or a logical (n_sets)-vector, where n_sets is the number of
                all cover sets
  MinSize   Minimum number of cover sets in an acceptable component
  Fal       Logical false vector for the cover sets
 
  Outputs:
  Components    Connected components, (n_comp x 1)-cell
  CompSize      Number of sets in the components, (n_comp x 1)-vector
 
  Changes from version 1.1.0 to 1.2.0, 4 Jun 2018:
  1) Corrected the CompSize for cases n = 2 and n = 3.
  '''
  Fal_i = Fal
  Components = defaultdict()
  CompSize = [] 
  if (len(Sub)==0) or (len(Nei)==0):
    return Components, CompSize
  if (len(Sub)<=3) and not isinstance(Sub[0], (np.bool_, bool)) and (Sub[0] > 0) :
    n = len(Sub)
    if n == 1:
      Components[1] = int(Sub)
      CompSize = [1]
    elif n == 2:
      I = Nei[Sub[0]] == Sub[1]
      if sum(I):
        Components[1] = int(Sub)
        CompSize = [2]
      else:
        Components[1] = int(Sub[0])
        Components[2] = int(Sub[1])
        CompSize = [1, 1]
    elif n == 2:
      I = Nei[Sub[0]] == Sub[1]
      J = Nei[Sub[0]] == Sub[2]
      K = Nei[Sub[1]] == Sub[2]
      if sum(I) + sum(J) + sum(K) >= 2:
        Components[1] = int(Sub)
        CompSize[3]
      elif sum(I):
        Components[1] = int(Sub[:2])
        Components[2] = int(Sub[2])
        CompSize = [2, 1]
      elif sum(J):
        Components[1] = int(Sub[[0,2]])
        Components[2] = int(Sub[1])
        CompSize = [2, 1]
      elif sum(K):
        Components[1] = int(Sub[[1,2]])
        Components[2] = int(Sub[0])
        CompSize = [2, 1]
      else:
        Components[1] = int(Sub[0])
        Components[2] = int(Sub[1])
        Components[3] = int(Sub[2])
        CompSize = [1, 1, 1]
  elif sum(Sub) or (len(Sub) == 1 and Sub[0] == 0):
    numB = len(Nei)
    numS = []
    if Fal_i == -1:
      Fal = np.zeros(numB, dtype='bool')
    if (len(Sub) == 1) and (Sub[0] ==0):
      # All the cover sets
      numS = numB
      if Fal_i == -1:
        Sub = np.ones(numB, dtype='bool')
      else:
        Sub = ~Fal
    elif not isinstance(Sub[0], (np.bool_, bool)):
      # Subset of cover sets
      numS = len(Sub)
      if Fal_i == -1:
        sub = np.zeros(numB, dtype='bool')
      else:
        sub = Fal
      sub[Sub] = True
      Sub = sub
    else:
      # Subset of cover sets
      numS = len(Sub[Sub])


    CompSize = np.zeros(numS, dtype=int)
    numC = 0
    m = 1
    while ~Sub[m-1]:
      m = m+1
    i = 0
    Comp = np.zeros(numS, dtype=int)

    while i < numS and m <= numB:
      Add = Nei[m]
      I = Sub[Add-1]
      a = len(Add)
      Comp[0] = m
      Sub[m-1] = False
      t = 1
      while a>0:
        Comp[t:t+a] = Add
        Sub[Add-1] = False
        t = t+a
        Add_temp = []
        for i in Add:
          Add_temp = np.concatenate((Add_temp, Nei[i]), axis=0) 
        Add = np.asarray(Add_temp).astype(int)
        print(f"Add: {Add}")
        I = Sub[Add]
        Add = Add[I]
        n = len(Add)
        if n > 2:
          I = np.ones(n, dtype='bool')
          for j in range(1,n+1):
            if ~Fal[Add[j-1]-1]:
              Fal[Add[j-1]-1] = True
            else:
              I[j-1] = False
          Fal[Add-1] = False
          Add = Add[I]
        elif n == 2:
          if Add[0] == Add[1]:
            Add = Add[0]
        a = len(Add)
      i = i +t
      if t>= MinSize:
        print(f"Comp: {Comp}")
        numC = numC +1
        print(Comp[:t+1])
        print(Comp)
        print(Components)
        print(CompSize)
        Components[numC] = int(Comp[:t])
        CompSize[numC] = t
        print(Components)
        print(CompSize)
      if i < numS:
        while m <= numB and Sub[m-1] == False:
          m = m + 1
    print(Components)
    print(CompSize)
    CompSize = CompSize[:numC]

  else:
    Components = {}
    CompSize = []
  
  return Components, CompSize
  
    

