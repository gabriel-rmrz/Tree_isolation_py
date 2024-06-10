DEBUG= True
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

def connected_components(Nei,Sub,MinSize,Fal=None):
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
  Fal_i = np.copy(Fal)
  if (len(Sub)==0) or (len(Nei)==0):
    return {}, [] 
  if DEBUG:
    print(Sub[0])
  if (len(Sub)<=3) and not isinstance(Sub[0], (np.bool_, bool)) and ((Sub[0] + 1) > 0) : # Added the + 1 to Sub[0] because we are using indices starting from 0.
    n = len(Sub)
    Components = {}
    if n == 1:
      Components[0] = np.array(Sub, dtype=np.uint32)
      CompSize = [1]
      return Components, CompSize
    elif n == 2:
      I = Nei[Sub[0]] == Sub[1]
      if np.any(I):
        Components[0] = np.array(Sub, dtype=np.uint32)
        CompSize = [2]
      else:
        Components[0] = np.array(Sub[0], dtype=np.uint32)
        Components[1] = np.array(Sub[1], dtype=np.uint32)
        CompSize = [1, 1]
      return Components, CompSize
    elif n == 3:
      I = Nei[Sub[0]] == Sub[1]
      J = Nei[Sub[0]] == Sub[2]
      K = Nei[Sub[1]] == Sub[2]
      if np.any(I) + np.any(J) + np.any(K) >= 2:
        Components[0] = np.array(Sub, dtype=np.uint32)
        CompSize = [3]
      elif np.any(I):
        Components[0] = np.asarray(Sub[:2], dtype=np.uint32)
        Components[1] = np.array(Sub[2], dtype=np.uint32)
        CompSize = [2, 1]
      elif np.any(J):
        Components[0] = np.asarray(Sub[[0,2]], dtype=np.uint32)
        Components[1] = np.array(Sub[1], dtype=np.uint32)
        CompSize = [2, 1]
      elif np.any(K):
        Components[0] = np.asarray(Sub[[1,2]], dtype=np.uint32)
        Components[1] = np.array(Sub[0], dtype=np.uint32)
        CompSize = [2, 1]
      else:
        Components[0] = np.array(Sub[0], dtype=np.uint32)
        Components[1] = np.array(Sub[1], dtype=np.uint32)
        Components[2] = np.array(Sub[2], dtype=np.uint32)
        CompSize = [1, 1, 1]
      return Components, CompSize
  elif np.any(Sub) or (len(Sub) == 1 and (Sub[0] + 1) == 0):

    numB = len(Nei)
    numS = 0 
    if Fal_i == None:
      Fal = np.zeros(numB, dtype='bool')
    if (len(Sub) == 1) and (Sub[0] + 1 ==0):
      # All the cover sets
      numS = numB
      if Fal_i == None:
        Sub = np.ones(numB, dtype='bool')
      else:
        Sub = ~Fal
    elif not isinstance(Sub[0], (np.bool_, bool)):
      # Subset of cover sets
      numS = len(Sub)
      if Fal_i == None:
        sub = np.zeros(numB, dtype='bool')
      else:
        sub = Fal
      sub[Sub] = True
      Sub = sub
    else:
      # Subset of cover sets
      numS = np.count_nonzero(Sub)
      if DEBUG:
        print(f"We are here!")
        print(f"Counting non zero values of sub")
        print(f"numS: {numS}")
        print(f"Sub: {Sub}")



    Components = {}#defaultdict()
    CompSize = np.zeros(numS, dtype=np.uint32)
    numC = 0
    m = 0
    while not Sub[m]:
      m = m+1
    i = 0
    Comp = np.zeros(numS, dtype=np.uint32)

    while i < numS and m < numB:
      Add = Nei[m]
      I = Sub[Add]
      Add = Add[I]
      if type(Add) != np.uint32:
        a = len(Add)
      else:
        a = 1
      Comp[0] = m
      Sub[m] = False
      t = 1
      while a>0:
        #print(f"len(Add): {len(Add)}")
        #print(f"len(Comp[t:t+a]): {len(Comp[t:t+a])}")
        Comp[t:t+a] = Add
        Sub[Add] = False
        t += a
        '''
        Add_temp = []
        for i in Add:
          Add_temp = np.concatenate((Add_temp, Nei[i]), axis=0) 
        Add = np.asarray(Add_temp).astype(int)
        '''
        if type(Add) != np.uint32:
          Add = np.concatenate([Nei[key] for key in Add])
        else:
          Add = Nei[Add] 
        '''
        if DEBUG:
          print(f"Add: {Add}")
        '''
        I = Sub[Add]
        Add = Add[I]
        if type(Add) != np.uint32:
          n = len(Add)
        else:
          n = 1
        if n > 2:
          I = np.ones(n, dtype='bool')
          for j in range(n):
            if not Fal[Add[j]]:
              Fal[Add[j]] = True
            else:
              I[j] = False
          Fal[Add] = False
          Add = Add[I] 
        elif n == 2:
          if Add[0] == Add[1]:
            Add = Add[0]
        if type(Add) != np.uint32:
          a = len(Add)
        else:
          a = 1


      i += t
      if t>= MinSize:
        Components[numC] = np.array(Comp[:t], dtype=np.uint32)
        CompSize[numC] = t
        numC += 1
      if i < numS:
        while m < numB and Sub[m] == False:
          m += 1

    '''
    if DEBUG:
      print(f"len(Components): {len(Components)}")
    Components = {key:Components[key] for key in range(numC)}
    '''
    CompSize = CompSize[:numC]
    return Components, CompSize

  else:
    Components = {}
    CompSize = []
  
  return Components, CompSize
  
    

