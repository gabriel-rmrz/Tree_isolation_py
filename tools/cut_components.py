DEBUG=False
import numpy as np
from tools.unique_elements import unique_elements
def cut_components(Nei,Cut,CutSize,Fal, isFal):
  Fal = np.copy(Fal)
  isFal = np.copy(isFal)
  Cut = np.copy(Cut)
  # Define the connected components of the Cut
  if CutSize == [1]:
    # Cut is connected and therfore Study is also
    CompSize = [1]
    Components = {0: np.array(Cut).astype(int)}
    return Components, CompSize
  elif CutSize == [2]:
    I = Nei[Cut[0]] = Cut[1]
    if np.any(I):
      Components = {0:np.array(Cut).astype(int)}
      CompSize = [1]
    else:
      Components = {0:np.array(Cut[0]).astype(int), 1:np.array(Cut[1]).astype(int)}
      CompSize = [1, 1]
    return Components, CompSize
  elif CutSize == 3:
    I = Nei[Cut[0]] == Cut[1]
    J = Nei[Cut[0]] == Cut[2]
    K = Nei[Cut[1]] == Cut[2]
    if np.any(I) + np.any(J) + np.any(K) >=2:
      CompSize = [1]
      Components = {0:np.array(Cut).astype(int)}
    elif np.any(I):
      Components = {0:np.asarray(Cut[:2]).astype(int), 1:np.array(Cut[2]).astype(int)}
      CompSize = [2,1]
    elif np.any(J):
      Components = {0:np.asarray(Cut[[0,2]]).astype(int), 1:np.array(Cut[1]).astype(int)}
      CompSize = [2,1]
    elif np.any(K):
      Components = {0:np.asarray(Cut[[1,2]]).astype(int), 1:np.array(Cut[0]).astype(int)}
      CompSize = [2,1]
    else:
      CompSize = [1, 1, 1]
      Components = {0:np.asarray(Cut[0]).astype(int), 1:np.array(Cut[1]).astype(int), 2:np.array(Cut[2]).astype(int)}
    return Components, CompSize
  else: 
    Components = {}
    CompSize = np.zeros(CutSize, dtype=np.int32)
    Comp = np.zeros(CutSize, dtype=np.int32)
    Fal[Cut] = True
    numC = 0
    m = Cut[0]
    i = 0
    while i < CutSize:
      Added = np.copy(Nei[m])
      I = Fal[Added]
      Added = Added[I]
      a = len(Added)
      Comp[0] = m
      Fal[m] = False
      t = 1
      while a > 0:
        #Comp[t+1-1:t+a] = Added
        '''
        if (len(Comp)+1 == a + t) and (len(Comp[t+1-1:t+a]) < len(Added)):
          Comp[t+1-1:t+a] = Added[:a-1]
          np.append(Comp,Added[a-1])
        else:
          Comp[t+1-1:t+a] = Added
        '''
        if len(Comp[t+1-1:t+a]) < len(Added):
          diff = len(Added) - len(Comp[t+1-1:t+a])
          Comp[t+1-1:t+a] = Added[:a-diff]
          for d in range(1, diff +1):
            np.append(Comp, Added[a-d])
        else:
          Comp[t+1-1:t+a] = Added
          
        Fal[Added] = False
        t = t + a
        Ext = np.concatenate([Nei[key] for key in Added])
        Ext = unique_elements(Ext,isFal)
        #Ext = np.unique(Ext)
        I = Fal[Ext]
        Added = Ext[I]
        a = len(Added)
      i += t
      numC +=1
      Components[numC-1] = Comp[:t]
      CompSize[numC-1] = t
      if i < CutSize:
        J = Fal[Cut]
        m = Cut[J]
        m = m[0]
    CompSize = CompSize[:numC]
    return Components, CompSize

