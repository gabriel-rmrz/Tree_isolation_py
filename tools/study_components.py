DEBUG = False
import numpy as np
from tools.unique_elements import unique_elements

def study_components(Nei, numS, Cut, CutComps, Forb, Fal, isFalse):
  Cut = np.copy(Cut)
  #CutComps = np.copy(CutComps)
  Forb = np.copy(Forb)
  Fal = np.copy(Fal)
  isFalse = np.copy(isFalse)
  Study = {}
  StudySize = np.zeros(numS,dtype=np.int32)
  Study[0] = Cut
  StudySize[0] = len(Cut)
  if numS >= 2:
    N = np.copy(Cut)
    i = 1
    while i < numS:
      Forb[N] = True
      N = np.concatenate([Nei[key] for key in N])
      N = unique_elements(N,Fal)
      I = Forb[N]
      N = N[~I]
      if len(N) > 0:
        i +=1
        Study[i-1] = N
        StudySize[i-1] = len(N)
      else:
        Study = {key: Study[key] for key in range(i)} 
        StudySize = [StudySize[key] for key in range(i)]
        i = numS +1
  
  # Define study as a vector
  numS = len(StudySize)
  studysize = sum(StudySize)
  study = np.concatenate([Study[key] for key in Study.keys()])

  # Determine the components of study
  numC = len(CutComps)
  i = 1 # index of cut component
  j = 0 # number of elements attributed to components
  k = 0 # number of study components
  Fal[study] = True
  Components = {}
  CompSize = np.zeros(numC, dtype=np.int32)
  Comp = np.zeros(studysize, dtype=np.int32)
  while i <= numC:
    C = CutComps[i-1]
    C = np.array(C)
    while j < studysize:
      if DEBUG:
        print(f"type(C): {type(C)}")
      if isinstance(C, np.ndarray):
        a = C.size
      else:
        a = 0
      Comp[:a] = C
      Fal[C] = False
      if DEBUG:
        print(f"C: {C}")
      if a>1:
        Add = unique_elements(np.concatenate([Nei[key] for key in C]), isFalse)
      else:
        Add =  Nei[C.item()] 
      '''
      elif a ==1 and len(C.shape) == 1:
        if DEBUG:
          print(f"C: {C}")
          print(f"C.item(): {C.item()}")
          print(f"type(C): {type(C)}")
          print(f"isinstance(C, np.ndarray): {isinstance(C, np.ndarray)}")
        Add =  Nei[C[0]] 
      '''
      t =a
      I = Fal[Add]
      Add = Add[I]
      a = len(Add)
      while a>0:
        if len(Comp[t+1-1:t+a]) < len(Add):
          diff = len(Add) - len(Comp[t+1-1:t+a])
          Comp[t+1-1:t+a] = Add[:a-diff]
          for d in range(1, diff +1):
            np.append(Comp, Add[a-d])
        else:
          Comp[t+1-1:t+a] = Add
        '''
        if (len(Comp)+1 == a + t) and (len(Comp[t+1-1:t+a]) < len(Add)):
          Comp[t+1-1:t+a] = Add[:a-1]

          np.append(Comp,Add[a-1])
        else:
          Comp[t+1-1:t+a] = Add
        '''
        #Comp[t+1-1:t+a] = Add

        Fal[Add] = False
        t += a
        Add = np.concatenate([Nei[key] for key in Add])
        Add = unique_elements(Add, isFalse)
        I = Fal[Add]
        Add = Add[I]
        a = len(Add)
      j += t
      k += 1
      Components[k-1] = Comp[:t]
      CompSize[k-1] = t
      if j < studysize:
        C = []
        while i < numC and len(C) ==0:
          i += 1
          C = CutComps[i-1]
          J = Fal[C]
          C = C[J]
        if i == numC and len(C) == 0:
          j = studysize
          i = numC + 1 
      else:
        i = numC +1
    Components = {key:Components[key] for key in range(k)}
    CompSize = CompSize[:k]
  
  # Determine BaseSize and Cont
  Cont = np.ones(k, dtype='bool')
  BaseSize = np.zeros(k, dtype=np.int32)
  #Bases = {key: [] for key in range(k)}
  Bases = {}
  if k > 1:
    Forb[study] = True
    Fal[study] = False
    Fal[Study[0]] = True
    for i in range(k):
      # Deteermine the size of the base of the components
      Set = unique_elements(np.concatenate([Components[i], Study[0]]),isFalse)
      isFalse[Components[i]] = True
      I = isFalse[Set] & Fal[Set]
      isFalse[Components[i]] = False
      Set = Set[I]
      Bases[i] = Set
      BaseSize[i] = len(Set)
    Fal[Study[0]] = False
    Fal[Study[numS-1]] = True
    Forb[study] = True
    for i in range(k):
      # Determine if the component ccan be extended
      Set = unique_elements(np.concatenate([Components[i], Study[numS-1]]), isFalse)
      isFalse[Components[i]] = True
      I = isFalse[Set] & Fal[Set]
      isFalse[Components[i]] = False
      Set = Set[I]
      if len(Set) > 0:
        N = np.concatenate([Nei[key] for key in Set])
        N = unique_elements(N,isFalse)
        I = Forb[N]
        N = N[~I]
        if len(N) == 0:
          Cont[i] = False
      else:
        Cont[i] = False
  return Components, Bases, CompSize, Cont, BaseSize
      

