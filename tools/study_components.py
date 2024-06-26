DEBUG = True
import numpy as np
from tools.unique_elements import unique_elements

def study_components(Nei, numS, Cut, CutComps, Forb, Fal, isFalse):
  Study = {}
  StudySize = np.zeros(numS,dtype=np.int32)
  Study[0] = Cut
  StudySize[0] = len(Cut)
  if DEBUG:
    print(f"StudySize: { StudySize}")
    print(f"Study: { Study}")
    print(f"numS: {numS}")

  if numS >= 2:
    N = np.copy(Cut)
    i = 1
    while i < numS:
      Forb[N] = True
      N = np.concatenate([Nei[key] for key in N])
      N = unique_elements(N,Fal)
      '''
      if DEBUG:
        print(f"N: {N}")
        print(f"len(N): {len(N)}")
      '''
      I = Forb[N]
      N = N[~I]
      '''
      if DEBUG:
        print(f"N: {N}")
        print(f"len(N): {len(N)}")
      '''
      if len(N) > 0:
        i +=1
        Study[i-1] = N
        StudySize[i-1] = len(N)
      else:
        Study = {key: Study[key] for key in range(i)} 
        StudySize = [StudySize[key] for key in range(i)]
        i = numS +1
  
  '''
  if DEBUG:
    print(f"StudySize: { StudySize}")
    print(f"Study: { Study}")
    print(f"len(Study): { len(Study)}")
    print(f"Study[0]: { Study[0]}")
    print(f"Study[1]: { Study[1]}")
  '''
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
      if isinstance(C, np.ndarray):
        a = len(C)
      elif isinstance(C, (np.int32, np.uint32)):
        a = 1 
      Comp[:a] = C
      Fal[C] = False
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
      if DEBUG:
        print(f"Add: {Add}")
        print(f"isFalse: {isFalse}")
        print(f"len(Add): {len(Add)}")
        print(f"len(isFalse): {len(isFalse)}")
        print(f"sum(isFalse): {sum(isFalse)}")
        print(f"C: {C}")
        for key in C:
          print(f"len(Nei[{key}]): {len(Nei[key])}")
          print(f"len(np.concatenate([Nei[key] for key in C])): {len(np.concatenate([Nei[key] for key in C]))}")
        exit()
      t =a
      I = Fal[Add]
      Add = Add[I]
      a = len(Add)
      while a>0:
        '''
        if len(Comp[t+1-1:t+a]) < len(Add):
          diff = len(Add) - len(Comp[t+1-1:t+a])
          Comp[t+1-1:t+a] = Add[:a-diff]
          for d in range(1, diff +1):
            np.append(Comp, Add[a-d])
        else:
          Comp[t+1-1:t+a] = Add
        '''
        Comp[t+1-1:t+a] = Add

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
  
  if DEBUG:
    print(f"Components: {Components}")
    print(f"len(Components): {len(Components)}")
    print(f"len(Components[0]): {len(Components[0])}")
    exit()
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
      

