import numpy as np
def shortest_paths(cover, Base, Forb):
  '''
  Computes the shortest path distances from every cover set to the base using
  Dijkstra's algorithm. Computes also the number of these paths going through 
  every set.
  '''
  numB = len(cover['ball'])
  NeiDis = cover['NeiDis']

  ## Determine the shortest paths
  Unvisited = np.zeros(numB, dtype=np.uint32)
  n = len(Base)
  Unvisited[:n] = Base
  UV = np.zeros(numB, dtype='bool')
  UV[Base] = True
  a = 1
  b = n
  J = 1
  PathDis = 1000*np.ones(numB, dtype=np.single)
  PathNei = np.zeros(numB, dtype=np.uint32)
  PathDis[Base] = 0
  C = Base[0]
  EndSet = np.zeros(numB, dtype=np.uint32)
  EndSet[Base] = Base
  while a <= b:
    N = cover['neighbor'][C+1]
    d = NeiDis[C+1]
    D = PathDis[C] + d
    I = (D < PathDis[N]) & ~Forb[N]
    N = N[I]
    PathDis[N] = D[I]
    PathNei[N] = C
    EndSet[N] = EndSet[C]
    if J > 1:
      Unvisited[a+J-1-1] = Unvisited[a-1]
    a +=1
    UV[C] = False
    I = UV[N]
    N = N[~I]
    Unvisited[b:b+len(N)] = N
    b = b + len(N)
    UV[N] = True
    if a<=b:
      J = np.argmin(PathDis[Unvisited[int(a-1):int(b)]])
      C = Unvisited[a -1+J-1]
    else:
      J = np.array([])
      C = np.array([])
  PathDis[PathDis ==1000] = np.max(PathDis[PathDis< 1000])+0.1

  ## Compute the number of shortest paths going through each set
  PathNum = np.ones(numB, dtype=int)
  for i in range(1, numB+1):
    N = PathNei[i-1]
    while N > 0:
      PathNum[N] = PathNum[N] + 1
      N = PathNei[N]


  return  PathNum, PathDis, EndSet
