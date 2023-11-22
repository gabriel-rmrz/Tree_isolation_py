
def shortest_paths(cover, Base, Forb):
  '''
  Computes the shortest path distances from every cover set to the base using
  Dijkstra's algorithm. Computes also the number of these paths going through 
  every set.
  '''
  numB = len(cover['ball'])
  NeiDis = cover['NeiDis']

  ## Determine the shortest paths
  Unvisited = np.zeros(numB, dtype=int)
  n = len(Base)
  Unvisited[:n] = Base
  UV = np.zeros(numB, dtype='bool')
  UV[Base] = True
  a = 1
  b = n
  J = 0
  PathDisT = 1000*np.ones(numB, dtype=int)
  PathNei = np.zeros(numB, dtype=int)
  PathDist[Base] = 0
  C = Base[0]
  EndSet = np.zeros(numB, dtype=int)
  EndSet[Base] = Base
  while a <= b:
    N = cover['neighbor'][C]
    d = NeiDist[C]
    D = PathDist[C-1] + d
    I = (D < PathDist[N]) and ~Forb[N]
    N = N[I]
    PathDist[N] = D[I]
    PathNei[N] = C
    EndSet[N] = EndSet[C]
    if J > 0:
      Unvisited[a-1+J-1] = Unvisited[a-1]
    a = a+1
    UV[C] = False
    I = UV[N]
    N = N[~I]
    Invisited[b:b+len(N)+1] = N
    b = b + len(N)
    UV[N] = True
    J = np.argmin(PathDist[Unvisited[a:b+1]])
    C = Unvisited[a -1+J-1]
  PathDist[PathDist ==1000] = np.max(Path[PathDist< 1000])+0.1

  ## Compute the number of shortest paths going through each set
  PathNum = np.ones(numB, dtype=int)
  for i in range(1, numB+1):
    N = PathNei[i-1]
    while N > 0:
      PathNum[N] = PathNum[N] + 1
      N = PathNei[N]


  return  PathNum, PathDist, EndSet
