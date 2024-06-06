DEBUG=False


import scipy.io
import numpy as np
from collections import defaultdict
from scipy.spatial import Delaunay

def compute_height(P, inputs):
  print('---------')
  print('Compute the height of the points...')
  ## Define the ground as points from each nonempty sq rectangle
  # Basic parameters
  numP = len(P[:,0])
  Min = P.min(axis=0).astype(np.double)
  Max = P.max(axis=0).astype(np.double)

  sq = inputs["HeightSquare"] # edge length of squares used in bottom definition

  # Define grid dimensions for the ground level

  # Number of rectangles in x,y directions
  #N = (Max-Min)/sq 
  Nx = (Max[0] - Min[0])/sq
  Ny = (Max[1] - Min[1])/sq
  

  if np.floor(Nx) == np.ceil(Nx):
    Nx = Nx + 1
  else:
    Nx = np.ceil(Nx)
  if np.floor(Ny) == np.ceil(Ny):
    Ny = Ny + 1
  else:
    Ny = np.ceil(Ny)

  Nx = int(Nx + 2)
  Ny = int(Ny + 2)

  #N = N.astype(np.int32)
  Bot = np.zeros((Nx, Ny)) + Min[2] + 100
  if DEBUG:
    print(f"Bot: {Bot}")

  #BL = defaultdict(list)
  BL = {}
  PointInd = np.asarray(range(numP), dtype=np.uint32)

  # Define Bot-values

  
  R = np.floor((P[:,:2]-Min[:2])/sq) + 2.
  LexOrd = (R[:,0] + (R[:,1] -1)*Nx).astype(np.int32)
  I = np.argsort(LexOrd)
  LexOrd = np.sort(LexOrd)
  PointInd = PointInd[I]


  q = 1
  m = numP
  while q <= m:
    t = 1
    while q+t <= m and LexOrd[q-1+t] == LexOrd[q-1]:
      t = t+1
    a = (P[PointInd[q-1: q+t-1],2]).min().astype(np.double)
    if Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')] > a:
      Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')] = a
    # Select the lowest 30 cm layer of points
    I = PointInd[q-1: q+t-1]
    J = P[PointInd[q-1:q+t-1],2] < (Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')]+0.3)
    BL[LexOrd[q-1]-1] = I[J]
    q = q + t
  if DEBUG:
    print(f"Bot: {Bot}")

  # Define Ground, a grid of ground level point as the mean of the 10%
  # lowest points of the lowest 30 cm
  t = 0
  Ground = np.zeros([Nx*Ny,3])
  if DEBUG:
    print(f"Ground: {Ground}")

  for i in range(2,Nx):
    for j in range(2,Ny):
      k=(j-1)*Nx+i
      if (np.asarray(BL[k-1])).any():
        I = np.argsort(P[BL[k-1],2])
        a = np.ceil(len(I)*0.1).astype(np.int32)
        Q = np.mean(P[BL[k-1][I[:a]],:],0)
        t = t+1
        Ground[t-1,:] = Q

        if i == 2:
          t = t+1
          Ground[t-1,:] = Q-sq*np.asarray([1, 0, 0])
        elif i == (Nx -1):
          t = t+1
          Ground[t-1,:] = Q+sq*np.asarray([1, 0, 0])
        if j == 2:
          t = t+1
          Ground[t-1,:] = Q-sq*np.asarray([0, 1, 0])
        elif j == (Ny -1):
          t = t+1
          Ground[t-1,:] = Q+sq*np.asarray([0, 1, 0])
  Ground = Ground[:t,:]
  if DEBUG:
    print(f"Ground: {Ground}")

  
  # Triangulate the ground points
  Tri = Delaunay(Ground[:,:2])
  Tri = Tri.simplices

  ## Generate more ground points from each triangle
  n = len(Tri[:,0])
  SQ = sq/40
  SQ = max(SQ, 0.15)
  a = np.ceil(sq/SQ*5)
  b = int(min(a**2*n,2e9))
  G = np.zeros([b,3])

  t = 0
  for i in range(1, n+1):
    T = (Tri[i-1,:]).astype(np.int32)
    Q = Ground[T,:]
    V1 = Q[1,:]- Q[0,:]
    V2 = Q[2,:]- Q[0,:]
    n1 = (np.ceil(np.linalg.norm(V1,1)/SQ*2)).astype(np.int32)
    n2 = (np.ceil(np.linalg.norm(V2,1)/SQ*2)).astype(np.int32)
    n1 = max(n1,n2)
    t = t+1
    G[t-1,:] = Q[0,:]
    for j in range(1, n1+1):
      V = Q[0,:] + j/n1*V2-Q[0,:]-j/n1*V1
      n2 = int(np.ceil(np.linalg.norm(V)/SQ*2))
      for k in range(n2+1):
        t = t+1
        G[t-1,:] = Q[0,:] + j/n1*V1+k/n2*V
    if t > 0.75*b:
      b = np.ceil(n/i)*b
      G[b,2] = 0
  G = G[:t,:]

  # Define height of the points
  Min = Ground[:,:2].min(axis=0).astype(np.double)
  Max = Ground[:,:2].max(axis=0).astype(np.double)
  nx = (Max[0]-Min[0])/SQ
  ny = (Max[1]-Min[1])/SQ

  if np.floor(nx) == np.ceil(nx):
    nx = nx + 1
  else:
    nx = np.ceil(nx)
  if np.floor(ny) == np.ceil(ny):
    ny = ny + 1
  else:
    ny = np.ceil(ny)
  
  nx = int(nx)
  ny = int(ny)

  
  Bot = np.zeros([nx,ny])
  R = np.floor((G[:,:2]-Min)/SQ)+1
  LexOrd = (R[:,0] + (R[:,1] -1) * nx).astype(np.uint32)
  LexOrd, I = np.unique(LexOrd, return_index=True, axis=0)
  G = G[I,:]
  Bot[np.unravel_index(LexOrd-1, Bot.shape,'F')] = G[:,2]


  Hei = np.zeros(numP, dtype=np.int32)
  p = 1
  PointInd = np.asarray(range(numP), dtype=np.uint32)
  R = np.floor((P[:,:2]-Min[:2])/SQ) + 1
  LexOrd = (R[:,0] + (R[:,1] -1) * nx).astype(np.uint32)
  SortOrd = np.argsort(LexOrd)
  LexOrd = np.sort(LexOrd)
  PointInd = PointInd[SortOrd]
  R = R[SortOrd,:]
  q = 1
  while q <= m:
    t =1
    while q+t <= m and LexOrd[q-1+t] == LexOrd[q-1]:
      t = t+1
    if Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape, 'F')] == 0:
      
      min1 = int(max(1,R[q-1,0]-a))
      max1 = int(min(nx,R[q-1,0]+a))
      min2 = int(max(1,R[q-1,1]-a))
      max2 = int(min(ny,R[q-1,1]+a))
      bot = Bot[min1:(max1+1), min2:(max2+1)]
      while not bot.any():
        a = a+1
        min1 = int(max(1,R[q-1,0]-a))
        max1 = int(min(nx,R[q-1,0]+a))
        min2 = int(max(1,R[q-1,1]-a))
        max2 = int(min(ny,R[q-1,1]+a))
        bot = Bot[min1:(max1+1), min2:(max2+1)]
      bot = np.mean(bot[bot != 0])
      Hei[PointInd[q-1:q+t-1]] = (P[PointInd[q-1:q+t-1],2] - bot)*100
      if DEBUG:
        print('here1')
        exit()
    else:
      Hei[PointInd[q-1:q+t-1]] = (P[PointInd[q-1:q+t-1],2] - Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')])*100
      if DEBUG:
        print('here2')
    if DEBUG:
      print(q)
      print(f"Hei: {Hei}")
      #print(f"(P[PointInd[q-1:q+t-1],2] - bot)*100: {(P[PointInd[q-1:q+t-1],2] - bot)*100}")
    q = q+t
  
  Hei[Hei < 0] = 0
  if DEBUG:
    print(f"len(Ground): {len(Ground)}")
    print(f"len(Tri): {len(Tri)}")
    np.savetxt("Hei.txt", Hei, fmt='%d')
    np.savetxt("SortOrd.txt", SortOrd, fmt='%d')
    np.savetxt("Bot.txt", Bot, fmt='%d')
    print(f"len(Hei): {(Hei)}")
    exit()

  return Hei, Ground, Tri








