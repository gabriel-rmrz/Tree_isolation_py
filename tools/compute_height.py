DEBUG=False

if DEBUG:
  dir_mat_files='debug/compute_height'
  dir_plots='debug/compute_height/plots'
  def print2darray(arr):
    for i in arr:
      for j in i: 
        print(f"{j:.4f}", end=' ')
      print()

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

  N = (Max-Min)/sq # Number of rectangles in x,y directions

  if np.floor(N[0]) == np.ceil(N[0]):
    N[0] = N[0] + 1
  else:
    N[0] = np.ceil(N[0])
  if np.floor(N[1]) == np.ceil(N[1]):
    N[1] = N[1] + 1
  else:
    N[1] = np.ceil(N[1])

  N[0] = N[0] + 2
  N[1] = N[1] + 2

  N = N.astype(np.int32)
  Bot = np.zeros(N[:2]) + Min[2] + 100

  BL = defaultdict(list)
  PointInd = np.asarray(range(numP), dtype=np.uint32)

  # Define Bot-values

  
  R = np.floor((P[:,:2]-Min[:2])/sq) + 2.
  LexOrd = (R[:,0] + (R[:,1] -1)*N[0]).astype(np.int32)
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
    K = BL[LexOrd[q-1]-1]
    J2 = P[K,2] < (Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')] + 0.3)
    if J2:
      BL[LexOrd[q-1]-1] = np.concatenate(K[J2],I[J])
    else:
      BL[LexOrd[q-1]-1] = I[J]
    q = q + t
  # Define Ground, a grid of ground level point as the mean of the 10%
  # lowest points of the lowest 30 cm
  t = 0
  Ground = np.zeros([N[0]*N[1],3])

  for i in range(2,N[0]):
    for j in range(2,N[1]):
      k=(j-1)*N[0]+i
      if (np.asarray(BL[k-1])).any():
        I = np.argsort(P[BL[k-1],2])
        a = np.ceil(len(I)*0.1).astype(np.int32)
        Q = np.mean(P[BL[k-1][I[:a]],:],0)
        t = t+1
        Ground[t-1,:] = Q


        if i == 2:
          t = t+1
          Ground[t-1,:] = Q-sq*np.asarray([1, 0, 0])
        elif i == (N[0] -1):
          t = t+1
          Ground[t-1,:] = Q+sq*np.asarray([1, 0, 0])
        if j == 2:
          t = t+1
          Ground[t-1,:] = Q-sq*np.asarray([0, 1, 0])
        elif j == (N[1] -1):
          t = t+1
          Ground[t-1,:] = Q+sq*np.asarray([0, 1, 0])
  Ground = Ground[:t,:]

  
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
    if DEBUG:
      print(f"n1 after: {n1}")
      print(f"b after: {b}")
      print(f"V: {V}")
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

  if DEBUG:
    print(f"nx: {nx}")
    print(f"ny: {ny}")
    print(f"Min: {Min}")
    print(f"Max: {Max}")
  
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
    else:
      Hei[PointInd[q-1:q+t-1]] = (P[PointInd[q-1:q+t-1],2] - Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')])*100


    q = q+t
  
  Hei[Hei < 0] = 0

  return Hei, Ground, Tri








