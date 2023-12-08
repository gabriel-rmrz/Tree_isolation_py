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
  if DEBUG:
    P_mat = np.asarray(scipy.io.loadmat('debug/compute_height/p.mat')['P'])
    import matplotlib.pyplot as plt
    plt.clf()
    plt.hist(P[:,0], bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/P_x.png")
    plt.clf()
    plt.hist(P_mat[:,0], bins=80,alpha=0.5)
    plt.savefig(f"{dir_plots}/P_x_mat.png")
    plt.clf()
    plt.hist(P[:,1], bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/P_y.png")
    plt.clf()
    plt.hist(P_mat[:,1], bins=80,alpha=0.5)
    plt.savefig(f"{dir_plots}/P_y_mat.png")
    plt.clf()
    plt.hist(P[:,2], bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/P_z.png")
    plt.clf()
    plt.hist(P_mat[:,2], bins=80,alpha=0.5)
    plt.savefig(f"{dir_plots}/P_z_mat.png")
    #P = P_mat

  ## Define the ground as points from each nonempty sq rectangle
  # Basic parameters
  numP = len(P[:,0])
  Min = P.min(axis=0).astype(float)
  Max = P.max(axis=0).astype(float)

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
  N = N.astype(int)

  Bot = np.zeros(N[:2]) + Min[2] + 100

  BL = defaultdict(list)
  PointInd = np.asarray(range(numP))

  # Define Bot-values

  
  R = np.floor((P[:,:2]-Min[:2])/sq) + 2
  LexOrd = (R[:,0] + (R[:,1] -1)*N[0]).astype(int)
  I = np.argsort(LexOrd)
  LexOrd = np.sort(LexOrd)
  PointInd = PointInd[I]

  if DEBUG:
    print(f"Bot:")
    print2darray(Bot)
    print(f"N: {N}")

  q = 1
  m = numP
  while q <= m:
    t = 1
    while q+t <= m and LexOrd[q-1+t] == LexOrd[q-1]:
      t = t+1
    a = (P[PointInd[q-1: q+t-1],2]).min().astype(float)
       
    if Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')] > a:
      Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')] = a
    '''
    if DEBUG:
      print(f"q: {q}")
      print(f"LexOrd[q-1]: {LexOrd[q-1]}")
      print(f"np.unravel_index(LexOrd[q-1],Bot.shape,'F'): {np.unravel_index(LexOrd[q-1],Bot.shape,'F')}")
      print(f"np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F'): {np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')}")
      #print(f"np.unravel_index(LexOrd[q-1],Bot.shape,'F')-(1,1): {np.unravel_index(LexOrd[q-1],Bot.shape,'F')-(1,1)}")
    '''

    # Select the lowest 30 cm layer of points
    I = PointInd[q-1: q+t-1]
    J = P[PointInd[q-1:q+t-1],2] < (Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')]+0.3)
    K = BL[LexOrd[q-1]]
    J2 = P[K,2] < (Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')] + 0.3)
    if J2:
      BL[LexOrd[q-1]] = np.concatenate(K[J2],I[J])
    else:
      BL[LexOrd[q-1]] = I[J]
    if DEBUG:
      print(f"K: {K}")
      print(f"I: {I}")
      print(f"J: {J}")
      print(f"J2: {J2}")
      print(f"a: {a}")
      print(f"len(I): {len(I)}")
      print(f"len(J): {len(J)}")
      print(f"I[J]: {I[J]}")
      print(f"len(I[J]): {len(I[J])}")
    q = q + t
  if DEBUG:
    print(f"Bot after while:")
    print2darray(Bot)
    print(f"BL: {BL}")
  # Define Ground, a grid of ground level point as the mean of the 10%
  # lowest points of the lowest 30 cm
  t = 0
  Ground = np.zeros([N[0]*N[1],3])

  for i in range(2,N[0]):
    for j in range(2,N[1]):
      k=(j-1)*N[0]+i
      '''
      if DEBUG:
        k2 = (np.ravel_multi_index([[i-1],[j-1]], Bot.shape,'raise','F'))[0] + 1
        print(f"k: {k}")
        print(f"k2: {k2}")
      '''
                                  
      if (np.asarray(BL[k])).any():
        I = np.argsort(P[BL[k],2])
        a = np.ceil(len(I)*0.1).astype(int)
        Q = np.mean(P[BL[k][I[:a]],:],0)
        t = t+1
        Ground[t-1,:] = Q
        '''
        if DEBUG:
          print(f"Ground[t-1]: {Ground[t-1]}")
          print(f"a: {a}")
          print(f"Q: {Q}")
          print(f"I[:a]: {I[:a]}")
          print(f"BL[k][I[:a]]: {BL[k][I[:a]]}")
        '''


        if i == 2:
          t = t+1
          Ground[t-1,:] = Q-sq*np.asarray([1, 0, 0])
          '''
          if DEBUG:
            print(f"Ground: {Ground}")
            print(f"Q: {Q}")
            print(f"Ground[t-1,:]: {Ground[t-1,:]}")
          '''
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

  if DEBUG:
    Ground_mat = np.asarray(scipy.io.loadmat('debug/compute_height/ground.mat')['Ground'])
    print(f"Ground:")
    print2darray(Ground)
    print(f"len(Ground): {len(Ground)}")
    print(f"len(BL): {len(BL)}")
    print(f"sq: {sq}")
    print(f"Ground:")
    import matplotlib.pyplot as plt
    plt.clf()
    plt.hist(Ground[:,0], bins=40, alpha=0.5)
    plt.savefig(f"{dir_plots}/Ground_x.png")
    plt.clf()
    plt.hist(Ground_mat[:,0], bins=40, alpha=0.5)
    plt.savefig(f"{dir_plots}/Ground_x_mat.png")
    plt.clf()
    plt.hist(Ground[:,1], bins=40, alpha=0.5)
    plt.savefig(f"{dir_plots}/Ground_y.png")
    plt.clf()
    plt.hist(Ground_mat[:,1], bins=40, alpha=0.5)
    plt.savefig(f"{dir_plots}/Ground_y_mat.png")
    plt.clf()
    plt.hist(Ground[:,2], bins=40, alpha=0.5)
    plt.savefig(f"{dir_plots}/Ground_z.png")
    plt.clf()
    plt.hist(Ground_mat[:,2], bins=40, alpha=0.5)
    plt.savefig(f"{dir_plots}/Ground_z_mat.png")

    print(f"Ground-ground_mat:")
    #print2darray(Ground-Ground_mat)
    '''
    for key in BL.keys():
      print(f"key: {key}")
      #print(f"BL[key]: {BL[key]}")
    '''
    #Ground = Ground_mat
  
  # Triangulate the ground points
  Tri = Delaunay(Ground[:,:2])
  Tri = Tri.simplices
  if DEBUG:
    Tri_mat = Delaunay(Ground_mat[:,:2])
    import matplotlib.pyplot as plt
    plt.clf()
    plt.triplot(Ground[:,0], Ground[:,1],Tri)
    plt.plot(Ground[:,0], Ground[:,1],'o')
    plt.savefig(f"{dir_plots}/Delaunay.png")
    plt.clf()
    plt.triplot(Ground_mat[:,0], Ground_mat[:,1],Tri_mat.simplices)
    plt.plot(Ground_mat[:,0], Ground_mat[:,1],'o')
    plt.savefig(f"{dir_plots}/Delaunay_mat.png")
  

  ## Generate more ground points from each triangle
  n = len(Tri[:,0])
  SQ = sq/40
  SQ = max(SQ, 0.15)
  a = np.ceil(sq/SQ*5)
  b = int(min(a**2*n,2e9))
  G = np.zeros([b,3])

  if DEBUG:
    print(f"Ground: {Ground}")
    print(f"len(Ground): {len(Ground[:,0])}")
    #print(f"Tri: {Tri.simplices}")
    ts = Tri
    tsm = Tri_mat.simplices
    print(f"len(Tri): {ts[ts[:,0] == 143]}")
    print(f"len(Tri_mat.simplices): {tsm[tsm[:,0] == 143]}")
    #print(f"len(Tri): {ts[:,0]}")
    #print(f"len(Tri_mat.simplices): {tsm[:,0]}")
    print(f"a: {a}")
    print(f"SQ: {SQ}")
    print2darray(ts)
  t = 0
  for i in range(1, n+1):
    T = (Tri[i-1,:]).astype(int)
    Q = Ground[T,:]
    if DEBUG:
      print(f"Q: {Q}")
      print(f"Q[0,:]: {Q[0,:]}")
      print(f"T: {T}")
    V1 = Q[1,:]- Q[0,:]
    V2 = Q[2,:]- Q[0,:]
    n1 = (np.ceil(np.linalg.norm(V1,1)/SQ*2)).astype(int)
    n2 = (np.ceil(np.linalg.norm(V2,1)/SQ*2)).astype(int)
    if DEBUG:
      print(f"V1: {V1}")
      print(f"V2: {V2}")
      print(f"n1: {n1}")
      print(f"n2: {n2}")
    n1 = max(n1,n2)
    if DEBUG:
      print(f"n1: {n1}")
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
  if DEBUG:
    G_mat = np.asarray(scipy.io.loadmat('debug/compute_height/g.mat')['G'])
    import matplotlib.pyplot as plt
    plt.clf()
    plt.hist(G[:,2], bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/G_z.png")
    plt.clf()
    plt.hist(G_mat[:,2], bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/G_mat_z.png")
    plt.hist(G[:,2], bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/G_both.png")

    print(f"G: {G} ")
    print(f"len(G): {len(G)} ")


  # Define height of the points
  Min = Ground[:,:2].min(axis=0).astype(float)
  Max = Ground[:,:2].max(axis=0).astype(float)
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
  LexOrd = (R[:,0] + (R[:,1] -1) * nx).astype(int)
  LexOrd, I = np.unique(LexOrd, return_index=True, axis=0)
  G = G[I,:]
  Bot[np.unravel_index(LexOrd-1, Bot.shape,'F')] = G[:,2]

  if DEBUG:
    Bot_mat = np.asarray(scipy.io.loadmat('debug/compute_height/bot.mat')['Bot'])
    print(f"G[:,2]: {G[:,2]}")
    print(f"len(G[:,2]): {len(G[:,2])}")
    print(f"Bot[np.unravel_index(LexOrd-1, Bot.shape,'F')]: {Bot[np.unravel_index(LexOrd-1, Bot.shape,'F')]}")
    np.savetxt('Bot.txt',Bot)
    import matplotlib.pyplot as plt
    Bot2 = Bot
    Bot2[Bot2 == 0.] = Bot2.min()
    dif = (Bot2.max() - Bot2.min())
    plt.clf()
    plt.imshow((255*(Bot2- Bot2.min())/dif), cmap ='Blues', interpolation ='bilinear')
    plt.savefig(f"{dir_plots}/Bot.png")
    Bot_mat[Bot_mat == 0.] = Bot_mat.min()
    dif = (Bot_mat.max() - Bot_mat.min())
    plt.clf()
    plt.imshow((255*(Bot_mat - Bot_mat.min())/dif), cmap ='Blues', interpolation ='bilinear')
    plt.savefig(f"{dir_plots}/Bot_mat.png")

  Hei = np.zeros(numP, dtype=np.int16)
  p = 1
  PointInd = np.asarray(range(numP))
  R = np.floor((P[:,:2]-Min[:2])/SQ) + 1
  LexOrd = (R[:,0] + (R[:,1] -1) * nx).astype(int)
  SortOrd = np.argsort(LexOrd)
  PointInd = PointInd[SortOrd]
  R = R[SortOrd,:]
  q = 1
  if DEBUG:
    print(f"LexOrd[265160]: {LexOrd[265160]}")
    print(Bot[np.unravel_index(LexOrd[265160], Bot.shape,'F')])
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
        if DEBUG:
          print("HERE")
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

      '''
      if DEBUG:
        #print(f"Here")
        hei = Hei[PointInd[q-1:q+t-1]]
        p = P[PointInd[q-1:q+t-1],2]
        if (hei[hei<0]).any():
          print(f"q: {q}")
          print(f"t: {t}")
          print("p =P[PointInd[q-1:q+t-1],2]")
          print("hei = Hei[PointInd[q-1:q+t-1]]")
          print(f"hei[hei<0]: {hei[hei<0]}")
          print(f"p[hei<0]: {p[hei<0]}")
          print(f"len(hei): {len(hei)}")
          print(f"Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')]: {Bot[np.unravel_index(LexOrd[q-1]-1,Bot.shape,'F')]}")
      '''

    q = q+t
  
  Hei[Hei < 0] = 0
  if DEBUG:
    print(f"Hei[Hei > 0].min(): {Hei[Hei > 0].min()}")
    print(f"Hei[Hei > 0].max(): {Hei[Hei > 0].max()}")
    import matplotlib.pyplot as plt
    plt.clf()
    plt.hist(Hei, bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/Hei.png")
    plt.clf()
    Hei_mat = np.asarray(scipy.io.loadmat('debug/compute_height/Hei.mat')['Hei'])
    print(f"Hei_mat[Hei_mat > 0].min(): {Hei_mat[Hei_mat > 0].min()}")
    print(f"Hei_mat[Hei_mat > 0].max(): {Hei_mat[Hei_mat > 0].max()}")
    plt.hist(Hei_mat, bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/Hei_mat.png")
    plt.clf()
    plt.hist(Hei_mat[Hei_mat > 100], bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/Hei_mat.png")
    plt.hist(Hei[Hei > 100], bins=80, alpha=0.5)
    plt.savefig(f"{dir_plots}/Hei_both.png")
    plt.clf()
    print(f"len(Hei): {len(Hei)}")
    print(f"Hei[Hei<0]: {Hei[Hei<0]}")
    print(f"len(Hei[Hei<0]): {len(Hei[Hei<0])}")

  return Hei, Ground, Tri








