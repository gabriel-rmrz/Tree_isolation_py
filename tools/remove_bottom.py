DEBUG=False
import numpy as np
def remove_bottom(P, Points, Hei, inputs):
  if DEBUG:
    import scipy.io
    print(f"type(P): {type(P)}")
    print(f"(P.shape): {(P.shape)}")
    print(f"type(Hei): {type(Hei)}")
    print(f"(Hei.shape): {(Hei.shape)}")
    print(f"type(Points): {type(Points)}")
    print(f"(Points.shape): {(Points.shape)}")

    #P = np.transpose(np.asarray(scipy.io.loadmat('debug/remove_bottom/P.mat')['P']))
    print("##############################################")
    print("##############################################")
    print(f"READING and setting pointcloud from MATLAB")
    print("##############################################")
    print("##############################################")

    #P = np.asarray(scipy.io.loadmat('debug/remove_bottom/P.mat')['P'])
    #Hei = np.asarray(scipy.io.loadmat('debug/remove_bottom/Hei.mat')['Hei'])
    #Hei = Hei[:,0].astype(int)
    #Points = np.asarray(scipy.io.loadmat('debug/remove_bottom/Points.mat')['Points'])
    #Points = Points[:,0].astype(bool)
    print(f"(P.shape) (technically P_mat): {(P.shape)}")
    print(f"(Hei.shape) (technically Hei_mat): {(Hei.shape)}")
    print(f"(Points.shape) (technically Points_mat): {(Points.shape)}")

  print('---------')
  print('Remove the bottom...')

  ## Determine the bottom fill
  # For each (BottomSquare m x BottomSquarem) square, remove all bottom points
  # below BottomHeight (m) if the trhird metter above the ground (2-3 m) is
  # empty of points
  numP = len(P[:,0])
  SQ = inputs['BottomSquare']
  Min = P[:,:2].min(axis=0).astype(float)
  Max = P[:,:2].max(axis=0).astype(float)
  nx = int(np.ceil((Max[0]-Min[0])/SQ))
  ny = int(np.ceil((Max[1]-Min[1])/SQ))
  n=3
  Filled = np.zeros((nx,ny,n), dtype='bool')
  PointInd = np.asarray(range(numP))
  I = Hei < n*100
  
  PointInd = PointInd[I]
  R = np.floor((P[PointInd,:2]-Min)/SQ) + 1
  J1 = R[:,0] > nx
  if DEBUG:
    print(f"Min: {Min}")
    print(f"Max: {Max}")
    print(f"nx: {nx}")
    print(f"ny: {ny}")
    print(f"sum(I): {sum(I)}")
    print(f"sum(J1): {sum(J1)}")
  if J1.any():
    R[J1,0] = nx
  J2 = R[:,1] > ny
  if J2.any():
    R[J2,1] = ny
  J3 = R <= 0
  if J3.any():
    R[J3] = 1 
  L = (np.floor((Hei[PointInd]).astype(float)/100.) + 1).astype(int)
  J4 = L > n
  if J4.any():
    L[J4] = n 
  J5 = (L <= 0)
  if DEBUG:
    print(f"sum(J5): {sum(J5)}")
    print(f"J5.any(): {J5.any()}")
    print(f"L[L<=0]: {L[L<=0]}")
    print(f"L[J5]: {L[J5]}")
    print(f"R: {R}")
    print(f"len(R): {len(R)}")
    print(f"Hei[PointInd][L<=0]: {Hei[PointInd][L<=0]}")
    print(f"Hei[Hei<0]: {Hei[Hei<0]}")
  if J5.any():
    L[J5] = 1
  LexOrd = (R[:,0] + (R[:,1]-1) *nx +(L-1)*nx*ny).astype(int)
  Filled[np.unravel_index(LexOrd-1,Filled.shape,'F')] = True

  ## Remove bottom points
  Pass = np.ones(numP, dtype='bool')
  PointInd = np.asarray(range(numP))
  p = 1

  I = Hei < 100* inputs['BottomHeight']
  PointInd = PointInd[I]
  R = (np.floor((P[PointInd,:2]-Min)/SQ) + 1).astype(int)
  J = R[:,0] > nx
  if J.any():
    J[J,0] = nx
  J = R[:,1] > ny
  if J.any():
    J[J,1] = ny
  J = R <= 0
  if J.any():
    R[J] = 1
  
  for j in range(1,len(PointInd)+1):
    if not Filled[R[j-1,0]-1, R[j-1,1]-1,n-1]:
      Pass[PointInd[j-1]] = False
  p = p+n
  
  # Define the outputs and summarize the results
  numP0 = len(Points)
  Ind = np.asarray(range(numP0))
  print(f"Points: {Points}")
  Ind = Ind[Points]
  print(f"Ind: {Ind}")
  print(f"len(Ind): {len(Ind)}")
  Points[Ind[~Pass]] = False
  H = np.zeros(numP0).astype(int)
  H[Points] = Hei[Pass]
  P = P[Pass,:]
  numP2 = len(P[:,0])
  print(f"\t Points before: {numP}")
  print(f"\t Filtered points: {numP-numP2}")
  print(f"\t Points left: {numP2}")



  return P, Points, H
  

