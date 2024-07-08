DEBUG= False
import numpy as np
def remove_bottom(P, Points, Hei, inputs):

  print('---------')
  print('Remove the bottom...')

  ## Determine the bottom fill
  # For each (BottomSquare m x BottomSquarem) square, remove all bottom points
  # below BottomHeight (m) if the third meter above the ground (2-3 m) is
  # empty of points
  numP = len(P[:,0])
  SQ = inputs['BottomSquare']
  Min = P[:,:2].min(axis=0).astype(np.double)
  Max = P[:,:2].max(axis=0).astype(np.double)
  nx = int(np.ceil((Max[0]-Min[0])/SQ))
  ny = int(np.ceil((Max[1]-Min[1])/SQ))
  n=3
  Filled = np.zeros((nx,ny,n), dtype='bool')
  PointInd = np.asarray(range(numP), dtype=np.uint32)
  I = Hei < n*100

  PointInd = PointInd[I]
  R = np.floor((P[PointInd,:2]-Min)/SQ) + 1
  J = R[:,0] > nx
  if np.any(J):
    R[J,0] = nx
  J = R[:,1] > ny
  if np.any(J):
    R[J,1] = ny
  J = R <= 0
  if np.any(J):
    R[J] = 1 
  L = (np.floor((Hei[PointInd])/100.) + 1).astype(np.uint32)
  J = L > n
  if np.any(J):
    L[J] = n 
  J = (L <= 0)
  if np.any(J):
    L[J] = 1
  LexOrd = (R[:,0] + (R[:,1]-1) *nx +(L-1)*nx*ny).astype(np.uint32)
  Filled[np.unravel_index(LexOrd-1,Filled.shape,'F')] = True


  ## Remove bottom points
  Pass = np.ones(numP, dtype='bool')
  PointInd = np.asarray(range(numP), dtype=np.uint32)
  p = 1

  I = Hei < 100* inputs['BottomHeight']
  PointInd = PointInd[I]
  R = (np.floor((P[PointInd,:2]-Min)/SQ) + 1).astype(np.int32)

  J = R[:,0] > nx
  if np.any(J):
    R[J,0] = nx
  J = R[:,1] > ny
  if np.any(J):
    R[J,1] = ny
  J = R <= 0
  if np.any(J):
    R[J] = 1
  
  for j in range(len(PointInd)):
    if not Filled[R[j,0]-1, R[j,1]-1,n-1]:
      Pass[PointInd[j]] = False
  p = p+n
  
  # Define the outputs and summarize the results
  numP0 = len(Points)
  Ind = np.asarray(range(numP0))
  Ind = Ind[Points]
  Points[Ind[~Pass]] = False
  H = np.zeros(numP0).astype(np.int32)
  H[Points] = Hei[Pass]
  P = P[Pass,:]
  numP2 = len(P[:,0])
  print(f"\t Points before: {numP}")
  print(f"\t Filtered points: {numP-numP2}")
  print(f"\t Points left: {numP2}")

  return P, Points, H
  

