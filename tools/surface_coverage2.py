import numpy as np
from tools.orthonormal_vectors import orthonormal_vectors
def surface_coverage2(Axis, Len, Vec, height, numL, numS):
  '''
  Computes surface coverage (number between 0 and 1) of points on cylinder
  surface defined by "Axis" and "Len". "Vec" are the vectors connecting
  points to the Axis and "height" are the heights of the points from
  the base of the cylinder
  '''

  U, W = orthonormal_vectors(Axis)
  Vec = Vec*np.column_stack((U,W))
  ang = np.arctan2((Vec[:,1], Vec[:,0])) + np.pi
  I = int(np.ceil(height/Len*numL))
  I[I==0] = 1
  I[I>numL] = numL
  J = int(np.ceil(ang/(2*np.pi)*numS))
  J[J==0] = 1
  K = np.column_stack((I, J-1))@np.transpose(((1,numL)))
  SurfCov = length(np.unique(K))/(numL*numS)
  return SurfCov
