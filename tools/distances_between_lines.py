import numpy as np
from tools.mat_vec_subtraction import mat_vec_subtraction


def distances_between_lines(PointRay,DirRay,PointLines,DirLines):
  '''
  
  Calculates the distances between a ray and lines
  
  PointRay      A point of the ray
  DirRay        Unit direction vector of the line
  PointLines    One point of every line
  DirLines      Unit direction vectors of the lines 
  
  '''
  PointLines = PointLines.astype(float)
  PointRay = PointRay.astype(float)
  DirLines = DirLines.astype(float)
  DirRay = DirRay.astype(float)

  # Calculate unit vectors N orthogonal to the ray and the lines
  N = np.column_stack((DirRay[1]*DirLines[:,2]-DirRay[2]*DirLines[:,1], DirRay[2]*DirLines[:,0]-DirRay[0]*DirLines[:,2], DirRay[0]*DirLines[:,1]-DirRay[1]*DirLines[:,0]))
  l = np.sqrt(np.sum(N*N,1))
  N = np.column_stack((1/l*N[:,0], 1/l*N[:,1], 1/l*N[:,2]))

  # Calculate the distances between the lines
  A = -mat_vec_subtraction(PointLines, PointRay)
  DistLines = np.sqrt(np.abs(np.sum(A*N,1))) # distances between the lines and the ray

  # Calculate the distances on ray and on lines
  b = DirLines@np.transpose(DirRay)
  d = A@np.transpose(DirRay)
  e = np.sum(A*DirLines,1)
  DistOnRay = (b*e-d)/(1-b**2) # Distances to PointRay from the closest point on the ray
  DistOnLines = (e-b*d)/(1-b**2) # Distances to PoinLines from the closest points on the lines

  return DistLines, DistOnRay, DistOnLines
