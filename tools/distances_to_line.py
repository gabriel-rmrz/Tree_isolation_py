import numpy as np

def distances_to_line(Q, LineDirec, LinePoint):
  '''
  Calculates the distances of the points, given in the rows of the
  matrix Q, to the line defined by one of its point and its direction.
  "LineDirec" must be a unit (1x3)-vector and LinePoint must be a (1x3)-vector.
  
  Last update 8 Oct 2021
  '''

  A = Q-LinePoint
  h = np.dot(A, np.transpose(LineDirec))
  B = np.array([LineDirec * a for a in A])
  V = A-B
  d = np.sqrt(np.sum(V*V,1))
  return d, V, h, B
