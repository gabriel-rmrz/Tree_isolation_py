import numpy as np

# Returns the rotation matriz for the given axis A and angle (in radians)

def rotation_matrix(A,angle):
  A = A/np.linalg.norm(A)
  R = np.zeros((3,3))
  c = np.cos(angle)
  s = np.sin(angle)
  R[0,:] = [A[0]**2+(1-A[0]**2)*c, A[0]*A[1]*(1-c)-A[2]*s, A[0]*A[2]*(1-c)+A[1]*s]
  R[1,:] = [A[0]*A[1]*(1-c)+A[2]*s, A[1]**2+(1-A[1]**2)*c, A[1]*A[2]*(1-c)-A[0]*s]
  R[2,:] = [A[0]*A[2]*(1-c)-A[1]*s, A[1]*A[2]*(1-c)+A[0]*s, A[2]**2+(1-A[2]**2)*c]
  
  return R
