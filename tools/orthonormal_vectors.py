import numpy as np

def orthonormal_vectors(U):  
  # Generate vectors V and W that are unit vectors orthogonal to themselves
  # and to the input vector U

  V = np.random.rand(3)
  V = np.cross(V,U)
  while np.linalg.norm(V) ==0:
    V = np.random.rand(3)
    V = np.cross(V,U)
  W = np.cross(V,U)
  W = W/np.linalg.norm(W)
  V = V/np.linalg.norm(V)

  return V, W



