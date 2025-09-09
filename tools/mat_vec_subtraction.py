import numpy as np

def mat_vec_subtraction(A, v):
  # Subtracts from each row of the matrix A the vector v
  # if A is (n x m)- matrix, then v needs to be m-vector


  for i in range(np.size(A,1)):
    A[:, i] = A[:, i] - v[i]
  return A
