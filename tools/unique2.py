DEBUG=False
import numpy as np
def unique2(Set):
  n = len(Set)
  if n>0:
    Set= np.sort(Set)
    d = Set[1:]-Set[:n-1]
    A = Set[1:]
    I = d>0
    if DEBUG:
      print(f"np.array([Set[0]]): {np.array([Set[0]])}")
      print(f"Set[0]: {Set[0]}")
      print(f"A[I]:  {A[I]}")
    return np.concatenate((np.array([Set[0]]),A[I]), axis=0)
  else:
    return Set


