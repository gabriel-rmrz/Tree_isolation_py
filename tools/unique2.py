DEBUG=False
import numpy as np
def unique2(Set):
  n = len(Set)
  if n>0:
    Set= np.sort(Set)
    d = Set[1:]-Set[:n-1]
    if DEBUG:
      print(f"len(d): {len(d)}")
      print(f"len(Set): {len(Set)}")
      print(f"len(Set[1:]): {len(Set[1:])}")
      print(f"len(Set[:n-1]): {len(Set[:n-1])}")
      exit()
    A = Set[1:]
    I = d>0
    return np.concatenate((np.array([Set[0]]),A[I]), axis=0)
  else:
    return Set


