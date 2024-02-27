import numpy as np
def unique_elements(Set, isFal):
  isFalse = np.copy(isFal)
  n = len(Set)
  if n > 2:
    I = np.ones(n, dtype='bool')
    for j in range(n):
      if ~isFalse[Set[j]]:
        isFalse[Set[j]] = True
      else:
        I[j] = False
    Set = Set[I]
  elif n ==2:
    if Set[0] == Set[1]:
      Set = Set[0]
  del isFalse
  return Set

            
