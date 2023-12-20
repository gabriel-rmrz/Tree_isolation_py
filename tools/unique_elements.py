import numpy as np
def unique_elements(Set, isFalse):
  n = len(Set)
  if n > 2:
    I = np.ones(n)
    for j in range(1,n+1):
      if ~isFalse[Set[j-1]]:
        isFalse[Set[j-1]] = True
      else:
        I[j-1] = False
    Set = Set[I]
  elif n ==2:
    if Set[0] == Set[1]:
      Set = Set[0]

  return Set

            
