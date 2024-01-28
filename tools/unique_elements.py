import numpy as np
def unique_elements(Set, isFal):
  isFalse = np.copy(isFal)
  print(f"sum(isFalse): {sum(isFalse)}")
  print(f"len(isFalse): {len(isFalse)}")
  n = len(Set)
  if n > 2:
    I = np.ones(n, dtype='bool')
    for j in range(n):
      '''
      print(f"j: {j}")
      print(f"Set[j]: {Set[j]}")
      print(f"isFalse[Set[j]]: {isFalse[Set[j]]}")
      '''
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

            
