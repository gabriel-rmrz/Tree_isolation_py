import numpy as np

def verticalcat(Dict):
  # Vertical concatenation of the given dictionary into a vector

  keys = sorted(Dict.keys())
  DictSize = [len(Dict[k]) for k in keys] # Determine the size of every dictionary
  numC = len(Dict) # Number of elements
  IndElements = np.zeros((numC, 2), dtype=int) # indices for elements in each cell
  IndElements[:,1] = np.cumsum(DictSize)
  IndElements[1:,0] = IndElements[1:,0] + IndElements[:-1, 1]
  #IndElements[:,0] = IndElements[:,0]
  Vector = np.zeros(np.sum(DictSize), dtype=int)
  for j in range(numC):
    Vector[IndElements[j,0]:IndElements[j,1]] = Dict[keys[j]]

  return Vector, IndElements
