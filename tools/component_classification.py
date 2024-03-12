DEBUG=False
import numpy as np

def component_classification(CompSize, Cont, BaseSize, CutSize):
  # Classifies study region components:
  # Class[i] == 0 continuation
  # Class[i] == 1 branch
  # Class[i] == 2 part of the continuation (small uncertain thing)

  numC = len(CompSize)
  StudySize = sum(CompSize)
  Class = np.ones(numC, dtype=np.int32) # true if a component is a branch to be further segmented
  ContiComp = 0
  # Simple initial classification
  for i in range(numC):
    if BaseSize[i] == CompSize[i] and ~Cont[i]:
      # component has no expansion, not a branch
      Class[i] = 2
    elif BaseSize[i] == 1 and CompSize[i] <= 2 and ~Cont[i]:
      # component has very small expansion, not a branch
      Class[i] = 2
    elif BaseSize[i]/CutSize < 0.05 and 2*BaseSize[i] >= CompSize[i] and ~Cont[i]:
      # component has very small expansion or is very small, not a branch
      Class[i] = 2
    elif CompSize[i] <= 3 and ~ Cont[i]:
      # very small component, not a branch
      Class[i] = 2
    elif BaseSize[i]/CutSize >= 0.7 or CompSize[i] >= 0.7*StudySize:
      # continuation of the segment
      Class[i] = 0
      ContiComp = i
    #else:
      # Component is probably a branch
  Branches = ( Class == 1)
  if ContiComp == 1 and np.any(Branches):
    Ind = np.arange(numC,dtype=np.uint32)
    if DEBUG:
      print(f"Branches: {Branches}")
      print(f"Ind[Branches]: {Ind[Branches]}")
    Branches = Ind[Branches]
    I = np.argmax(CompSize[Branches])
    Class[Branches[I]] = 0
  
  return Class

