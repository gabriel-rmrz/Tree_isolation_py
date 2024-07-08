DEBUG = False 
import numpy as np
from tools.unique_elements import unique_elements
def define_cut(Nei, CutPre, Forb, Fal):
  CutPre = np.copy(CutPre)
  Forb = np.copy(Forb)
  Fal = np.copy(Fal)
  # Defines the "Cut" region
  if CutPre.size > 1:
    Cut = np.concatenate([Nei[cut] for cut in CutPre])
  elif len(CutPre) ==0:
    return []
  else:
    Cut = np.copy(Nei[CutPre.item()])
  Cut = unique_elements(Cut,Fal)
  I = Forb[Cut]
  Cut = Cut[~I]
  return Cut
