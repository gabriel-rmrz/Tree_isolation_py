import numpy as np
from tools.unique_elements import unique_elements
def define_cut(Nei, CutPre, Forb, Fal):
  # Defines the "Cut" region
  Cut = np.concatenate([Nei[cut] for cut in CutPre])
  Cut = unique_elements(Cut,Fal)
  I = Forb[Cut]
  Cut = Cut[~I]

  return Cut
