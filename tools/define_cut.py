DEBUG = True
import numpy as np
from tools.unique_elements import unique_elements
def define_cut(Nei, CutPre, Forb, Fal):
  # Defines the "Cut" region
  Cut = np.concatenate([Nei[cut] for cut in CutPre])
  if DEBUG:
    print(f"len(Fal): {len(Fal)}")
    print(f"type(Fal): {type(Fal)}")
    print(f"type(Fal[0]): {type(Fal[0])}")
    print(f"len(Cut): {len(Cut)}")
    print(f"type(Cut): {type(Cut)}")
    print(f"type(Cut[0]): {type(Cut[0])}")
    print(f"len(np.unique(Cut)): {len(np.unique(Cut))}")
  Cut = unique_elements(Cut,Fal)
  #Cut = np.unique(Cut)
  #Cut = Cut[Fal[Cut]]
  if DEBUG:
    print(f"len(Cut) afte unique_elements: {len(Cut)}")
  I = Forb[Cut]
  Cut = Cut[~I]
  if DEBUG:
    print(f"len(Cut) After forb: {len(Cut)}")

  return Cut
