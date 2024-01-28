import numpy as np

def neighbor_distances(P, cover):
  '''
  Compute neighbor distances "NeiDis"
  '''
  Ce = P[cover["center"]]
  nb = len(cover["ball"])
  NeiDis = {}
  for key in cover["ball"].keys():
    Bot = cover["neighbor"][key] # the neighbors of the ball "i"
    V = Ce[key] - Ce[Bot]
    V = np.power(V,2)
    NeiDis[key]=  np.sqrt(V.sum(axis=1))
  return NeiDis

