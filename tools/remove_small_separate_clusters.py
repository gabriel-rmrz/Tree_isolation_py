import numpy as np
from collections import defaultdict
from cover_sets_plot import cover_sets_plot
from tools.connected_components import connected_components
def remove_small_separate_clusters(P, Points, inputs):
  print('--------')
  print('Remove small separate clusters...')

  # Generate a cover
  inputs['PatchDiam'] = inputs['filt']['PatchDiam']
  inputs['BallRad'] = inputs['filt']['BallRad']
  inputs['nmin'] = inputs['filt']['nmin']
  inputs['NCubes'] = inputs['filt']['NCubes']

  cover = cover_sets_plot(P,inputs)
  
  # Determine the components
  Comps, CompSize = connected_components(cover['neighbor'], [], inputs['FiltClusterSize'])
  Comps_temp=[]
  for key in Comps.keys():
    print(Comps[key])
    Comps_temp = np.concatenate((Comps_temp,Comps[key]),axis=0) 
  Comps = np.asarray(Comps_temp)
  numB = len(cover['ball'])

  Rem = np.ones(numB, dtype='bool')
  if Comps:
     Rem[Comps] = False
  K = np.asarray(list(cover['ball'].keys()))
  K = K[Rem]
  Remove_temp = np.asarray([])
  for key in K:
    Ball = cover['ball'][key]
    Remove_temp = np.concatenate((Ball, Remove_temp),axis=0)
  Remove = np.asarray(Remove_temp).astype(np.uint32)
  numP = len(Points)
  Ind = np.asarray(range(numP)).astype(np.uint32)
  Ind = Ind[Points]
  Points[Ind[Remove]] = False

  # Display filtering result
  print(f"\t Points before: {len(P[:,0])}")
  print(f"\t Filtered points: {sum(Remove)}")
  print(f"\t Points left: {sum(Points)}")
  




  return Points

  
