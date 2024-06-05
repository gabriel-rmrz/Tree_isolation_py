import numpy as np

def restrict_size(P,inputs):
  numP = len(P[:,0])
  bound = inputs["bound"]
  print(len(bound[0]))
  Points = np.zeros(numP, dtype='bool') # The points passing the filtering
  if len(bound[0]) == 1:
    print("---------")
    print(f"Restrict the point cloud into {2*bound[0][0]}m  diameter plot")
    #Points = np.zeros(numP, dtype='bool') # The points passing the filtering
    # Filter the points in blocks of 1e7 points

    Points = (np.power(P,2)).sum(axis=1) < bound[0][0]**2
  else:
    # TODO: This part is not running in the original code, discuss with Pasi
    Points = (np.power(P,2)).sum(axis=1) < bound[0]
  
  # Display results

  print(f"\t {numP} points originally")
  print(f"\t {Points.sum()} points after restriction")
    
  #Points = np.asarray(range(numP))
  return Points
