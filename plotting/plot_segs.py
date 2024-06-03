import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True


def plot_segs(P, comps, ms, Bal=None, prefix=''):
  fig = plt.figure(figsize=(12,8))
  """
  Plots the point cloud segments given in the dictionary "comps".
  If 4 inputs, dictionary values contain the point indexes. If 5 inputs, dictionary values contain
  the indexes of the cover sets given by "Bal".
  "fig" is the figure number and "ms" is the marker size.
  """
  if not comps:
      return

  col = [
    [0.00, 0.00, 1.00],
    [0.00, 0.50, 0.00],
    [1.00, 0.00, 0.00],
    [0.00, 0.75, 0.75],
    [0.75, 0.00, 0.75],
    [0.75, 0.75, 0.00],
    [0.25, 0.25, 0.25],
    [0.75, 0.25, 0.25],
    [0.95, 0.95, 0.00],
    [0.25, 0.25, 0.75],
    [0.75, 0.75, 0.75],
    [0.00, 1.00, 0.00],
    [0.76, 0.57, 0.17],
    [0.54, 0.63, 0.22],
    [0.34, 0.57, 0.92],
    [1.00, 0.10, 0.60],
    [0.88, 0.75, 0.73],
    [0.10, 0.49, 0.47],
    [0.66, 0.34, 0.65],
    [0.99, 0.41, 0.23]
  ]

  n = len(comps)
  if n < 21:
    col = col * (n // 20 + 1)
  else:
    col = np.random.rand(n, 3)
    col[:, 2] *= 0.2  # Adjusting the blue channel

  S = comps[0]
  if isinstance(S, dict):
    n = len(comps)
    for i in range(n):
      S = comps[i]
      if len(S)>0:
        S = np.concatenate([S[key] for key in S.keys()])
        comps[i] = S
      else:
        comps[i] = []

  plt.clf()
  #plt.figure(fig)

  ax = fig.add_subplot(projection='3d')



  if Bal is None:
    # Plot the segments without Bal
    for i, (key, C) in enumerate(comps.items()):
      #ax.scatter(P[C, 0], P[C, 1], P[C, 2], color=col[i % len(col)])
      ax.scatter(P[C,0], P[C,1], P[C,2], color=col[i % len(col)], marker='.', s=1)
  else:
    # Plot the segments with Bal
    num_points = P.shape[0]
    D = np.zeros(num_points, dtype='bool')
    for i, (key, comp) in enumerate(comps.items()):
      if np.any(comp):
        C = (np.unique(np.concatenate([Bal[c] for c in comp if not isinstance(Bal[c],np.int32) and Bal[c].size >0] ))).astype(np.int32)
        I = D[C]
        C = C[~I]
        D[C] = True
        
        #plt.scatter(P[C,0], P[C,1], P[C,2], color=col[i % len(col)], marker='.', s=1)
        ax.scatter(P[C,0], P[C,1], P[C,2], color=col[i % len(col)], marker='.', s=1)
        #ax.scatter(P[C,0], P[C,1], P[C,2], color='b', marker='.', s=1)
  plt.axis('equal')
  plt.savefig(f"plots/{prefix}_segs.png", dpi=500)
  #plt.show()
