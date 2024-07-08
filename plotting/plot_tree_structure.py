import numpy as np
import matplotlib.pyplot as plt

def set_diff(a,b):
  c = []
  for i in a:
    if i not in b:
      c.append(i)
  return c

def plot_tree_structure(P, cover, segment, ms, B0, segind, prefix):
  fig = plt.figure(figsize=(12,8))
  """
  Plots the branch-segmented tree point cloud so that each branching order
  has its own color Blue = trunk, green = 1st-order branches,
  red = 2nd-order branches, etc.
  
  Inputs
  P         Point cloud
  cover     Cover sets structure
  Segs      Segments structure
  fig       Figure number
  ms        Marker size
  BO        How many branching orders are plotted. 0 = all orders
  segind    Index of the segment where the plotting of tree structure
                    starts. If segnum = 1 and BO = 0, then plots the whole
                    tree. If segnum = 1 and B0 = 2, then plots the stem and
                    the 1st-order branches. If segnum = 2 and BO = 0, then
                    plots the branch whose index is 2 and all its sub-branches.
  """
  Segs = segment["segments"]
  SChi = segment["children"]

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
    [0.99, 0.41, 0.23],
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
    [0.99, 0.41, 0.23],
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
    [0.99, 0.41, 0.23],
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
    [0.99, 0.41, 0.23],
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

  n = len(Segs)
  '''
  if n < 21:
    col = col * (n // 20 + 1)
  else:
    col = np.random.rand(n, 3)
    col[:, 2] *= 0.2  # Adjusting the blue channel
  '''

  if isinstance(Segs[0], dict):
    #n = np.max([len(s[key]) for key in S.keys()])
    n = len(Segs)
    Seg = {}
    for i in range(n):
      m = len(Segs[i])
      S = []
      for j in range(m):
        s = Segs[i][j]
        S.append(s)
      Seg[i] = np.concatenate(S)
  else:
    Seg = Segs

  if B0 == 0:
    B0 =1000

  S = np.concatenate([cover["ball"][key] for key in Seg[segind]])
  #if len(Seg.keys()) == 1:
  #  S = cover["ball"][Seg[segind]]
  #elif len(Seg.keys()) > 1:
  #  S = np.concatenate([cover["ball"][key] for key in Seg[segind]])

  plt.clf()
  ax = fig.add_subplot(projection='3d')
  ax.scatter(P[S,0], P[S,1], P[S,2], color=col[0], marker='.', s=1)
  forb = S
  if B0 > 1:
    c = SChi[segind]
    i = 1
    while (i < B0) and (len(c)>0):
      C = np.concatenate([cover["ball"][key] for key in np.unique(np.concatenate([Seg[k] for k in c]))])
      C = set_diff(C,forb)
      ax.scatter(P[C,0], P[C,1], P[C,2], color=col[i % len(col)], marker='.', s=0.03)
      '''
      print("******************************************")
      print("******************************************")
      print("******************************************")
      print("******************************************")
      print("******************************************")
      print(f"c: {c}")
      '''
      c = np.unique(np.concatenate([SChi[k] for k in c]))
      i += 1
      forb = np.union1d(forb, C)
  
  plt.axis('equal')
  plt.savefig(f"plots/{prefix}_tree.png")
