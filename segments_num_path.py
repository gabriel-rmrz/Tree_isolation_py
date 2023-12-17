def sements_num_path(cover, Base, Forb, PathNum):
  '''
  ---------------------------------------------------------------------
  SEGMENTS_NUM_PATH.M        Segments the covered point cloud into branches
                              following the shortest paths to the tree base

  Copyright (C) 2013-2023 Pasi Raumonen
  ---------------------------------------------------------------------

  Segments the tree into branches and records their parent-child-relations. 
  Bifurcations are recognized by studying connectivity of a "study"
  region moving along the tree. In case of multiple connected components 
  in "study", the components are classified as the continuation of the 
  current segment and branches. Uses the shortest paths to decide where the
  segment continues in a bifurcation: the continuation is the component
  that has most shortest paths going trough (heuristic being that there are 
  more woody material there as there are more sets). 
  

  Inputs:
  cover         Cover sets
  Base          Base of the tree
  Forb          Cover sets not part of the tree
  PathNum       Number of shortest going through each set 

  Outputs:
  segment       Structure array containing the followin fields:
    segments      Segments found, (n_seg x 1)-cell, each cell contains 
                    a cell array the cover sets
    parent        Parent segment of each segment, (n_seg x 1)-vector,
                    equals to zero if no parent segment
    children      Children segments of each segment, (n_seg x 1)-cell
  '''

  Nei = cover["neighbor"]
  numB = len(Nei)                              # the number of cover sets
  a = max(100000 numB/100)                     # Estimate the maximum number of segments
  SBas = {}                                    # The segment bases found
  Segs = {}                                    # The segment found
  SPar = np.zeros((a,2), dtype=np.uint32)      # The parent segment of each segment
  SChi = {}                                    # The childred segmen of each segment

  # Initialize SChi
  SChi[1] = np.zeros(5000, dtype=np.uint32)
  C = np.zeros(200)
  for i in range(2, a+1):
    SChi[i] = C
  NChi = np.zeros(a)  # Number of child segments found for each segment

  Fal = np.zeros(numB, dtype='bool')  # Logical false-vector for cover sets
  s=1                                 # The index of the segment under expansion
  b=s                                 # The index of the latest found base

  ForbAll = Fal      # The forbiden sets
  ForbAll[Forb] = True
  ForbAll[Base] = True
  Forb = ForbAll     # The bornidden sets for the segment under expansion

  Continue = True   # True as long as the component can be semented further
  NewSeg = True     # True if the first Cut for the current segment
  nl = l            # The number of cover set layers currently in the segment

  # Segmenting stops when there are no more segments to be found






