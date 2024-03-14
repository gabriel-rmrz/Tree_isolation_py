DEBUG= False
import numpy as np
from tools.define_cut import define_cut
from tools.cut_components import cut_components
from tools.study_components import study_components
from tools.component_classification import component_classification
def segments_num_path(cover, Base, Forb, PathNum):
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
  a = max(200000 , numB/100)                     # Estimate the maximum number of segments
  SBas = {}                                    # The segment bases found
  Segs = {}                                    # The segment found
  SPar = np.zeros((a,2), dtype=np.uint32)      # The parent segment of each segment
  SChi = {}                                    # The childred segmen of each segment

  # Initialize SChi
  SChi[0] = np.zeros(5000, dtype=np.uint32)
  C = np.zeros(200)
  for i in range(1, a):
    SChi[i] = C
  NChi = np.zeros(a, dtype=np.uint32)  # Number of child segments found for each segment

  Fal = np.zeros(numB, dtype='bool')  # Logical false-vector for cover sets
  s=0                                 # The index of the segment under expansion
  b=s                                 # The index of the latest found base

  SBas[s] = Base
  Seg = {}             # The cover set layers in the current segment
  Seg[0] = Base

  ForbAll = np.copy(Fal)      # The forbiden sets
  ForbAll[Forb] = True
  ForbAll[Base] = True
  Forb = ForbAll     # The forbidden sets for the segment under expansion

  Continue = True   # True as long as the component can be semented further
  NewSeg = True     # True if the first Cut for the current segment
  numL = 1            # The number of cover set layers currently in the segment

  # Segmenting stops when there are no more segments to be found

  while Continue and (b < numB-1):
    # Update the forbiden sets
    Forb[Seg[numL-1]] = True

    # Define the study
    Cut = define_cut(Nei, Seg[numL-1], Forb, Fal)
    CutSize=len(Cut)

    if NewSeg:
      NewSeg = False
      numS = min(CutSize, 6)

    # Define the components of cut and study regions
    #numC = 0
    if CutSize > 0:
      CutComps, _ = cut_components(Nei, Cut, CutSize, Fal, Fal)
      if DEBUG:
        print(f"CutComps: {CutComps}")
        print(f"len(CutComps): {len(CutComps)}")
      numC = len(CutComps)
      if numC >1:
        StudyComps, Bases, CompSize, Cont, BaseSize = study_components(Nei, numS, Cut, CutComps, Forb, Fal, Fal)

        if DEBUG:
          print(f"############################")
          print(f"Study components")
          print(f"############################")
          print(f"numC: {numC}")
          print(f"len(Cont) {len(Cont)}")
        numC = len(Cont)
    else:
      numC = 0
    if DEBUG:
      print(f"numC: {numC}")
      

    # Classify study  region components
    if numC == 1:
      # One component, continue expansion of the current segment
      numL += 1
      if len(Cut.shape)>1:
        if len(Cut.shape[1]) > 1:
          Seg[numL-1] = np.transpose(Cut)
      else:
        Seg[numL-1] = Cut
    elif numC>1:
      # Classify the components of the Study region
      Class = component_classification(CompSize, Cont, BaseSize, CutSize)
      if DEBUG:
        print(f"Class: {Class}")

      # Use the bumber of paths to decide which component is the continuation
      N = np.zeros(numC, dtype=np.int32)
      for i in range(numC):
        if len(Bases[i]) > 0:
          N[i] = np.max(PathNum[Bases[i]])
      I = np.argmax(N)
      m = Class[I]
      if m > 0:
        J = Class == 0
        K = Class == 2
        Class[I] = 0
        Class[K] = 0
        Class[J] = 1
      
      # Update information based on the components classification
      for i in range(numC):
        if Class[i] == 1:
          Base = Bases[i]
          ForbAll[Base] = True
          Forb[StudyComps[i]] = True
          J = Forb[Cut]
          Cut = Cut[~J]
          b = b+1
          SBas[b] = Base
          SPar[b,:] = np.array([s, numL-1])
          NChi[s] = NChi[s] + 1
          if DEBUG:
            print(f"s: {s}")
            print(f"NChi[s]: {NChi[s]}")
          SChi[s][NChi[s]] = b
      
      # Define the new cut.
      # If the cut is empty, determine de new base
      if len(Cut) == 0:
        Segs[s] = {key:Seg[key] for key in range(numL)}
        S = np.concatenate([Seg[key] for key in range(numL)])
        ForbAll[S] = True
        if s < b:
          s += 1
          Seg[0] = SBas[s]
          Forb = ForbAll
          NewSeg = True
          numL = 1 
        else:
          Continue = False
      else:
        if len(Cut.shape) > 1:
          Cut = np.transpose(Cut)
        numL += 1
        Seg[numL-1] = Cut
    else:
      # If the study region has zero size, then the current segment is 
      # complete and determine the base of the next segment
      Segs[s] = {key:Seg[key] for key in range(numL)}
      S = np.concatenate([Seg[key] for key in range(numL)])
      ForbAll[S] = True

      if s < b:
        s += 1
        Seg[0] = SBas[s]
        Forb = ForbAll
        NewSeg = True
        numL = 1
      else:
        Continue = False
  Segs = {key:Segs[key] for key in range(b)}
  SPar = SPar[:b,:]
  schi = {key:SChi[key] for key in range(b)}

  # Define output
  SChi = {}
  for i in range(b):
    if NChi[i] >0:
      SChi[i] = (schi[i][:NChi[i]]).astype(np.uint32)
    else:
      SChi[i] = np.asarray([], dtype=np.uint32)
    S = Segs[i]
    for j in range(len(S)):
      S[j] = (S[j]).astype(np.uint32)
    Segs[i] = S
  segment = {}
  segment['segments'] = Segs
  segment['parent'] = SPar
  segment['children'] = SChi

  return segment

