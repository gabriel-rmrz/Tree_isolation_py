DEBUG=False
import numpy as np
from scipy.spatial.distance import cdist, squareform
from collections import defaultdict
from tools.cubical_partition import cubical_partition
from tools.connected_components import connected_components
from tools.unique2 import unique2

def shortest_paths_height(P, cover, Hei, Base, BaseDist, inputs, Forb=None):
  # Determines the shortest paths from every point to the base

  nei = cover['neighbor'] # the neighbors of the patches
  NeiDis = cover['NeiDis'] # distances between the neighbors
  numB = len(cover['ball']) # number of patches
  Ce = P[cover['center'],:] # the centers of the patches
  ind = np.asarray(range(numB)).astype(int) # indices from 0 to numB-1
  Hei = np.asarray(Hei[cover['center']]/100.).astype(float) # the heights of the patches (in m)
  GD = inputs['GD0'] # Intial gap distance that can be brigded over with a new link
  DHRel = inputs['DHRel0'] # Initial maximum allowed path_lenght/height relation

  # If not forbidden patches given, then define all patches allowable
  if Forb ==None:
    Forb = np.zeros(numB, dtype='bool')
  
  ## Determine the shoretest paths for as many sets as possible.
  # Paths are not computed for the sets in different components than the base
  # sets and sets whose paths distance compared to height is to large
  # Sort the base according to their distance:

  I = np.argsort(BaseDist)
  BaseDist = BaseDist[I]
  Base = Base[I]
  C = Base[0] # The current set whose patch distance is to be determined
  Unvisited = np.zeros(numB, dtype=np.uint32) # patchees/sets not yet visited by paths
  n = len(Base)
  Unvisited[:n] = Base # Base set have been visited
  UV = np.zeros(numB, dtype='bool') # Unvisited sets in logical vercor also
  UV[Base] = True
  a = 1 # Counts how many unvisited sets have been visited by a shortest path
  b = n # how many unvisited sets there currently are
  J = 1 # index of the set C (a+J-1)
  PathLen = 1000*np.ones(numB,dtype=np.single)
  PathNei = np.zeros(numB, dtype=np.uint32)
  EndSet = np.zeros(numB, dtype=np.uint32)
  PathLen[Base] = BaseDist
  EndSet[Base] = Base
  while a <= b:
    if DEBUG: 
      print(f"C: {C}")
      print(f"nei[C]: {nei[C]}")
    N = nei[C] # neighbors of the set C
    d = NeiDis[C] # the distances of the the neighbors of C
    L = d + PathLen[C]  # Path lenght to the neighbors of C
    # Accept the path via C to N if the path lenghts are shoreter, if the 
    # length/height ratio is small enough, if the sets are not forbidden.
    I = (L < PathLen[N]) & (L / Hei[N] < DHRel) & ~Forb[N]
    N = N[I]
    if len(N)>0:
      # Update the length, neighbor and endset for the sets N:
      if DEBUG:
        print(f"PathLen: {PathLen}")
        print(f"PathLen[C]: {PathLen[C]}")
        print(f"C: {C}")
        print(f"d: {d}")
        print(f"L: {L}")
        print(f"type(L): {type(L)}")
        print(f"I: {I}")
        print(f"type(I): {type(I)}")
        print(f"N: {N}")
      PathLen[N] = L[I]
      PathNei[N] = C
      EndSet[N] = EndSet[C]

    # Define C as visited and N as unvisited
    if J>1: 
      Unvisited[a-J-1-1] = Unvisited[a-1]

    a += 1
    UV[C] = False
    if len(N) > 0:
      I = UV[N]
      N = N[~I] 
      if len(N) > 0:
        Unvisited[b+1-1: b+len(N)] = N
        b = b + len(N)
        UV[N] = True

    # Determine the next set "C" form the unvisited sets that has the minimum
    # path length
      if DEBUG:
        print(f"PathLen: {PathLen}")
        print(f"Unvisited: {Unvisited}")
        print(f"a: {a}")
        print(f"b: {b}")
        print(f"PathLen[Unvisited[int(a-1):int(b)]]: {PathLen[Unvisited[int(a-1):int(b)]]}")
    if len(PathLen[Unvisited[int(a-1):int(b)]] > 0):
      J = np.argmin(PathLen[Unvisited[int(a-1):int(b)]])
      C = Unvisited[a + J -1 -1]
  
  ## Iteratively expand trees base on the shortest paths
  # Partition of the cover sets into search space

  Partition, CC, info = cubical_partition(Ce, GD, nargout=3)
  Partition = defaultdict(np.ndarray, Partition)
  CC = CC.astype(np.double)
  Voxels = CC[:,0] + info[3]*(CC[:,1] -1) + info[3]*info[4]*(CC[:,2] -1)
  CompCorresp = False
  TreeSets = PathLen < 1000 # set belonging to trees, their path length has been determined

  while np.any(~TreeSets) and DHRel <= inputs["MaxDHRel"]:
    ## Determine the sets close to the tree sets
    # Define "Tree" as a spatial distribution of "TreeSets" in voxel space
    cc = (np.unique(CC[TreeSets,:], axis=0)).astype(int)
    if DEBUG:
      print(f"info: {info}")
      print(f"info[3:6]: {info[3:6]}")
    Tree = np.zeros(info[3:6].astype(int), dtype='bool') # tree sets in space/voxelization
    n = len(cc[:,0])
    for i  in range(n):
      Tree[cc[i,0]-1-1:cc[i,0]+1,cc[i,1]-1-1:cc[i,1]+1,cc[i,2]-1-1:cc[i,2]+1] = True

    #Determine the components of the other sets (non-tree sets)
    Other = ind[~TreeSets]
    if DEBUG:
      print(f"Voxels[Other]: {Voxels[Other]}")
      print(f"(Voxels[Other]).astype(int): {(Voxels[Other]).astype(int)}")
    VoxOtherSets= (Voxels[Other]).astype(int)
    I = Tree[np.unravel_index(VoxOtherSets-1,Tree.shape,'F')]
    Other = Other[I]
    Comps, CS = connected_components(nei, Other, 1)
    if DEBUG:
      print(f"Connected components done")

    I = np.argsort(CS)[::-1]
    CS = CS[I]
    if DEBUG:
      print(f"Comps: {Comps}")
    #Comps = Comps[I]
    numC = len(Comps)

    # Determine the correspondence between the new and the previous components
    #The following part is not implemented in the original code.
    '''
    if CompCorresp:
      CoS0 = np.zeros(numB, dtype=np.int32)
      for i in range(1,len(Comps0)):
    '''

    numT = np.count_nonzero(TreeSets)
    numO = np.count_nonzero(~TreeSets)

    ## Check if component can be joined to a tree
    TSS = np.zeros(numC, dtype=np.int32)
    
    for i in I:#range(1,numC+1): # Used I+1 instead of range(1,numC+1) to keep the decreasing order of CS
      Comp = Comps[i] # sets in the component

      if np.any(~TreeSets[Comp]):
        # Determine the closeby tree sets:
        V = unique2(Voxels[Comp])
        n = len(V)
        if DEBUG:
          print(f"V: {V}")
        a = info[3]*info[4]
        Vox = np.zeros(9*n)
        for j in range(1, n+1):
          v = V[j-1]
          Vox[(j-1)*9:j*9] = [v, v-1, v+1, v-a, v-a-1, v-a+1, v+a, v+a-1, v+a+1]

        voxels = (unique2(Vox)).astype(int)
        if DEBUG:
          for key in voxels:
            print(f"voxel.key: {key}")

        TSS[i-1] = np.sum([len(Partition[np.unravel_index(key-1,Tree.shape,'F')]) for key in voxels if np.unravel_index(key-1,Tree.shape,'F') in Partition])
        

        if DEBUG:
          print(f"TSS[i-1]: {TSS[i-1]}")
          #print(f"Partition.keys(): {Partition.keys()}")

        if ~CompCorresp: # or (CompCor[i-1] ==0) or (CompCor[i-1] > 0 and (TSS[i-1] != TSS0(CompCor[i-1])) :  This line have been partially commented because the variables are not defined
          sets = [Partition[np.unravel_index(key-1,Tree.shape,'F')] for key in voxels if np.unravel_index(key-1,Tree.shape,'F') in Partition]
          if len(sets)>0:
            sets = np.concatenate(sets)
            treesets = sets[TreeSets[sets]]
          else:
            sets = np.array([])
            treesets = np.array([])
        else:
          treesets = np.array([])
  

        if DEBUG:
          print(f"sets: {sets}")
  
        if len(treesets) > 0:
          ## Determine the possible connections between the component and tree
          # Compute the distances between the component and tree
          if len(Comp)*len(treesets) < 1e8:
            L = cdist(Ce[Comp,:],Ce[treesets,:])

            # Select the closest set from tree set for each set in the componen
            D1 = np.min(L,axis=1)
            I = np.argmin(L,axis=1)
          else:
            m = len(Comp)
            I = np.zeros(m)
            D1 = I
            for j in range(1, int(np.ceil(m/1000.))+1):
              a = (j-1)*1000 + 1
              b = min(j*1000,m)
              L = cdist(Ce[Comp[a-1:b],:], Ce[treesets,:])
              # Select the closest set from tree set for each set in the component
              D2 = np.min(L, axis=1)
              I2 = np.argmin(L, axis=1)
              D1[a-1:b] = D2
              I1[a-1:b] = I2
          TreeSet = treesets[I]
          if DEBUG:
            print(f"L: {L}")
            print(f"L.shape: {L.shape}")
            print(f"Comp.shape: {Comp.shape}")

          # Compute the new path distances
          L = PathLen[TreeSet] + D1
          # Compute the ratio of path distance and heigh of the sets
          DH = L/Hei[Comp]

          # Sort the sets from the closest to the furthest
          SOrd = np.argsort(D1)
          D1 = D1[SOrd]
          if (D1[0]< GD) and np.any(DH < DHRel):
            n = len(D1)
            TSA = np.zeros(n)
            CSA = np.zeros(n)
            for j in range(1, n+1):
              if (D1[j-1] < GD) and (DH[SOrd[j-1]] < DHRel) and ~TreeSets[Comp[SOrd[j-1]]]:
                ## Join the sets and expand
                CS = Comp[SOrd[j-1]] # the set in the component
                TS = TreeSet[SOrd[j-1]] # the set in the treeset
                dist = D1[j-1] # the distance between these sets
                TSA[j-1] = TS
                CSA[j-1] = CS

                # Join the component set to the tree, update neighbors and 
                # other information:
                if DEBUG:
                  print(f"TS: {TS}")
                  print(f"CS: {CS}")
                  print(f"nei[TS]: {nei[TS]}")
                  print(f"nei[CS]: {nei[CS]}")

                nei[TS] = np.concatenate([nei[TS], [CS]])
                nei[CS] = np.concatenate([nei[CS], [TS]])
                if DEBUG:
                  print(f"TS: {TS}")
                  print(f"CS: {CS}")
                  print(f"nei[TS]: {nei[TS]}")
                  print(f"nei[CS]: {nei[CS]}")

                NeiDis[TS] = np.concatenate([NeiDis[TS], [dist]])
                NeiDis[CS] = np.concatenate([NeiDis[CS], [dist]])
                PathLen[CS] = PathLen[TS] + dist
                PathNei[CS] = TS
                EndSet[CS] = EndSet[TS]

                # Expand as much as possibl;e
                unvisited = np.setdiff1d(nei[CS], TS)
                C = CS
                m = len(unvisited)
                Unvisited[:m] = unvisited
                UV[unvisited] = True
                a = 1
                b = m 
                J= 1
                while a <=b:
                  N = nei[C]
                  d = NeiDis[C]
                  L = PathLen[C] + d
                  if DEBUG:
                    print(f"N : {N}")
                    print(f"Hei[N] : {Hei[N]}")
                    print(f"PathLen[N] : {PathLen[N]}")
                    print(f"Forb[N] : {Forb[N]}")
                    print(f"DHRel : {DHRel}")
                  I = (L<PathLen[N]) & (L/Hei[N] < DHRel) & ~Forb[N]
                  N = N[I]
                  if len(N) >0:
                    PathLen[N] = L[I]
                    PathNei[N] = C
                    EndSet[N] = EndSet[C]
                  
                  if J >1:
                    Unvisited[a+J-1-1] = Unvisited[a-1]
                  if DEBUG:
                    print(f"a : {a}")
                    print(f"b: {b}")
                  
                  a +=1
                  UV[C] = False
                  if len(N)>0:
                    I = UV[N]
                    N = N[~I]
                    if len(N)>0:
                      Unvisited[b+1-1:b+len(N)] = N
                      b = b+len(N)
                      UV[N] = True

                  if DEBUG:
                    print(f"a (after): {a}")
                    print(f"b (after): {b}")
                  if a <=b:
                    J = np.argmin(PathLen[Unvisited[int(a-1):int(b)]])
                    C = Unvisited[a+J-1-1]
                  else:
                    J = np.array([])
                    C = np.array([])
    Comps0 = Comps              
    CompCorresp = True

    numO0 = numO
    numO = np.count_nonzero(~TreeSets)

    ## Increase maximum gap distance and path distance height ratio
    if numO == numO0:
      GD = GD + inputs["GD"]
      DHRel = DHRel + inputs["DHRel"]
      CompCorresp = False

      ## Partition of cover sets into search space
      Partition, CC, info = cubical_partition(Ce, GD)
      CC = CC.astype(float)
      Voxels = CC[:,0] + info[3]*(CC[:,1] -1) + info[3]*info[4]*(CC[:,2] -1)

      ## Determine the shortest paths again
      C = Base[0]
      Unvisited = np.zeros(numB, dtype=np.uint32)
      n = len(Base)
      Unvisited[:n] = Base
      UV = np.zeros(numB, dtype='bool')
      UV[Base] = True
      a = 1
      b = n
      J = 1
      PathLen = 1000 * np.ones(numB, dtype=np.single) # The path distances
      PathNei = np.zeros(numB, dtype=np.uint32)
      EndSet = np.zeros(numB, dtype=np.uint32)
      PathLen[Base] = BaseDist
      EndSet[Base] = Base
      while a <=b:
        N = nei[C]
        d = NeiDis[C]
        L = PathLen[C] + d
        if DEBUG:
          print(f"N : {N}")
        I = (L<PathLen[N]) & (L/Hei[N] < DHRel) & ~Forb[N]
        N = N[I]
        if len(N) >0:
          PathLen[N] = L[I]
          PathNei[N] = C
          EndSet[N] = EndSet[C]
        
        if J >1:
          Unvisited[a+J-1-1] = Unvisited[a-1]
        if DEBUG:
          print(f"a : {a}")
          print(f"b: {b}")
        
        a +=1
        UV[C] = False
        if len(N)>0:
          I = UV[N]
          N = N[~I]
          if len(N)>0:
            Unvisited[b+1-1:b+len(N)] = N
            b = b+len(N)
            UV[N] = True

        if DEBUG:
          print(f"a (after): {a}")
          print(f"b (after): {b}")
        if a <=b:
          J = np.argmin(PathLen[Unvisited[int(a-1):int(b)]])
          C = Unvisited[a+J-1-1]
        else:
          J = np.array([])
          C = np.array([])
  
  # Update the neighbor relation
  cover["neighbor"] = nei
  cover["NeiDis"] = NeiDis

  return  EndSet, cover, PathLen

