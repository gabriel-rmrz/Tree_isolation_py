DEBUG=True
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
  Hei = np.asarray(Hei[cover['center']]/100).astype(float) # the heights of the patches (in m)
  GD = inputs['GD0'] # Intial gap distance that can be brigded over with a new link
  DHRel = inputs['DHRel0'] # Initial maximum allowed path_lenght/height relation

  # If not forbidden patches given, then define all patches allowable
  if Forb ==None:
    For = np.zeros(numB, dtype='bool')
  
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
    if Forb != None:
      I = (L < PathLen[N]) & (L / Hei[N] < DHRel) & ~Forb[N]
    else:
      I = (L < PathLen[N]) & (L / Hei[N] < DHRel) 
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
    
    for i in I+1: # Used I+1 instead of range(1,numC+1) to keep the decreasing order of CS
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

        TSS[i] = np.sum([len(Partition[np.unravel_index(key-1,Tree.shape,'F')]) for key in voxels if np.unravel_index(key-1,Tree.shape,'F') in Partition])
        

        if DEBUG:
          print(f"TSS[i]: {TSS[i]}")
          #print(f"Partition.keys(): {Partition.keys()}")

        if ~CompCorresp: # or (CompCor[i-1] ==0) or (CompCor[i-1] > 0 and (TSS[i-1] != TSS0(CompCor[i-1])) :  This line have been partially commented because the variables are not defined
          sets = np.concatenate([Partition[np.unravel_index(key-1,Tree.shape,'F')] for key in voxels if np.unravel_index(key-1,Tree.shape,'F') in Partition])
          treesets = sets[TreeSets[sets]]
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


        if DEBUG:
          print(f"L: {L}")
          print(f"np.argmin(L,axis=1): {np.argmin(L,axis=1)}")
          print(f"np.min(L,axis=1): {np.min(L,axis=1)}")
          print(f"len(np.argmin(L,axis=1)): {len(np.argmin(L,axis=1))}")
          print(f"len(np.min(L,axis=1)): {len(np.min(L,axis=1))}")
          print(f"L.shape: {L.shape}")
          exit()



  EndSet = 0
  PathLen = 0
  return  EndSet, cover, PathLen

