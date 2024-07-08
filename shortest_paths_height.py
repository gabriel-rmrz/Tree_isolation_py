DEBUG=False
import numpy as np
from scipy.spatial.distance import cdist, squareform
from collections import defaultdict
from tools.cubical_partition import cubical_partition
from tools.connected_components import connected_components
from tools.unique2 import unique2
if DEBUG:
  from tabulate import tabulate
  def print_boolean_matrix(matrix):
    # Convert boolean values to int for a cleaner output
    int_matrix = [[int(cell) for cell in row] for row in matrix]
    print(tabulate(int_matrix, tablefmt="plain"))
  def print_boolean_array(array):
    # Convert the array to a numpy array if it isn't already
    np_array = np.array(array, dtype=int)
    # Convert the array to a string without any spaces between elements and print
    print(''.join(map(str, np_array)))

def shortest_paths_height(P, cover, Hei, Base, BaseDist, inputs, Forb=None):
  # Determines the shortest paths from every point to the base

  nei = cover['neighbor'] # the neighbors of the patches
  NeiDis = cover['NeiDis'] # distances between the neighbors
  numB = len(cover['ball']) # number of patches
  Ce = P[cover['center'],:] # the centers of the patches
  ind = np.asarray(range(numB)).astype(np.uint32) # indices from 0 to numB-1
  Hei = np.asarray(Hei[cover['center']]/100.).astype(np.double) # the heights of the patches (in m)
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
  UV = np.zeros(numB, dtype='bool') # Unvisited sets in logical vector also
  UV[Base] = True
  a = 1 # Counts how many unvisited sets have been visited by a shortest path
  b = n # how many unvisited sets there currently are
  J = 0 # index of the set C (a+J-1)
  PathLen = 1000*np.ones(numB,dtype=np.double)
  PathNei = np.zeros(numB, dtype=np.uint32)
  EndSet = np.zeros(numB, dtype=np.uint32)
  PathLen[Base] = BaseDist
  EndSet[Base] = Base
  while a <= b:
    N = cover['neighbor'][C] # neighbors of the set C
    d = NeiDis[C] # the distances of the the neighbors of C
    L = d + PathLen[C]  # Path lenght to the neighbors of C
    # Accept the path via C to N if the path lenghts are shoreter, if the 
    # length/height ratio is small enough, if the sets are not forbidden.
    I = (L < PathLen[N]) & ((L / Hei[N]) < DHRel) & ~Forb[N]
    N = N[I]
    if len(N)>0:
      # Update the length, neighbor and endset for the sets N:
      PathLen[N] = L[I]
      PathNei[N] = C
      EndSet[N] = EndSet[C]

    # Define C as visited and N as unvisited
    if J>0: 
      Unvisited[a+J-1] = Unvisited[a-1]

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
    if int(a-1)<int(b):
      J = np.argmin(PathLen[Unvisited[int(a-1):int(b)]]) 
      C = Unvisited[a + J -1]
  
  ## Iteratively expand trees base on the shortest paths
  # Partition of the cover sets into search space

  Partition, CC, info = cubical_partition(Ce, GD, method=1)
  if DEBUG:
    print(f"CC.shape: {CC.shape}")
    print(f"Partition.keys(): {Partition.keys()}")
  #Partition = defaultdict(tuple, Partition)
  CC = CC.astype(np.double)
  Voxels = CC[:,0] + info[3]*(CC[:,1] -1) + info[3]*info[4]*(CC[:,2] -1)
  CompCorresp = False
  TreeSets = np.array(PathLen < 1000, dtype='bool') # set belonging to trees, their path length has been determined
  if DEBUG:
    print(f"len(PathLen): {len(PathLen)}")
    print(f"len(TreeSets): {len(TreeSets)}")
    print(f"sum(TreeSets): {sum(TreeSets)}")

  Comps0={}
  while np.any(~TreeSets) and DHRel <= inputs["MaxDHRel"]:
    ## Determine the sets close to the tree sets
    # Define "Tree" as a spatial distribution of "TreeSets" in voxel space
    cc = (np.unique(CC[TreeSets,:], axis=0)).astype(np.int32)
    if DEBUG:
      print(f"cc: {cc}")
      print(f"len(cc): {len(cc)}")

    Tree = np.zeros(info[3:6].astype(int), dtype='bool') # tree sets in space/voxelization
    n = len(cc[:,0])
    for i  in range(n):
      Tree[cc[i,0]-1-1:cc[i,0]+1,cc[i,1]-1-1:cc[i,1]+1,cc[i,2]-1-1:cc[i,2]+1] = True # this should be ok

    if DEBUG:
      #print(f"Tree[:,:,20]: {Tree[:,:,20]}")
      print(f"Tree.shape: {Tree.shape}")
      #print_boolean_matrix(Tree[:,:,20])

    #Determine the components of the other sets (non-tree sets)
    Other = ind[~TreeSets]
    if DEBUG:
      print(f"len(Other): {len(Other)}")
      print(f"Other: {Other}")
    VoxOtherSets= (Voxels[Other]).astype(np.uint32)
      
    
    cubes_voxOtherSets = np.unravel_index(VoxOtherSets-1,Tree.shape,'F')
    #cubes_voxOtherSets = tuple(arr -1 for arr in cubes_voxOtherSets) This does not seem to be necesary

    I = Tree[cubes_voxOtherSets] # This should be ok
    Other = Other[I]
    if DEBUG:
      print(f"Other: {Other}")
      print(f"VoxOtherSets: {VoxOtherSets}")
      print(f"len(Other): {len(Other)}")
      print(f"len(VoxOtherSets): {len(VoxOtherSets)}")
    Comps, CS = connected_components(nei, np.copy(Other), 1)
    #I = np.argsort(CS)[::-1]
    I0 = np.argsort(-CS)
    CS = CS[I0]

    Comps = {I0[key]: value for key, value in Comps.items()}
    #Comps = { key: Comps[key] for key in I0}
    numC = len(Comps)
    if DEBUG:
      print(f"len(Comps): {len(Comps)}")


    # Determine the correspondence between the new and the previous components
    if CompCorresp:
      CoS0 = np.zeros(numB, dtype=np.int32)
      for i in range(len(Comps0)):
        CoS0[Comps0[i]] = i
      CompCor = np.zeros(numC)
      for i in range(numC):
        C = CoS0[Comps[i]]
        if C[0] > 0 and np.all(C==C[0]):
          C = C[0]
          if len(Comps0[C]) == len(Comps[i]):
            CompCor[i] = C
      TSS0 = TSS
    numT = np.count_nonzero(TreeSets)
    numO = np.count_nonzero(~TreeSets)
    if DEBUG:
      print(f"numT: {numT}")
      print(f"numO: {numO}")

    ## Check if component can be joined to a tree
    TSS = np.zeros(numC, dtype=np.int32)
    
    for i in range(numC):#range(1,numC+1): # Used I+1 instead of range(1,numC+1) to keep the decreasing order of CS
      Comp = Comps[i] # sets in the component

      if np.any(~TreeSets[Comp]):
        # Determine the closeby tree sets:
        V = unique2(Voxels[Comp])
        n = len(V)
        a = info[3]*info[4]
        Vox = np.zeros(9*n)
        for j in range(n):
          v = V[j]
          Vox[j*9:(j+1)*9] = [v, v-1, v+1, v-a, v-a-1, v-a+1, v+a, v+a-1, v+a+1] 
        if DEBUG:
          print(f"Vox: {Vox}")
          print(f"len(Vox): {len(Vox)}")

        voxels = (unique2(Vox)).astype(np.uint32)
        par_keys = []
        for vo in voxels:
          k_vox =np.unravel_index(vo-1,Tree.shape,'F')
          k_vox = tuple([k_vox[0]+1, k_vox[1]+1, k_vox[2] + 1])
          par_keys.append(k_vox)

        TSS[i] = np.sum([len(Partition[key]) for key in par_keys])
        if ~CompCorresp or (CompCor[i] ==0) or (CompCor[i] > 0 and (TSS[i] != TSS0(CompCor[i]))) :  
          sets = np.concatenate([Partition[key] for key in par_keys]).astype(np.int32)
          if DEBUG:
            print(f"sets: {sets}")
            print(f"len(sets): {len(sets)}")
          treesets = sets[TreeSets[sets]]
          if DEBUG:
            print(f"treesets: {treesets}")
          '''
          if len(sets)>0:
            sets = np.concatenate(sets)
            treesets = sets[TreeSets[sets]]
          else:
            sets = np.array([])
            treesets = np.array([])
          '''
        else:
          treesets = np.array([])
  
        if len(treesets) > 0:
          ## Determine the possible connections between the component and tree
          # Compute the distances between the component and tree
          if len(Comp)*len(treesets) < 1e8:
            L = cdist(Ce[Comp,:],Ce[treesets,:])
            if DEBUG:
              print(f"L: {L}")
              print(f"len(L): {len(L)}")
    

            # Select the closest set from tree set for each set in the componen
            D1 = np.min(L,axis=1)
            I = np.argmin(L,axis=1)
            #ID1 = np.argsort(D1)
            #D1 = D1[ID1]
            #I = I[ID1]
            if DEBUG:
              print(f"D1: {D1}")
              print(f"len(D1): {len(D1)}")
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
              #I2 = np.argmin(L, axis=1)
              #ID2 = np.argsort(D2)
              #D2 = D2[ID2]
              #I2 = I2[ID2]
              D1[a-1:b] = D2
              I[a-1:b] = I2
          TreeSet = treesets[I]
          if DEBUG:
            print(f"D1: {D1}")
            print(f"treesets: {treesets}")
            print(f"TreeSet: {TreeSet}")
            print(f"len(TreeSet): {len(TreeSet)}")

          # Compute the new path distances
          L = PathLen[TreeSet] + D1
          if DEBUG:
            print(f"L: {L}")
          # Compute the ratio of path distance and heigh of the sets
          DH = L/Hei[Comp]

          # Sort the sets from the closest to the furthest
          SOrd = np.argsort(D1)
          D1 = D1[SOrd]
          if (D1[0]< GD) and np.any(DH < DHRel):
            n = len(D1)
            TSA = np.zeros(n)
            CSA = np.zeros(n)
            for j in range(n):
              if (D1[j] < GD) and (DH[SOrd[j]] < DHRel) and ~TreeSets[Comp[SOrd[j]]]:
                ## Join the sets and expand
                CS = Comp[SOrd[j]] # the set in the component
                TS = TreeSet[SOrd[j]] # the set in the treeset
                dist = D1[j] # the distance between these sets
                TSA[j] = TS
                CSA[j] = CS

                # Join the component set to the tree, update neighbors and 
                # other information:
                nei[TS] = np.concatenate([nei[TS], [CS]])
                nei[CS] = np.concatenate([nei[CS], [TS]])
                NeiDis[TS] = np.concatenate([NeiDis[TS], [dist]])
                NeiDis[CS] = np.concatenate([NeiDis[CS], [dist]])
                PathLen[CS] = PathLen[TS] + dist
                PathNei[CS] = TS
                EndSet[CS] = EndSet[TS]

                # Expand as much as possibl;e
                unvisited = np.setdiff1d(nei[CS], TS)
                if DEBUG:
                  print(f"unvisited: {unvisited}")


                C = CS
                m = len(unvisited)
                Unvisited[:m] = unvisited
                UV[unvisited] = True
                a = 1
                b = m 
                J=0 
                while a <=b:
                  N = nei[C]
                  d = NeiDis[C]
                  L = PathLen[C] + d
                  I = (L<PathLen[N]) & (L/Hei[N] < DHRel) & ~Forb[N]
                  N = N[I]
                  if len(N) >0:
                    PathLen[N] = L[I]
                    PathNei[N] = C
                    EndSet[N] = EndSet[C]
                    TreeSets[N] = True
                  
                  if J >0:
                    Unvisited[a+J-1] = Unvisited[a-1]
                  
                  a +=1
                  UV[C] = False
                  if len(N)>0:
                    I = UV[N]
                    N = N[~I]
                    if len(N)>0:
                      Unvisited[b+1-1:b+len(N)] = N
                      b = b+len(N)
                      UV[N] = True
                  if a <=b:
                    J = np.argmin(PathLen[Unvisited[int(a-1):int(b)]]) 
                    C = Unvisited[a+J-1]
    Comps0 = Comps
    CompCorresp = True

    numO0 = numO
    numO = np.count_nonzero(~TreeSets)
    if DEBUG:
      print(f"Comps: {Comps}")
      print(f"len(Comps): {len(Comps)}")
      print(f"numO: {numO}")
      print(f"numO0: {numO0}")

    ## Increase maximum gap distance and path distance height ratio
    if numO == numO0:
      GD = GD + inputs["GD"]
      DHRel = DHRel + inputs["DHRel"]
      CompCorresp = False

      ## Partition of cover sets into search space
      Partition, CC, info = cubical_partition(Ce, GD, method=1)
      CC = CC.astype(np.double)
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
      J = 0
      PathLen = 1000 * np.ones(numB, dtype=np.single) # The path distances
      PathNei = np.zeros(numB, dtype=np.uint32)
      EndSet = np.zeros(numB, dtype=np.uint32)
      PathLen[Base] = BaseDist
      EndSet[Base] = Base
      while a <=b:
        N = nei[C]
        d = NeiDis[C]
        L = PathLen[C] + d
        I = (L<PathLen[N]) & (L/Hei[N] < DHRel) & ~Forb[N]
        N = N[I]
        if len(N) >0:
          PathLen[N] = L[I]
          PathNei[N] = C
          EndSet[N] = EndSet[C]
        
        if J >0:
          Unvisited[a+J-1] = Unvisited[a-1]
        a +=1
        UV[C] = False
        if len(N)>0:
          I = UV[N]
          N = N[~I]
          if len(N)>0:
            Unvisited[b+1-1:b+len(N)] = N
            b = b+len(N)
            UV[N] = True
        if a <=b:
          J = np.argmin(PathLen[Unvisited[int(a-1):int(b)]])
          C = Unvisited[a+J-1]
  
  # Update the neighbor relation
  cover["neighbor"] = nei
  cover["NeiDis"] = NeiDis

  return  EndSet, cover, PathLen

