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


