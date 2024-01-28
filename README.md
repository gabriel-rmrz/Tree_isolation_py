# Tree\_isolation
Ported from Pasi Raumonen code

# NOTES
 - For some reason my implementetion of unique_elements (that is use in define_cut-> segments_num_path) is not working propertly. I used np.unique as substitute, but no 100 percent is doing the  same job.
 - The part "Determine the correspondence between the new and the previous components" starting from the line 107 of hte shortest_paths_height.m script is been implemented with variables that doesn't exist and it's not being used.
 - The result after applying the cubical downsamplig is still different. !!! Locate the source
   - Test if the using the operations by blocks could be the cause.
 - Why is the LexOrd vector defined by multiplying C * [1, N(1), N(1) * N(2)] and not [1, N(1), N(2) * N(2)]
 - Line 490 ( R = floor((double(P(PointInd,1:2))-Min(1:2))/SQ)+1;) of the origina method compute\_height might not be correct as the value of Min has been redefined from P.Min() to Ground.Min(), it issues negative values of the indices.
 - We have to keep in mind that we ask for the type of Sub variable in tools/connected_components.py to make decisions.


# TODO
 - Check why the Components value in connected_componens.m:80 and segments_paths_height:244 are different, even if they seem to implement the same solution.
 - See if we can use a generic method instead of tools/unique_elements.py
 - Change the objects in the script compute height with the ones comming from matlab to see where the differences come from.


 - Ask Pasi about the definition in line 442 in filtering_plot.m script
 - Implement lenght method including int and arrays
 - Still missing the implementation of plot_point_cloud called in isolate_trees
 - Still missing the implementation of plot_segs called in isolate_trees
 - Still missing the implementation of define_stem_sections called in isolate_trees

 - There are places where Pasi uses the method any() and we are using sum(), maybe there is a more efficient way in python to do that without performing the full sum.
 - add the parameter of inputs that included in the different scrips (like shortest_path_height.py, at the end, as 'Restriction and iteration parameters  for the shortest path computation'. There are other in the rest of the scripts. They should be defined directly in config/inputs.yaml.
 - Polygonal restriction not working in either version

# DONE
 - remove\_bottoms method is deleting all the entries, which is not consistent with the orinal code. The problem could also be in the compute\_height method. There was a logic error since cubical_down_sampling method. Corrected now.
