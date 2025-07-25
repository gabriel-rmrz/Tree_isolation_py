DEBUG=False
import numpy as np
import pandas as pd
import yaml
import laspy
import open3d as o3d
from memory_profiler import memory_usage
from cover_sets_plot import cover_sets_plot
from filtering_plot import filtering_plot
from isolate_trees import isolate_trees
from tools.remove_bottom import remove_bottom
from tools.compute_height import compute_height 

import subprocess
import matplotlib.pyplot as plt
import scipy.io
from plotting.plot_segs import plot_segs


def get_pointcloud(las, isRGB=False):
  point_data = np.stack([las.X, las.Y, las.Z], axis=0).transpose((1,0))
  geom = o3d.geometry.PointCloud()
  geom.points = o3d.utility.Vector3dVector(point_data)

  if isRGB:
    point_color = np.stack([las.red, las.green, las.blue], axis=0).transpose((1,0))
    point_color  = point_color/point_color.max()
    geom.colors = o3d.utility.Vector3dVector(point_color)
  return geom

def run_subprocess():
    input_file_name = 'extracted_points.las'
    matlab_dir = "~/Documents/MATLAB/Matlab_tree_isolation"
    matlab_dir = "/Users/garamire/Work/RemoteSensing/Tree_isolation"
    command = f''' matlab -sd {matlab_dir} -batch "addpath('{matlab_dir}/lasdata', '{matlab_dir}/tools','{matlab_dir}/plotting'); P_temp = lasdata('{input_file_name}'); P = 0.1*[P_temp.x P_temp.y P_temp.z];isolate_trees(P);exit" '''
    subprocess.run(command, shell=True)

def run_isolation():
  with open('configs/inputs.yaml', 'r') as file:
    inputs = yaml.safe_load(file)
  las = laspy.read('extracted_points.las')      
  #las = laspy.read('/home/garamire/work/remote_sensing/treevol/outputs/Region_9.las')      
  #las = laspy.read('/home/garamire/work/remote_sensing/treevol/Area_2_LAS_15.las')      

  P = 0.1 * np.stack([las.x, las.y, las.z], axis=0).transpose((1,0))
  isolate_trees(P)

def main():
  if DEBUG:
    print(f"type(P): {type(P)}")
    exit()
  #subprocess.run(command, shell=True)
  mem_usage = memory_usage(run_subprocess)
  print(f"Memory usage for matlab script: {max(mem_usage)} MB")
  mem_usage = memory_usage(run_isolation)
  print(f"Memory usage for python script: {max(mem_usage)} MB")
  #isolate_trees(P)



if __name__ == '__main__':
  main()

  
