import numpy as np
import pandas as pd
import yaml
import laspy
import open3d as o3d
from cover_sets_plot import cover_sets_plot
from filtering_plot import filtering_plot
from isolate_trees import isolate_trees
from tools.remove_bottom import remove_bottom
from tools.compute_height import compute_height 

def get_pointcloud(las, isRGB=False):
  point_data = np.stack([las.X, las.Y, las.Z], axis=0).transpose((1,0))
  geom = o3d.geometry.PointCloud()
  geom.points = o3d.utility.Vector3dVector(point_data)

  if isRGB:
    point_color = np.stack([las.red, las.green, las.blue], axis=0).transpose((1,0))
    point_color  = point_color/point_color.max()
    geom.colors = o3d.utility.Vector3dVector(point_color)
  return geom

def main():
  with open('configs/inputs.yaml', 'r') as file:
    inputs = yaml.safe_load(file)
  las = laspy.read('extracted_points.las')      
  #print("offset")
  #print(las.header.x_offset)
  #print("scale factor")
  #print(las.header.scales)
  #P = get_pointcloud(las, isRGB=True)
  #P = 0.001 * (np.asarray(P.points).transpose()).astype(float)
  P = 0.001 * np.stack([las.X, las.Y, las.Z], axis=0).transpose((1,0))

  #print(inputs['inputs'])


  #cover_sets_plot(P,inputs['inputs'])
  #filtering_plot(P,inputs['inputs'])
  isolate_trees(P)
  #remove_bottom(P,np.array([0]),np.array([0]),inputs['inputs'])
  #compute_height(P,inputs['inputs'])

if __name__ == '__main__':
  main()

  
