import open3d as o3d
#import pandas as pd
#import ellipsis as el
import matplotlib.pyplot as plt
import numpy as np


def plot_point_cloud(P):
  geom = o3d.geometry.PointCloud()
  geom.points = o3d.utility.Vector3dVector(P)

  viewer = o3d.visualization.Visualizer()
  viewer.create_window()
  viewer.add_geometry(geom)
  img = viewer.capture_screen_float_buffer(True)
  plt.imshow(np.asarray(img))
  plt.savefig('test_o3d.png')
  viewer.destroy_window()
