import open3d as o3d
#import pandas as pd
#import ellipsis as el
import matplotlib.pyplot as plt
import numpy as np


def plot_point_cloud(P, prefix=""):
  '''
  Creates 3D plots of a point cloud from different perspectives
  '''
  #geom = o3d.geometry.PointCloud()
  #geom.points = o3d.utility.Vector3dVector(P)

  #viewer = o3d.visualization.Visualizer()
  #viewer.create_window()
  #viewer.add_geometry(geom)
  #img = viewer.capture_screen_float_buffer(True)
  #plt.imshow(np.asarray(img))
  fig = plt.figure(figsize=(12,8))

  # Perspective view
  ax1 = fig.add_subplot(221, projection='3d')
  ax1.scatter(P[:,0], P[:,1], P[:,2], c='b', marker='.', s=1)
  ax1.set_title('Perspective View')

  #Top view
  ax2 = fig.add_subplot(222)
  ax2.scatter(P[:,0], P[:,1], c='b', marker='.', s=1)
  ax2.set_title('Top View')
  ax2.set_xlabel('X-axis')
  ax2.set_ylabel('Y-axis')
  ax2.axis('equal')

  #Side view
  ax3 = fig.add_subplot(223)
  ax3.scatter(P[:,0], P[:,2], c='b', marker='.', s=1)
  ax3.set_title('Side View')
  ax3.set_xlabel('X-axis')
  ax3.set_ylabel('Z-axis')
  ax3.axis('equal')

  #From view
  ax3 = fig.add_subplot(224)
  ax3.scatter(P[:,1], P[:,2], c='b', marker='.', s=1)
  ax3.set_title('Front View')
  ax3.set_xlabel('Y-axis')
  ax3.set_ylabel('Z-axis')
  ax3.axis('equal')


  plt.tight_layout()

  plt.savefig(f"plots/{prefix}_point_cloud.png")
  #viewer.destroy_window()
