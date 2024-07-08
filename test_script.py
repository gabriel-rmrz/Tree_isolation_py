TEST_ISO= True
TEST_FILT=False
TEST_COV= False
DEBUG=False
import numpy as np
import pandas as pd
import yaml
import laspy
import open3d as o3d
from cover_sets_plot import cover_sets_plot
from filtering_plot import filtering_plot
from isolate_trees import isolate_trees
from test_isolate_trees import test_isolate_trees
from tools.remove_bottom import remove_bottom
from tools.compute_height import compute_height 
#from pycallgraph import PyCallGraph
#from pycallgraph.output import GraphvizOutput

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

def main():
  with open('configs/inputs.yaml', 'r') as file:
    inputs = yaml.safe_load(file)
  #las = laspy.read('extracted_points.las')      
  las = laspy.read('Area_2_LAS_15.las')      
  #print("offset")
  #print(las.header.x_offset)
  #print("scale factor")
  #print(las.header.scales)
  #P = get_pointcloud(las, isRGB=True)
  #P = 0.001 * (np.asarray(P.points).transpose()).astype(float)
  P = 0.001 * np.stack([las.X, las.Y, las.Z], axis=0).transpose((1,0))
  if DEBUG:
    print(f"type(P): {type(P)}")
    exit()

  if TEST_ISO:
    #isolate_trees(P)
    test_isolate_trees(P)
  exit()


  if TEST_FILT:
    Pass, Hei, Ground, Tri = filtering_plot(P,inputs['inputs'])
    #Pass, Hei, Ground, Tri = filtering_plot(np.copy(P),inputs)
    nPass = np.count_nonzero(Pass)
    print(f"nPass: {nPass}")
    #counts, bins = np.histogram(Hei)    
    #plt.stairs(counts, bins)
    plt.hist(Hei,bins=100)
    plt.savefig('plots/Hei.png')
    plt.clf()
    Hei_mat = scipy.io.loadmat('Hei.mat')
    print(Hei_mat['Hei'])
    plt.hist(Hei_mat['Hei'].flatten(),bins=100)
    plt.savefig('plots/Hei_mat.png')

  if TEST_COV:
    Pass, Hei, Ground, Tri = filtering_plot(P,inputs['inputs'])
    P = P[Pass,:]
    cover = cover_sets_plot(P,inputs['inputs'])
    print(f"Plotting...")
    plot_segs(P,cover["neighbor"],1,cover["ball"], 'cover')
    '''
    P_mat = scipy.io.loadmat('P.mat')
    P_mat = P_mat['P']-1 
    #cover_mat = cover_sets_plot(P_mat,inputs['inputs'],test=True)
    cover_mat = cover_sets_plot(P_mat,inputs['inputs'])
    plot_segs(P_mat,cover_mat["neighbor"],1,cover_mat["ball"], 'cover_mat')
    '''


  #print(inputs['inputs'])


  #remove_bottom(P,np.array([0]),np.array([0]),inputs['inputs'])
  #compute_height(P,inputs['inputs'])

  '''
  # Configuration of the output graph
  # graphviz = GraphvizOutput()
  # graphviz.output_file = 'graph.png'
  #with PyCallGraph(output=graphviz):
  isolate_trees(P)
  '''


if __name__ == '__main__':
  main()

  
