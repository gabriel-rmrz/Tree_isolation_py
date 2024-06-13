import scipy.io
import numpy as np

def mat_struct_to_dict(mat_struct):
  if isinstance(mat_struct, np.ndarray) and mat_struct.dtype.names:
    return{field: mat_struct_to_dict(mat_struct[field][0,0]) for field in mat_struct.dtype.names}
    
  elif isinstance(mat_struct, np.ndarray):
    if mat_struct.size == 1:
      if isinstance(mat_struct.item(), np.ndarray):
        return mat_struct.item().flatten()
      else:
        return mat_struct.item()
    else:
      return [mat_struct_to_dict(element) for element in mat_struct]
  else:
    return mat_struct
def matstruct_to_pydict(filename, structname):
  mat_data = scipy.io.loadmat(filename)
  
  mat_struct = mat_data[structname]
  
  
  py_dict = mat_struct_to_dict(mat_struct)
  #print(f"cover['ball']: {cover['ball']}")
  #print(f"py_dict['ball'][0]: {py_dict['ball'][0]}")
  #print(f"py_dict['ball'][1]: {py_dict['ball'][1]}")
  #print(f"py_dict['ball'][2]: {py_dict['ball'][2]}")
  return py_dict

'''
def main():
  matstruct_to_pydict('cover.mat', 'cover')

if __name__ == '__main__':
  main()

'''
