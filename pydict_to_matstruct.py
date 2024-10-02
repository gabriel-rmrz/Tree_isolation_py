import scipy.io
import numpy as np

def dict_to_mat_struct(d):
  mat_struct = {}
  for key, value in d.items():
    if isinstance(value, dict):
      mat_struct[key] = dict_to_mat_struct(value)
    else:
      mat_struct[key] = np.array(value)

  return mat_struct
def pydict_to_matstruct(my_struct_dict, filename, structname):
  mat_struct = dict_to_mat_struct(my_struct_dict)
  scipy.io.savemat(filename, {structname:mat_struct})

'''
def main():
  my_struct_dict = {
    'name': 'John Doe',
    'age': 30,
    'scores': [95, 85, 88],
    'details': {
      'address': '123 Main St',
      'zip': 12345
      }
    }
  pydict_to_matstruct(my_struct_dict, 'test.mat', 'test')

if __name__ == '__main__':
  main()
'''
