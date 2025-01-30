import sys

sys.path.append('../../')

from read_2d_data import read_binary_file
from write_2d_data import write_binary_file
from manipulate_2d_data import h2d_to_f2d

data, meta = read_binary_file('c2dtorus_h001.fld')

meta.nsave = 1

write_binary_file('c2dtorus_h_1_001.fld', data, meta)

data, meta = read_binary_file('c2dtorus_h_1_001.fld')

h2d_to_f2d('c2dtorus_h_1_', '../../geom/mesh_test/m2dtorus001.fld', outpattern='c2dtorus')