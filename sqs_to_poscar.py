'''
This program reads bestsqs file from user and generates POSCAR file 
To run use "python sqs_to_poscar.py bestsqs.out"
'''


import sys
import numpy as np


# Reading bestsqs.out
data = []
file_name = sys.argv[1]
with open(file_name, 'r') as reader:
    data = reader.read().rstrip().split("\n")


def extract_vals(s_i, e_i):
    vec = []
    for i in range(s_i, e_i):
        vec.append(data[i].split())
    return vec


unit_cell = np.array(extract_vals(0, 3)).astype(np.float)
lattice_vectors = np.array(extract_vals(3, 6)).astype(np.float)
atomic_positions = np.array(extract_vals(6, len(data)))

# Arrange atoms according to elements
elem_arr = np.unique(atomic_positions[:, 3])
a_pos = np.array([['', '', '', '']])
num_atoms_arr = []
for elem in elem_arr:
    tmp_arr = atomic_positions[atomic_positions[:, 3] == elem]
    num_atoms_arr.append(len(tmp_arr))
    a_pos = np.concatenate(
        (a_pos, tmp_arr), axis=0)
a_pos = a_pos[1:]


# multiply to obtain lattice vector and atomic positions
l_vec = np.matmul(lattice_vectors, unit_cell)
a_pos = np.matmul(a_pos[:, 0:3].astype(np.float), unit_cell)


# Converting data to POSCAR
str_arr = [
    'POSCAR',
    '1.0'
]

for i in range(0, 3):
    str_arr.append(
        f'  {l_vec[i,0]:.10f}  {l_vec[i,1]:.10f}  {l_vec[i,2]:.10f}',)

s = '  '
for elem in elem_arr:
    s += (elem+' ')
str_arr.append(s)

s = '  '
for num_atoms in num_atoms_arr:
    s += (str(num_atoms)+' ')
str_arr.append(s)

str_arr.append('Cartesian')

for a in a_pos:
    str_arr.append(
        f'  {a[0]:.10f}  {a[1]:.10f}  {a[2]:.10f}',)

# Save to file
with open('POSCAR', 'w') as f:
    f.write('\n'.join(str_arr))
