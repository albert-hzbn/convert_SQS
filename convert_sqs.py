'''
This program reads bestsqs file from user and generates structure file 
for quantum espresso or VASP format
To run use "python3 convert_sqs.py bestsqs.out"
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


unit_cell = np.array(extract_vals(0, 3)).astype(float)
lattice_vectors = np.array(extract_vals(3, 6)).astype(float)
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
a_pos = np.matmul(a_pos[:, 0:3].astype(float), unit_cell)

out_file = ""
str_arr = []
str_out = int(input("Enter 1 for QE format and 2 for VASP format: "))
if(str_out == 1):
    # Converting data to QE format
    out_file = "qe_inp.in"
    str_arr = [
        '&CONTROL',
        '/',
        '&SYSTEM',
        '  ibrav = 0',
        f'  nat = {len(a_pos)}',
        f'  ntyp = {len(elem_arr)}',
        '/',
        '&ELECTRONS',
        '/',
        'ATOMIC_SPECIES'
    ]
    for elem in elem_arr:
        str_arr.append(f"  {elem} xx xxxx")
    
    str_arr.append("CELL_PARAMETERS angstrom")
    for i in range(0, 3):
        str_arr.append(
            f'  {l_vec[i,0]:.10f}  {l_vec[i,1]:.10f}  {l_vec[i,2]:.10f}',)

    str_arr.append('ATOMIC_POSITIONS angstrom')
    ind = 0
    arr_ind = 0
    for a in a_pos:
        if(ind == num_atoms_arr[arr_ind]):
            arr_ind += 1
            ind = 0
        str_arr.append(
            f'{elem_arr[arr_ind]}  {a[0]:.10f}  {a[1]:.10f}  {a[2]:.10f}',)
        ind += 1
    str_arr.append("K_POINTS automatic")
    str_arr.append("  1 1 1 0 0 0")
    

elif(str_out == 2):
    # Converting data to POSCAR
    out_file = "POSCAR"

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
with open(out_file, 'w') as f:
    f.write('\n'.join(str_arr))
