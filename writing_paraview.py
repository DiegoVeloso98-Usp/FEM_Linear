import os

import numpy as np


def appending_tec(SIGMA, DISP, NODES, nnodes, nelems, title):
    var1 = '"displacement_x"'
    var2 = '"displacement_y"'
    var3 = '"stress_x"'
    var4 = '"stress_y"'
    var5 = '"stress_xy"'
    zone = '"zone"'
    posproc = []

    posproc.append(f'TITLE = {title}')
    posproc.append(f'VARIABLES ="x","y",{var1},{var2},{var3},{var4},{var5}')
    posproc.append(f'ZONE T = {zone}')
    posproc.append(f'NODES={nnodes}, ELEMENTS={nelems}, F=FEPOINT, ET=Trilateral')
    
    for i, sig in enumerate(SIGMA):
        disp_x = DISP[2 * i]     # x-component of displacement
        disp_y = DISP[2 * i + 1] # y-component of displacement
        posproc.append(f'{NODES[i][0]} {NODES[i][1]} {disp_x} {disp_y} {sig[0]} {sig[1]} {sig[2]}')
    
    return posproc


def appending_vtk_sigma(SIGMA, DISP, ELEMS, NODES, nnodes, nelems, title, direction, npe):
    posproc = []
    posproc.append(f'# vtk DataFile Version 3.0')
    posproc.append(f'SomeLib output')
    posproc.append(f'ASCII')
    posproc.append(f'DATASET UNSTRUCTURED_GRID')
    posproc.append(f'POINTS {nnodes} float')
    
    for node in NODES:
        posproc.append(f'{node[0]} {node[1]} 0')
    
    if npe == 3:
        posproc.append(f'CELLS {nelems} {nelems * 4}')
        for elem in ELEMS:
            posproc.append(f'{3} {int(elem[0] - 1)} {int(elem[1] - 1)} {int(elem[2] - 1)}')
    elif npe == 6:
        posproc.append(f'CELLS {nelems} {nelems * 7}')
        for elem in ELEMS:
            posproc.append(f'{6} {int(elem[0] - 1)} {int(elem[1] - 1)} {int(elem[2] - 1)} {int(elem[3] - 1)} {int(elem[4] - 1)} {int(elem[5] - 1)}')
    elif npe == 10:
        posproc.append(f'CELLS {nelems} {nelems * 11}')
        for elem in ELEMS:
            posproc.append(f'{10} {int(elem[0] - 1)} {int(elem[1] - 1)} {int(elem[2] - 1)} {int(elem[3] - 1)} {int(elem[4] - 1)} {int(elem[5] - 1)} {int(elem[6] - 1)} {int(elem[7] - 1)} {int(elem[8] - 1)} {int(elem[9] - 1)}')
    
    posproc.append(f'CELL_TYPES {nelems}')
    for i in range(nelems):
        posproc.append(f'5')

    '''
    posproc.append(f'CELL_TYPES {nelems}')
    if npe == 3:
        cell_type = 5  # VTK_TRIANGLE
        for i in range(nelems):
            posproc.append(f'{cell_type}')
    elif npe == 6:
        cell_type = 22  # VTK_QUADRATIC_TRIANGLE
        for i in range(nelems):
            posproc.append(f'{cell_type}')
    elif npe == 10:
        cell_type = 24  # VTK_QUADRATIC_TRIANGLE
        for i in range(nelems):
            posproc.append(f'{cell_type}')

    '''
    
    # Add displacement as vectors
    posproc.append(f'POINT_DATA {nnodes}')
    posproc.append(f'VECTORS displacement float')
    for i in range(nnodes):
        disp_x = DISP[2 * i]     # x-component of displacement
        disp_y = DISP[2 * i + 1] # y-component of displacement
        posproc.append(f'{disp_x} {disp_y} 0')
    
    # Add stress data
    posproc.append(f'SCALARS {title} float')
    posproc.append(f'LOOKUP_TABLE default')
    for sig in SIGMA:
        if direction == 'sigmaX':
            posproc.append(f'{sig[0]}')
        elif direction == 'sigmaY':
            posproc.append(f'{sig[1]}')
        elif direction == 'sigmaXY':
            posproc.append(f'{sig[2]}')
    
    posproc = np.array(posproc)
    return posproc



def exporting_paraview(posproc,path):

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as file:
        file.write('\n'.join(map(str, posproc)))
    print(f"File '{path}' created successfully.")