import numpy as np
from shape import pascal_triangle, shapeFunc

def fem(KGLOBAL, FGLOBAL):
    U = np.linalg.solve(KGLOBAL, FGLOBAL)
    return U

def particleDisplacement(U,SIGMA_full, PARTICLES, KSIETA, ELEMS, gaprox, npe):
    # gaprox and npe for the elements, not the particles
    FFORMA, _, _, _, _ = shapeFunc(gaprox, npe)
    fi, _, _ = pascal_triangle(FFORMA, gaprox)
    UPART = []  # Initialize a list to store UPART_i
    SIGPART=[]  # Initialize a list to store SIGMA_i

    #=REMOVING LAST COLUMN OF SIGMA AND TRANSFORMING INTO A VECTOR=#
    SIGMA_M=np.array(SIGMA_full[:,:-1])
    SIGMA_V=SIGMA_M.flatten()
    for part in PARTICLES:
        U_elem = []  # Stores the displacement vector for each element
        SIG_elem=[]  # Stores the stress vector for each element
        FI_nod = np.zeros((3 * 2, 3 * 2 * npe))  # Shape (6, 6*npe)
        FI_sigma=np.zeros((3 * 3 , 3 * 3 * npe))

        for nod in range(3):  # Loop through 3 nodes for each particle
            elem = int(KSIETA[int(part[nod] - 1)][0]-1)  # Get the corresponding element

            # Create vector with displacements U and V for each node of the element
            for node in range(npe):
                current_node = int(ELEMS[elem][node] - 1) # 0-based index
                U_elem.extend(U[2*current_node : 2*current_node + 2])
                SIG_elem.extend(SIGMA_V[3*current_node : 3*current_node + 3])

            # Get ksi and eta coordinates
            ksi = KSIETA[int(part[nod] - 1)][1]
            eta = KSIETA[int(part[nod] - 1)][2]
  
            #__ Fill the FI matrix for calculating U and SIGMA __#
            
            for dir in range(2):
                line = nod * 2 + dir
                for n in range(npe):
                    column = nod * 2 * npe + 2 * n + dir
                    FI_nod[line][column] = fi(ksi, eta)[n]  # Use shape function value

            for component in range(3):
                line=3*nod+component
                for n in range(npe):
                    column=nod*3*npe+3*n+component
                    FI_sigma[line][column]=fi(ksi,eta)[n]


        # Compute UPART_i = FI * U_elem
        U_elem = np.array(U_elem)  
        UPART_i = np.dot(FI_nod, U_elem)

        UPART.append(UPART_i)  


        # exit()

        SIGPART_i=np.dot(FI_sigma,SIG_elem)
        SIGPART.append(SIGPART_i)

    # After the loop, UPART will be a list of arrays. Convert it to a NumPy array.
    U_PART = np.array(UPART)
    SIGMA_P=np.array(SIGPART)
    # Initialize an empty list for reshaped UPART
    # Reshape and accumulate UPART in a single step


    return U_PART, SIGMA_P

