import numpy as np
import os

NODES = np.array([[0,0],[1,0],[1,1],[0,1]])
ELEMS = np.array([[1,0,1,2],[2,0,2,3]])
U = np.array([0,0,0,0,0,0,0,0])
SIGMA = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

nelems = 2
nnodes = 4
npe = 3
gaprox = 1
def appending(SIGMA,ELEMS,NODES,U,SIGMA_P,PARTICLES,NODESP,U_PART, nelems, nnodes,npe, gaprox):
    posproc = []
    posproc.append('arquivo de entrada')
    posproc.append('n.nos n.elems n.listas')
    posproc.append('#')
    posproc.append(f'{nnodes+len(NODESP)} {nelems+len(PARTICLES)} 5')
    posproc.append('coordx coordy coordz deslx desly deslz')
    posproc.append('#')
    for i in NODES:
        posproc.append(f'{i[0]} {i[1]} {0} {0} {0} {0}')
    for i in NODESP:
        posproc.append(f'{i[0]} {i[1]} {0} {0} {0} {0}')
    
    
    posproc.append('tpelem (1-barra/2-triang/3-quad) grauaprox nó1 nó2...nó_n group \n')
    posproc.append('#')
    for i in ELEMS:
        i = i.tolist()  # Convert numpy array to Python list
        list=i[0:npe]
        list.append(0)
        list.insert(0, gaprox)
        list.insert(0,2)
        posproc.append(' '.join(map(str,map(int,list))))

    for i in PARTICLES:
        i = i.tolist()
        list = [x + nnodes for x in i[0:3]]
        list.append(1)
        list.insert(0, 1)
        list.insert(0, 2)
        posproc.append(' '.join(map(str,map(int,list))))

    posproc.append('listas \n')

    posproc.append('#')
    posproc.append('desl.x')
    for i in range(nnodes):
        posproc.append(f'{U[2*i]} {U[2*i+1]} {0} {U[2*i]}')
    print(U_PART)
    U_PART_V=np.ravel(U_PART)
    print(U_PART_V)
    # exit()

    for i in range(len(NODESP)):
        posproc.append(f'{U_PART_V[2*i]} {U_PART_V[2*i+1]} {0} {U_PART_V[2*i]}')
    posproc.append('#')
    posproc.append('desl.y')
    for i in range(nnodes):
        posproc.append(f'{U[2*i]} {U[2*i+1]} {0} {U[2*i+1]}')
    
    for i in range(len(NODESP)):
        posproc.append(f'{U_PART_V[2*i]} {U_PART_V[2*i+1]} {0} {U_PART_V[2*i+1]}')

    posproc.append('#')
    posproc.append('sigma.x')

    for i in range(nnodes):
        posproc.append(f'{U[2*i]} {U[2*i+1]} {0} {SIGMA[i][0]/SIGMA[i][3]}')

    for i in range(len(NODESP)):
        posproc.append(f'{U_PART_V[2*i]} {U_PART_V[2*i+1]} {0} {SIGMA_P[i][0]}')  

    posproc.append('#')
    posproc.append('sigma.y')
    for i in range(nnodes):
        posproc.append(f'{U[2*i]} {U[2*i+1]} {0} {SIGMA[i][1]/SIGMA[i][3]}')

    for i in range(len(NODESP)):
        posproc.append(f'{U_PART_V[2*i]} {U_PART_V[2*i+1]} {0} {SIGMA_P[i][1]}') 

    posproc.append('#')
    posproc.append('tal.xy')
    for i in range(nnodes):
        posproc.append(f'{U[2*i]} {U[2*i+1]} {0} {SIGMA[i][2]/SIGMA[i][3]}')
    
    for i in range(len(NODESP)):
        posproc.append(f'{U_PART_V[2*i]} {U_PART_V[2*i+1]} {0} {SIGMA_P[i][2]}')  

    return posproc
def exporting(posproc, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as file:
        file.write('\n'.join(map(str, posproc)))
    print(f"File '{path}' created successfully.")


if __name__ == '__main__':
    posproc=appending(SIGMA,ELEMS, NODES,U,nelems, nnodes, npe, gaprox)
    filename = 'G:\\Meu Drive\\Mestrado\\2º Semestre\\Estratégias de Programação\\Semana_3\\V.4\\posproc.ogl'
    exporting(posproc, filename)

