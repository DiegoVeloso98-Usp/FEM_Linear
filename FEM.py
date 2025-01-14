import numpy as np
import importlib

nph=int(12)
hammer = np.array([
   [0.501426509658179, 0.249286745170910, 0.249286745170910, 
    0.873821971016996, 0.063089014491502, 0.063089014491502, 
    0.053145049844816, 0.310352451033785, 0.636502499121399, 
    0.310352451033785, 0.636502499121399, 0.053145049844816
    ],
   [0.249286745170910, 0.249286745170910, 0.501426509658179, 
    0.063089014491502, 0.063089014491502, 0.873821971016996, 
    0.310352451033785, 0.636502499121399, 0.053145049844816, 
    0.053145049844816, 0.310352451033785, 0.636502499121399
    ],
   [0.116786275726379/2.0, 0.116786275726379/2.0, 
    0.116786275726379/2.0, 0.050844906370207/2.0, 
    0.050844906370207/2.0, 0.050844906370207/2.0, 
    0.082851075618374/2.0, 0.082851075618374/2.0, 
    0.082851075618374/2.0, 0.082851075618374/2.0, 
    0.082851075618374/2.0, 0.082851075618374/2.0
    ]
   ])

        ###################################
        #### READ DATA FROM INPUT FILE ####
        ###################################

file_name = 'tensao.txt'  
with open(file_name, 'r') as file:
    dados = file.readlines()
dados = [line.strip().split() for line in dados]
count=0
nnodes=int(dados[count][1])
count+=1
nmats=int(dados[count][1])
count+=1
nthick=int(dados[count][1])
count+=1
nelem3=int(dados[count][1])
count+=1
nelem6=int(dados[count][1])
count+=1
nelem10=int(dados[count][1])

nelems=nelem3+nelem6+nelem10

count+=1
nforces=int(dados[count][1])
count+=1
npressures=int(dados[count][1])
count+=1
ndisplas=int(dados[count][1])
count+=3

NODES=np.zeros((nnodes,2))

for i in range(nnodes):
    NODES[i][0]=dados[count][1]
    NODES[i][1]=dados[count][2]
    count+=1
count+=2
MATS=np.zeros((nmats,2))            # column 1: E, column 2: ni

#################################
### PROPERTIES MATRIX (E, ni) ###
#################################

for i in range(nmats):              ### dados[count] ainda tem uma terceira coluna, de densidade ###
    MATS[i][0]=dados[count][1]
    MATS[i][1]=dados[count][2]
    count+=1

#################################
### THICKNESS VECTOR (h, mat) ###
#################################

count+=2
THICKS=np.zeros(int(nthick))     

for i in range(nthick):
    THICKS[i]=dados[count][1]
    count+=1

#################################
#### ELEMENTS MATRIX (nodes) ####
#################################
count+=2
if nelem3!=0:
    npe=0           ## ACREDITO QUE SEJA 3, MAS NO MAT ESTÁ 0 !!!!!!!
elif nelem6!=0:
    npe=6
else:
    npe=10
print(npe)
ELEMS=np.zeros((nelems,npe+5))      
for iel in range(nelems):                   # Lopping over the elements
    for ino in range(npe):                  # Looping over the nodes of the element 
        ELEMS[iel][ino]=dados[count][ino+1] # (Last 5 columns are the material properties)
    
    ## Filling the last 5 columns with the material properties ##
    for i in range(5):
        if i<2:
            ELEMS[iel][npe+i]=MATS[int(dados[count][npe+i+1])-1][i]
        elif i==2:
            ELEMS[iel][npe+i]=THICKS[int(dados[count][npe+i])-1]   # CONFERIR INDECES AQUI
        else:
            ELEMS[iel][npe+i]=0
    count+=1
count+=2

#######################################
#### ESSA PARTE AQUI ESTA ESTRANHA ####
#######################################
LOAD = []

for ino in range(nforces):
    # If the third column of 'dados' is non-zero, append specific load information
    if dados[count][2] != 0:  # Adjust for 0-based indexing in Python
        LOAD.append([dados[count][1], 1, dados[count][2]])
    # Otherwise, if the element in 'dados[count]' is non-zero, append a different set of load data
    if dados[count][3] != 0:  # Adjust for 0-based indexing
        LOAD.append([dados[count][1], 2, dados[count][3]])
    count+=1

#######################################
## ESSA PARTE AQUI TBM ESTA ESTRANHA ##
#######################################

count+=2
VINC=[]

for ino in range(ndisplas):
    if dados[count][2]=="'X'":
        VINC.append([int(dados[count][1]),1,float(dados[count][3]),1])
    if dados[count][2]=="'Y'":
        VINC.append([int(dados[count][1]),2,float(dados[count][3]),1])
    if dados[count][2]=="'BOTH'":
        VINC.append([int(dados[count][1]),1,float(dados[count][3]),1])
        VINC.append([int(dados[count][1]),2,float(dados[count][3]),1])
    count+=1

#######################################
#######################################

gaprox=int(3)
ep=0

#### COORDINATES KSI AND ETA ####

COORD = np.empty((0, 2))
for il in range(1,(gaprox+1)+1):                            # loop over the lines of the element 
    for ic in range(1,(gaprox+1-il+1)+1):                   # loop over the columns of the element
        ksi=1/gaprox*(ic-1)                                 # ksi coordinate (x) of the node of the element
        eta=1/gaprox*(il-1)                                 # eta coordinate (y) of the node of the element
        COORD = np.append(COORD, [[ksi, eta]], axis=0)      # append the coordinates to the array coord as a row

#### COEFFICIENTS AND DERIVATIVES ####

MAT=np.zeros((npe, npe))
counter=0
COEF=[]
DCOEFDSKI=[]
DCOEFDETA=[]
for il in range(1,(gaprox+1)+1):                            # loop over the lines of the TRIANGULAR element
    for ic in range(1,(gaprox+1-il+1)+1):                   # loop over the columns of the TRIANGULAR element
        COEF.append(f'ksi**{ic-1}*eta**{il-1}')
        DCOEFDSKI.append(f'{ic-1}*ksi**{ic-2}*eta**{il-1}')
        DCOEFDETA.append(f'ksi**{ic-1}*{il-1}*eta**{il-2}')
        counter+=1
        for ino in range(int(npe)):                         # For gaprox=1, the nested for's il and ic happens 3 times, so ino goes through all 3x3 positions
            ksi=COORD[ino][0]
            eta=COORD[ino][1]
            termo1=1 if ic==1 else ksi**(ic-1)              ### NÃO ENTENDI ESSA PARTE
            termo2=1 if il==1 else eta**(il-1)
            termo=termo1*termo2
            MAT[[ino],[counter-1]]=termo
#### SHAPE FUNCTIONS ####

FFORMA=np.empty((0, npe))

for if_ in range(npe):
    vec=np.zeros(npe)
    vec[if_]=1
    FFORMA = np.append(FFORMA, [np.linalg.solve(MAT, vec)], axis=0)

fi = lambda KSI, ETA: np.dot(FFORMA,[eval(c)for c in COEF])
dfidski = lambda KSI, ETA: np.dot(FFORMA,[eval(d)for d in DCOEFDSKI])
dfideta = lambda KSI, ETA: np.dot(FFORMA,[eval(d)for d in DCOEFDETA])

FGLOBAL=np.zeros(2*nnodes)
KGLOBAL=np.zeros((2*nnodes,2*nnodes))
#### MATERIAL PROPERTIES ####
for iel in range(int(nelems)):
    cx=[]
    cy=[]
    for ine in range(int(npe)):
        cx.append(NODES[int(ELEMS[iel][ine])-1][0])
        cy.append(NODES[int(ELEMS[iel][ine])-1][1])
    for ih in range(nph):
        ksi=hammer[0,ih]
        eta=hammer[1,ih]
        peso=hammer[2,ih]

        dxdksi=np.dot(dfidski(ksi,eta),cx)
        dxdeta=np.dot(dfideta(ksi,eta),cx)
        dydksi=np.dot(dfidski(ksi,eta),cy)
        dydeta=np.dot(dfideta(ksi,eta),cy)
        vfi=fi(ksi,eta)

        J=np.array([[dxdksi, dxdeta],[dydksi, dydeta]])

        detjac=np.linalg.det(J)
        JINV=np.linalg.inv(J)
        dksidx=JINV[0,0]
        dksidy=JINV[0,1]
        detadx=JINV[1,0]
        detady=JINV[1,1]

        DX=np.array([[dksidx,0,detadx,0],[0,dksidx,0,detadx]])
        DY=np.array([[dksidy,0,detady,0],[0,dksidy,0,detady]])

        DFI=np.zeros((4,2*npe))
        
        for ine in range(npe):            #AQUI O ERRO CABAÇO
            DFI[0,2*(ine-1)+1]=dfidski(ksi,eta)[ine]
            DFI[1,2*(ine-1)+2]=dfidski(ksi,eta)[ine]
            DFI[2,2*(ine-1)+1]=dfideta(ksi,eta)[ine]
            DFI[3,2*(ine-1)+2]=dfideta(ksi,eta)[ine]

        young=ELEMS[iel][npe]
        ni=ELEMS[iel][npe+1]
        h=ELEMS[iel][npe+2]
        bx=ELEMS[iel][npe+3]        #Body forces bx and by
        by=ELEMS[iel][npe+4]        # !!! ELEMS[iel][4 & 5] are being assigned 0 !!!

        ## TRANSFORMING BODY FORCES INTO NODAL FORCES ##
        FLOCAL=np.zeros(2*npe)
        for ine in range(npe):
            FLOCAL[2*ine]=h*bx*detjac*peso*vfi[ine]
            FLOCAL[2*ine+1]=h*by*detjac*peso*vfi[ine]

        #################################################

        if ep==0:
            d11=young/(1-ni**2)
            d12=d11*ni
            d21=d12
            d22=d11
            d33=d11*(1-ni)
            d44=d33
        else:
            d11=young/(1+ni)/(1-2*ni)*(1-ni)
            d12=d11*ni
            d21=d12
            d22=d11
            d33=d11*(1-2*ni)
            d44=d33
        MXX=np.array([[d11,0],[0,d33]])
        MXY=np.array([[0,d12],[d33,0]])
        MYX=np.array([[0,d44],[d21,0]])
        MYY=np.array([[d44,0],[0,d22]])

        KLOCAL=h*DFI.T @(DX.T @ MXX@DX +DX.T@MXY@DY+DY.T@MYX@DX+DY.T@MYY@DY)@DFI*detjac*peso

        ############################################################
        ##### GLOBAL FORCES VECTOR AND GLOBAL STIFFNESS MATRIX #####
        ############################################################

        for in_ in range(npe):
            for idir in range(2):
                FGLOBAL[2*int(ELEMS[iel][in_]-1)+idir]=+FLOCAL[2*(ine)+idir]
                for jn in range(npe):
                    for jdir in range(2):
                        KGLOBAL[2*int(ELEMS[iel][in_]-1)+idir,2*int(ELEMS[iel][jn]-1)+jdir]+=KLOCAL[2*(in_-1)+idir,2*(jn-1)+jdir]

for i in range(len(LOAD)):
    ino=int(LOAD[i][0])
    idir=int(LOAD[i][1])
    value=float(LOAD[i][2])
    FGLOBAL[2*(ino-1)+idir-1]+=value  

for i in range(len(VINC)):
    ino=int(VINC[i][0])
    idir=int(VINC[i][1])
    value=float(VINC[i][2])
    type_=int(VINC[i][3])
    if type_==1:
        FGLOBAL-=value*KGLOBAL[:,2*(ino-1)+idir-1]      # FGLOBAL=FGLOBAL-value*KGLOBAL[columns of the fixed node]
        KGLOBAL[:,2*(ino-1)+idir-1]=0                   # Nulling the columns of the fixed node
        KGLOBAL[2*(ino-1)+idir-1,:]=0                   # Nulling the rows of the fixed node
        KGLOBAL[2*(ino-1)+idir-1,2*(ino-1)+idir-1]=1    # Assigning 1 to the diagonal of the fixed node
        FGLOBAL[2*(ino-1)+idir-1]=value                 
    else:
        KGLOBAL[2*(ino-1)+idir-1,2*(ino-1)+idir-1]+=value

U=np.linalg.solve(KGLOBAL, FGLOBAL)
u=np.copy(U)

displa = np.copy(NODES)
for in_node in range(nnodes):
    displa[in_node, 0] += 100 * u[2 * in_node]
    displa[in_node, 1] += 100 * u[2 * in_node + 1]

# Plot the deformed structure
import matplotlib.pyplot as plt
plt.fill(displa[:, 0], displa[:, 1], 'b', alpha=0.5)  # Fill the area between the points with blue color and 50% transparency
plt.scatter(displa[:, 0], displa[:, 1], s=10)  # Set the size of the points to 10
plt.show()




