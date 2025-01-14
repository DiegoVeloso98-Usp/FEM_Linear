import numpy as np

from properties import materialProperties
from shape import pascal_triangle

#==========================================================================================================#
#==========================================================================================================#

def localStiffness(ELEMS,npe, hammer, ep,cx,cy,ih,iel,fi, dfideta,dfidski):
    # fi, dfidski, dfideta = pascal_triangle(FFORMA, gaprox)

    ksi=hammer[0,ih]
    eta=hammer[1,ih]
    peso=hammer[2,ih]

    dxdksi=np.dot(dfidski(ksi,eta),cx)
    dxdeta=np.dot(dfideta(ksi,eta),cx)
    dydksi=np.dot(dfidski(ksi,eta),cy)
    dydeta=np.dot(dfideta(ksi,eta),cy)
    vfi=fi(ksi,eta)
    # print(f'ksi: {ksi}, eta: {eta}')
    # print(f'fi,{fi(ksi,eta)} ')
    # print(f' dfidski: {dfidski(ksi,eta)}')
    # print(f'dfideta: {dfideta(ksi,eta)}')



    J=np.array([[dxdksi, dxdeta],[dydksi, dydeta]])

    # print('J:',J)
    # exit()

    detjac=np.linalg.det(J)

    # print('detjac:',detjac)
    # exit()

    JINV=np.linalg.inv(J)

    # print('Jinv:',JINV)
    # exit()
    
    dksidx=JINV[0,0]
    dksidy=JINV[0,1]
    detadx=JINV[1,0]
    detady=JINV[1,1]

    DX=np.array([[dksidx,0,detadx,0],[0,dksidx,0,detadx]])
    DY=np.array([[dksidy,0,detady,0],[0,dksidy,0,detady]])

    # print('DX:',DX)
    # print('DY:',DY)
    # exit()
    DFI=np.zeros((4,2*npe))
    
    for ine in range(npe):            
        DFI[0,2*(ine)]=dfidski(ksi,eta)[ine]
        DFI[1,2*(ine)+1]=dfidski(ksi,eta)[ine]
        DFI[2,2*(ine)]=dfideta(ksi,eta)[ine]
        DFI[3,2*(ine)+1]=dfideta(ksi,eta)[ine]

    d11, d12, d21, d22, d33, d44, h, bx, by=materialProperties(ELEMS, iel, npe, ep)

    # print('DFI:',DFI)
    # exit()

    #== TRANSFORMING BODY FORCES INTO NODAL FORCES ==#
    FLOCAL=np.zeros(2*npe)
    for ine in range(npe):
        FLOCAL[2*ine]=h*bx*detjac*peso*vfi[ine]
        FLOCAL[2*ine+1]=h*by*detjac*peso*vfi[ine]

    #================================================#

    MXX=np.array([[d11,0],[0,d33]])
    MXY=np.array([[0,d12],[d33,0]])
    MYX=np.array([[0,d44],[d21,0]])
    MYY=np.array([[d44,0],[0,d22]])

    # print('MXX:',MXX)
    # print('MXY:',MXY)
    # print('MYX:',MYX)
    # print('MYY:',MYY)
    # exit()

    KLOCAL=h*DFI.T @(DX.T @ MXX@DX +DX.T@MXY@DY+DY.T@MYX@DX+DY.T@MYY@DY)@DFI*detjac*peso

    return KLOCAL, FLOCAL


#==========================================================================================================#
#==========================================================================================================#

def globalStiffness(ELEMS, NODES, LOAD, VINC, FFORMA, gaprox,npe, nelems, nnodes, nph, hammer, ep):
    FGLOBAL=np.zeros(2*nnodes)
    KGLOBAL=np.zeros((2*nnodes,2*nnodes))
    #-Since shape functions dont change, I dont calculate them in the loop-#
    fi, dfidski, dfideta = pascal_triangle(FFORMA, gaprox)
    for iel in range(int(nelems)):
        print(f'Elemento {iel}')
        cx=[]
        cy=[]
        for ine in range(int(npe)):
            cx.append(NODES[int(ELEMS[iel][ine])-1][0])
            cy.append(NODES[int(ELEMS[iel][ine])-1][1])
        # print('cx:',cx)
        # print('cy:',cy)
        # exit()
        klocal_test=np.zeros((2*npe,2*npe))
        for ih in range(nph):
            # print(f'CX: {cx}')
            KLOCAL,FLOCAL=localStiffness(ELEMS,npe, hammer, ep,cx,cy,ih,iel,fi, dfideta,dfidski)  
            # print('KLOCAL:',KLOCAL)
            # exit()
            klocal_test=klocal_test+KLOCAL
                       
            for in_ in range(npe):
                for idir in range(2):
                    FGLOBAL[2*int(ELEMS[iel][in_]-1)+idir]=+FLOCAL[2*(in_)+idir]
                    for jn in range(npe):
                        for jdir in range(2):
                            KGLOBAL[2*int(ELEMS[iel][in_]-1)+idir,2*int(ELEMS[iel][jn]-1)+jdir]+=KLOCAL[2*(in_)+idir,2*(jn)+jdir]
                            #print('sum KGLOBAL',sum(sum(row) for row in KGLOBAL))
                            #print('KLOCAL',KLOCAL)
                            #print('sum KLOCAL',sum(sum(row) for row in KLOCAL))
                            #exit()

        # print('KLOCAL:',klocal_test)
        # exit()
    FGLOBAL=globalForce(LOAD, FGLOBAL)
    KGLOBAL, FGLOBAL=boundaryCondition(VINC, KGLOBAL, FGLOBAL)
    return KGLOBAL, FGLOBAL, dfideta, dfidski


def globalParticleStiffness(KGLOBAL,FGLOBAL,ELEMS,PARTICLES, NODESP, KSIETA, LOAD, VINC, FFORMA,FFORMA_P,gaprox, gaprox_p,npe,npe_p, npart, nph, hammer, ep):
    fi_p,dfidski_p, dfideta_p=pascal_triangle(FFORMA_P, gaprox_p)
    fi, dfidski, dfideta = pascal_triangle(FFORMA, gaprox)
    for iel in range(int(npart)):
        print(f'Particle {iel}')
        cx=[]
        cy=[]
        for ine in range(int(npe_p)):
            cx.append(NODESP[int(PARTICLES[iel][ine])-1][0])
            cy.append(NODESP[int(PARTICLES[iel][ine])-1][1])
        for ih in range(nph):
            KLOCAL,FLOCAL=localStiffness(PARTICLES,npe_p, hammer, ep,cx,cy,ih,iel,fi_p, dfideta_p,dfidski_p)

            #=== VFI_node_i = FI ( KSI_node_i, ETA_node_i ) ===#


            vfi1=fi(KSIETA[int(PARTICLES[iel][0])-1][1],KSIETA[int(PARTICLES[iel][0])-1][2])
            vfi2=fi(KSIETA[int(PARTICLES[iel][1])-1][1],KSIETA[int(PARTICLES[iel][1])-1][2])
            vfi3=fi(KSIETA[int(PARTICLES[iel][2])-1][1],KSIETA[int(PARTICLES[iel][2])-1][2])
            MFI=np.zeros((6,3*2*npe))
            indexes=[]
            for i in range(npe_p):
                for j in range(npe):
                    for  k in range(2):

                        #                   {node     [ elemento (           node         )        ]}
                        indexes.append(2*int(ELEMS[int(KSIETA[int(PARTICLES[iel][i]-1)][0])-1][j]-1)+k)

            for i in range(npe):
                MFI[0][2*i]=vfi1[i]
                MFI[1][2*i+1]=vfi1[i]
                MFI[2][2*npe+2*i]=vfi2[i]
                MFI[3][2*npe+2*i+1]=vfi2[i]
                MFI[4][4*npe+2*i]=vfi3[i]
                MFI[5][4*npe+2*i+1]=vfi3[i]
            KG=MFI.T@KLOCAL@MFI
            print('indexes:',indexes)
            for i in range(npe_p*2*npe):
                for j in range(npe_p*2*npe):
                    KGLOBAL[indexes[i]][indexes[j]]+=KG[i][j]
                    
    # FGLOBAL=globalForce(LOAD, FGLOBAL)
    KGLOBAL, FGLOBAL=boundaryCondition(VINC, KGLOBAL, FGLOBAL)
    return KGLOBAL, FGLOBAL, fi,dfideta, dfidski


def globalForce(LOAD, FGLOBAL):
    for load in LOAD:
        ino=int(load[0])
        idir=int(load[1])
        value=float(load[2])
        FGLOBAL[2*(ino-1)+idir-1]+=value  
    return FGLOBAL

def boundaryCondition(VINC, KGLOBAL, FGLOBAL):
    for vinc in VINC:
        ino=int(vinc[0])
        idir=int(vinc[1])
        value=float(vinc[2])
        type_=int(vinc[3])
        if type_==1:
            FGLOBAL-=value*KGLOBAL[:,2*(ino-1)+idir-1]      # FGLOBAL=FGLOBAL-value*KGLOBAL[columns of the fixed node]
            KGLOBAL[:,2*(ino-1)+idir-1]=0                   # Nulling the columns of the fixed node
            KGLOBAL[2*(ino-1)+idir-1,:]=0                   # Nulling the rows of the fixed node
            KGLOBAL[2*(ino-1)+idir-1,2*(ino-1)+idir-1]=1    # Assigning 1 to the diagonal of the fixed node
            FGLOBAL[2*(ino-1)+idir-1]=value                 
        else:
            KGLOBAL[2*(ino-1)+idir-1,2*(ino-1)+idir-1]+=value

    return KGLOBAL, FGLOBAL