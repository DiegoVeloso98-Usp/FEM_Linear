import numpy as np
from shape import pascal_triangle
from shape import shapeFunc
from properties import materialProperties

#========================================#
#========================================#

def sigma(ELEMS, NODES, U, nnodes, nelems, npe, ep, gaprox):

    FFORMA, COEF, DCOEFDSKI, DCOEFDETA, COORD=shapeFunc(gaprox,npe)

    fi, dfidski, dfideta = pascal_triangle(FFORMA, gaprox)

    SIGMA=np.zeros((nnodes,4))

    for iel in range(nelems):
        cx,cy,du,dv=[],[],[],[]

        for ine in range(npe):
            cx.append(NODES[int(ELEMS[iel][ine]-1)][0])
            cy.append(NODES[int(ELEMS[iel][ine]-1)][1])
            du.append(U[2*int(ELEMS[iel][ine]-1)])
            dv.append(U[2*int(ELEMS[iel][ine]-1)+1])
 
        for ih in range(npe):
            ksi=COORD[ih,0]
            eta=COORD[ih,1]

            dxdksi=np.dot(dfidski(ksi,eta),cx)
            dxdeta=np.dot(dfideta(ksi,eta),cx)
            dydksi=np.dot(dfidski(ksi,eta),cy)
            dydeta=np.dot(dfideta(ksi,eta),cy)

            dudksi=np.dot(dfidski(ksi,eta),du)
            dudeta=np.dot(dfideta(ksi,eta),du)
            dvdksi=np.dot(dfidski(ksi,eta),dv)
            dvdeta=np.dot(dfideta(ksi,eta),dv)

            J=np.array([[dxdksi, dxdeta],[dydksi, dydeta]])

            # if np.linalg.det(J)>=0.75:
            #     exit(f'Volumetric deformation in element {iel+1} equal to {np.linalg.det(J)} is too high.')
                
            JINV=np.linalg.inv(J)

            dksidx=JINV[0,0]
            dksidy=JINV[0,1]
            detadx=JINV[1,0]
            detady=JINV[1,1]

            if (iel+1)%10==0:
                print("_______________")
                print("ELEM:",iel)
                print(f'J: {J}, JINV: {JINV}')
                print(f'dksidx: {dksidx}, dksidy: {dksidy}, detadx: {detadx}, detady: {detady}')
                print(f'dudksi: {dudksi}, dudeta: {dudeta}, dvdksi: {dvdksi}, dvdeta: {dvdeta}')
                print(f'epsonx: {epsonx}, epsony: {epsony}, epsonxy: {epsonxy}')
                print(f'sigmax: {sigmax}, sigmay: {sigmay}, talxy: {talxy}')
                # exit()
            d11, d12, d21, d22, d33, d44, h, bx, by=materialProperties(ELEMS, iel, npe, ep)

            #===============#
            #==== SIGMA ====#
            #===============#

            epsonx=dudksi*dksidx+dudeta*detadx
            epsony=dvdksi*dksidy+dvdeta*detady
            epsonxy=0.5*(dudksi*dksidy+dudeta*detady+dvdksi*dksidx+dvdeta*detadx)

            
            sigmax=d11*epsonx+d12*epsony
            sigmay=d21*epsonx+d22*epsony
            talxy=d33*epsonxy
            # print(f'sigmax: {sigmax}, sigmay: {sigmay}, talxy: {talxy}')
            # exit()
            SIGMA[int(ELEMS[iel][ih])-1,:]+=[sigmax,sigmay,talxy,1]
            
    return SIGMA

