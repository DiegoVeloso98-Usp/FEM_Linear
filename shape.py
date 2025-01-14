import numpy as np

import hammerpoints

#========================================#
#========================================#

def pascalTriangleCoef(gaprox):
    COEF=[]
    DCOEFDSKI=[]
    DCOEFDETA=[]
    for il in range(1, (gaprox + 1) + 1):  # loop over the lines of the TRIANGULAR element
        for ic in range(1, (gaprox + 1 - il + 1) + 1):  # loop over the columns of the TRIANGULAR element

            # Avoid division by zero when computing powers of ksi and eta
            if ic > 1:
                dcoef_dksi_term = f'{ic-1}*ksi**{ic-2}*eta**{il-1}'
            else:
                dcoef_dksi_term = f'0*eta**{il-1}'  # when ic == 1, derivative is 0 w.r.t ksi
            
            if il > 1:
                dcoef_deta_term = f'ksi**{ic-1}*{il-1}*eta**{il-2}'
            else:
                dcoef_deta_term = f'ksi**{ic-1}*0'  # when il == 1, derivative is 0 w.r.t eta
            
            # Append the coefficients
            COEF.append(f'ksi**{ic-1}*eta**{il-1}')
            DCOEFDSKI.append(dcoef_dksi_term)
            DCOEFDETA.append(dcoef_deta_term)
    return COEF, DCOEFDSKI, DCOEFDETA

#========================================#
#========================================#

def pascal_triangle(FFORMA, gaprox):

    COEF, DCOEFDSKI, DCOEFDETA = pascalTriangleCoef(gaprox)

    fi = lambda KSI, ETA: np.dot(FFORMA,[eval(c, {'ksi': KSI, 'eta': ETA}) for c in COEF])
    dfidski = lambda KSI, ETA: np.dot(FFORMA, [eval(d, {'ksi': KSI, 'eta': ETA}) for d in DCOEFDSKI])
    dfideta = lambda KSI, ETA: np.dot(FFORMA, [eval(d, {'ksi': KSI, 'eta': ETA}) for d in DCOEFDETA])

    return fi, dfidski, dfideta

#========================================#
#========================================#

def shapeFunc(gaprox,npe):
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

    for il in range(1,(gaprox+1)+1):                            # loop over the lines of the TRIANGULAR element
        for ic in range(1,(gaprox+1-il+1)+1):                   # loop over the columns of the TRIANGULAR element
            COEF, DCOEFDSKI, DCOEFDETA = pascalTriangleCoef(gaprox)

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
    print("MAT:")
    print(MAT)
    # exit()
    for if_ in range(npe):
        vec=np.zeros(npe)
        vec[if_]=1
        FFORMA = np.append(FFORMA, [np.linalg.solve(MAT, vec)], axis=0)
    print("FFORMA:")
    print(FFORMA)
    # exit()

    return FFORMA, COEF, DCOEFDSKI, DCOEFDETA, COORD

#========================================#
#========================================#

if __name__ == '__main__':
    gaprox = 3
    npe = 10
    FFORMA, COEF, DCOEFDSKI, DCOEFDETA, COORD = shapeFunc(gaprox, npe)
    print(FFORMA)
    print(COEF)
    print(DCOEFDSKI)
    print(DCOEFDETA)
    print(COORD)
    nph,hammer=hammerpoints.hammer()
    for ih in range(nph):
        fi, dfidski, dfideta = pascal_triangle(FFORMA, gaprox)
        ksi=hammer[0,ih]
        eta=hammer[1,ih]
        # print(f'fi({ksi},{eta}):', fi(ksi, eta))
        # print(f'dfidski({ksi},{eta}):', dfidski(ksi, eta))

    ksi=0.501426509658179
    eta=0.249286745170910
    print("COORD:", COORD)
    print('FFORMA:', FFORMA)
    print(f'fi({ksi},{eta}):', fi(ksi, eta))
    print(f'dfidski({ksi},{eta}):', dfidski(ksi, eta))
    print(f'dfideta({ksi},{eta}):', dfideta(ksi, eta))

    print('Done')
