import numpy as np


def inputProperties(dados,nnodes,nmats,nthick,nelem3, nelem6, nelem10,nelems, nforces,npressures,ndisplas,count):
    NODES=np.zeros((nnodes,2))

    for i in range(nnodes):
        NODES[i][0]=dados[count][1]
        NODES[i][1]=dados[count][2]
        count+=1
    count+=2

    #####################################################
    ############# PROPERTIES MATRIX (E, ni) #############
    #####################################################

    MATS=np.zeros((nmats,2))            # column 1: E, column 2: ni
    for i in range(nmats):              ### dados[count] ainda tem uma terceira coluna, de densidade ###
        MATS[i][0]=dados[count][1]
        MATS[i][1]=dados[count][2]
        count+=1

    #####################################################
    ############# THICKNESS VECTOR (h, mat) #############
    #####################################################

    count+=2
    THICKS=np.zeros(int(nthick))     

    for i in range(nthick):
        THICKS[i]=dados[count][1]
        count+=1

    ##########################################################
    ################## ELEMENTS MATRIX (nodes) ###############
    ##########################################################
    count+=2
    if nelem3!=0:
        npe=3           ## ACREDITO QUE SEJA 3, MAS NO MAT ESTÁ 0 !!!!!!!
        gaprox=1
    elif nelem6!=0:
        npe=6
        gaprox=2
    else:
        npe=10
        gaprox=3

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

    ###########################################################
    ###################### LOADS MATRIX #######################
    ###########################################################
    LOAD = []

    for ino in range(nforces):
        
        if dados[count][2] != 0: 
            LOAD.append([dados[count][1], 1, dados[count][2]])
        
        if dados[count][3] != 0:  
            LOAD.append([dados[count][1], 2, dados[count][3]])
        count+=1

    ############################################################
    ################### RESTRICTIONS MATRIX ####################
    ############################################################

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

    ################################################################
    ################################################################


    return NODES,MATS,THICKS,ELEMS,LOAD, VINC, gaprox, npe, count

def  materialProperties(ELEMS, elem, node, ep):
    #elem IS THE INDEX OF THE FOR THAT WILL GO THROUGH THE NUMBER nelems OF ELEMENTS (LINES OF ELEMS)
    #node IS THE INDEX OF THE FOR THAT WILL GO THROUGH THE NUMBER npe OF NODES(COLUMNS OF ELEMS)

    young=ELEMS[elem][node]
    ni=ELEMS[elem][node+1]
    h=ELEMS[elem][node+2]
    bx=ELEMS[elem][node+3]        #Body forces bx and by
    by=ELEMS[elem][node+4]        # !!!!!!! ELEMS[iel][4 & 5] are being assigned 0 !!!!!!!!
    if ep==0:
        d11=young/(1-ni**2)
        d12=d11*ni
        d21=d12
        d22=d11
        d33=d11*(1-ni)/2
        d44=d33
    else:
        d11=young/(1+ni)/(1-2*ni)*(1-ni)
        d12=d11*ni
        d21=d12
        d22=d11
        d33=d11*(1-ni)
        d44=d33

    return d11, d12, d21, d22, d33, d44, h, bx, by