import os
import time

import numpy as np

import assembly
import hammerpoints
import particle_generator
import plotting
import properties
import reading
import shape
import sigma
import solver
import writing_acadview
import writing_paraview

#=============================================================================================================================#
#========================================================= LINEAR FEM ========================================================#
#=============================================== by *Diego Dias Veloso*, 2024 ================================================#
#=============================================================================================================================#


def main():
    start_time = time.time()

    #====================================================#
    #================= READ INPUT FILE ==================#
    #====================================================#

    title='tensao.txt'
    input=reading.read(title)
    dados,nnodes,nmats,nthick, nelem3, nelem6, nelem10,nelems,nforces,npressures,ndisplas,count=input

    #====================================================#
    #========= EXTRACT PROPERTIES FROM RAW DATA =========#
    #====================================================#
    NODES,MATS,THICKS,ELEMS,LOAD, VINC, gaprox, npe, count=properties.inputProperties(*input)
    nph, hammer=hammerpoints.hammer()
    ep=0

    #=====================================================#
    #============ CALCULATE SHAPE FUNCTIONS ==============#
    #=====================================================#

    FFORMA, COEF, DCOEFDSKI, DCOEFDETA, COORD=shape.shapeFunc(gaprox,npe)

    #=====================================================#
    #================= GENERATE PARTICLES ================#
    #=====================================================#

    fv=0.0001
    radius=0.1
    npe_p=3
    young_p,ni_p,h_p,bx_p,by_p=1000,0.25,1,0,0
    NODESP, PARTICLES, npoints, nparticles,KSIETA=particle_generator.particle_location(NODES,ELEMS, npe, gaprox,radius,fv,young_p,ni_p,h_p,bx_p,by_p)


    # NODESP=[[0.89346,0.171142],
    #         [0.754277, 0.06805],
    #         [0.91314,-0.0001]
    #         ]
    # PARTICLES=[[1,2,3,young_p,ni_p,h_p,bx_p,by_p]]
    # KSIETA=[[2,0.12232,0.77114],[2,0.08622,0.66806],[2,0.3141,0.59905]]
    # nparticles=1
    #=====================================================#
    #= ASSEMBLY GLOBAL STIFFNESS MATRIX AND FORCE VECTOR =#
    #=====================================================#

        #= Elements global stiffness matrix =#
    KGLOBAL,FGLOBAL, dfidski, dfideta=assembly.globalStiffness(ELEMS, NODES, LOAD, VINC, FFORMA, gaprox,npe, nelems, nnodes, nph, hammer, ep)
  
        #= Elements and particles Global stiffness matrix =#
    gaprox_p,npe_p=1,3
    FFORMA_P, COEF_P, DCOEFDSKI_P, DCOEFDETA_P, COORD_P=shape.shapeFunc(gaprox_p,npe_p)
    # KGLOBAL, FGLOBAL, fi,dfideta, dfidski=assembly.globalParticleStiffness(KGLOBAL,FGLOBAL,ELEMS,PARTICLES, NODESP, KSIETA, LOAD, VINC,FFORMA, FFORMA_P,gaprox, gaprox_p,npe,npe_p, nparticles, nph, hammer, ep)
    
 
    #=====================================================#
    #======== SOLVE FOR THE UNKNOWN DISPLACEMENTS ========#
    #=====================================================#
    
    print('FGLOBAL:',FGLOBAL)
    print('KGLOBAL:',KGLOBAL)
    U=solver.fem(KGLOBAL,FGLOBAL)
    # print('U:',U)
    # exit()

    #=====================================================#
    #================= CALCULATE STRESSES ================#
    #=====================================================#

    SIGMA=sigma.sigma(ELEMS, NODES, U, nnodes, nelems, npe, ep, gaprox)


    U_PART,SIGMA_P_m=solver.particleDisplacement(U,SIGMA,PARTICLES,KSIETA,ELEMS,gaprox,npe)


    # Transforming matrix SIGMA_P_m (9 columns) into a matrix SIGMA_P (3 columns)
    SIGMA_P = SIGMA_P_m.reshape(-1, 3)


    end_time = time.time()
    print(f"Execution time: {end_time - start_time} seconds")

#=============================================================================================================================#
#======================================================= POST-PROCESSING =====================================================#
#=============================================================================================================================#

    
    #=====================================================#
    #===================== PLOT RESULTS ==================#
    #=====================================================#
    # Group pairs using zip and slicing

    a=10   #- Amplification factor -#

    plotting.undeformedPlotting(ELEMS,NODES)
    plotting.deformedPlotting(ELEMS,NODES,U,nnodes,a)


    # #_Transforming vectors U
    # UPART_matrix=list(zip(U_PART[::2], U_PART[1::2]))
    # U_matrix = list(zip(U[::2], U[1::2]))
    # plotting.displacedParticles(NODESP,U_PART,NODES,U_matrix,a)
    #=====================================================#
    #============ EXPORT RESULTS TO ACADVIEW =============#
    #=====================================================#

    # Get the directory of the current script
    folder = os.path.dirname(os.path.abspath(__file__))

    acadview_name = title.split('.')[0]

    posproc1=writing_acadview.appending(SIGMA,ELEMS,NODES,U,SIGMA_P,PARTICLES,NODESP,U_PART, nelems, nnodes,npe, gaprox)

    path = f'{folder}\\{acadview_name}.ogl'
    writing_acadview.exporting(posproc1, path)

    #=====================================================#
    #============ EXPORT RESULTS TO PARAVIEW =============#
    #======== (CURRENTLY ONLY FOR 3-NODED ELEMS) =========#
    #=====================================================#

    paraview_name = title.split('.')[0]
    for direction in ['sigmaX', 'sigmaY', 'sigmaXY']:
        path = f'{folder}\\{paraview_name}_{direction}.vtk.49'
        npe=3
        posproc2 = writing_paraview.appending_vtk_sigma(SIGMA, U, ELEMS, NODES, nnodes, nelems, title, direction, npe)
        writing_paraview.exporting_paraview(posproc2, path)
    
    #=========================================================#
    #=========================================================#
    #=========================================================#

if __name__ == '__main__':
    main()

#========================================================================================================================================#
#========================================================================================================================================#
#========================================================================================================================================#

