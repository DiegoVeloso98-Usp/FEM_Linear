
import numpy as np

def generate_particles(NODES,radius,fv,young,ni,h,bx,by):
 
    xi=min(NODES[:,0])
    yi=min(NODES[:,1])

    xf=max(NODES[:,0])
    yf=max(NODES[:,1])

    area=(xf-xi)*(yf-yi)

    xi+=radius
    yi+=radius

    xf-=radius
    yf-=radius

    dx,dy=[(xf-xi),(yf-yi)]
    areap=3*(radius**2)*(3**0.5)/4
    fraction=0
    NODESP=[]
    PARTICLES=[]
    while fraction<fv*area:
        xc=xi+dx*np.random.rand()
        yc=yi+dy*np.random.rand()
        alfa=(2*np.pi/3)*np.random.rand()
        
        NODESP.append([xc+radius*np.cos(alfa),yc+radius*np.sin(alfa)])
        NODESP.append([xc+radius*np.cos(alfa+2*np.pi/3),yc+radius*np.sin(alfa+2*np.pi/3)])
        NODESP.append([xc+radius*np.cos(alfa+4*np.pi/3),yc+radius*np.sin(alfa+4*np.pi/3)])
        PARTICLES.append([len(NODESP)-2,len(NODESP)-1,len(NODESP),young,ni,h,bx,by])

        fraction+=areap
    NODESP=np.array(NODESP)
    PARTICLES=np.array(PARTICLES)
    return NODESP, PARTICLES


def particle_location(NODES,ELEMS, npe, gaprox,radius,fv,young,ni,h,bx,by):  

    NODESP,PARTICLES=generate_particles(NODES,radius,fv,young,ni,h,bx,by)

    ksieta=[]

    for nodes in NODESP:
        for i,elem in enumerate(ELEMS):

            xp=nodes[0]
            yp=nodes[1]
            ##### Location of the vertices nodes of the element #####
            x1,y1=NODES[int(elem[0]-1)][0]  ,  NODES[int(elem[0]-1)][1]
            
            x2,y2=NODES[int(elem[gaprox]-1)][0]  ,  NODES[int(elem[gaprox]-1)][1]
            
            x3,y3=NODES[int(elem[npe-1]-1)][0]  ,  NODES[int(elem[npe-1]-1)][1]
            
            ksi=-((x3*y1-xp*y1-x1*y3+xp*y3+x1*yp-x3*yp)/(x2*y1-x3*y1-x1*y2+x3*y2+x1*y3-x2*y3))
            eta=-((x2*y1-xp*y1-x1*y2+xp*y2+x1*yp-x2*yp)/(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3))

            if ksi>=0 and ksi<=1 and eta>=0 and eta<=1 and (1-ksi-eta)>=0 and (1-ksi-eta)<=1:
                ksieta.append([i+1,ksi,eta])
                break

    ksieta=np.array(ksieta)
    #====   N_PARTICLES = N_NODES/3 (3 NODES PER PARTICLE)  ====#
    return NODESP, PARTICLES, len(NODESP), int(len(NODESP)/3),ksieta 
'''
class ParticleGenerator:
    def __init__(self, filename):
        self.filename = filename
        self.reading = read(filename)
        self.particles = self.reading.get_particles()

    def get_particles(self):
        return self.particles
'''