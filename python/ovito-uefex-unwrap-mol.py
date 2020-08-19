from ovito.data import *
import numpy as np

def get_wrapped_pos(upos,e0,e1,e2,N0,N1,N2,wpos):
    wpos[2] = upos[2] % e2[2]
    N2 = upos[2] / e2[2]
    uc=np.zeros(3)
    for i in range(3):
        uc[i]=upos[i]-N2*e2[i]

    wpos[1] = uc[1] % e1[1]
    N1 = uc[1] / e1[1]
    ucb=np.zeros(3)
    for i in range(3):
        ucb[i]=uc[i]-N1*e1[i]

    wpos[0] = ucb[0] % e0[0]
    N0 = ucb[0] / e0[0]

def color_rainbow(i,max):
    f=float(i)/float(max)
    r=0
    g=0
    b=0
    if f < 0.25:
        r=1.0
        g=4*f
        b=0.0
    elif f < 0.5:   
        r=1.0-4*(f-0.25)
        g=1.0
        b=0.0
    elif f < 0.75:
        r=0.0
        g=1.0
        b=4*(f-0.5)
    else:

        r=0.0
        g=1.0-4*(f-0.75)
        b=1.0
    return [r,g,b]

def vec_rot(rot,invec,outvec):
    outvec[0]=invec[0]*rot[0,0] +invec[1]*rot[1,0] +invec[2]*rot[2,0]
    outvec[1]=invec[0]*rot[0,1] +invec[1]*rot[1,1] +invec[2]*rot[2,1]
    outvec[2]=invec[0]*rot[0,2] +invec[1]*rot[1,2] +invec[2]*rot[2,2]


def modify(frame, data, odata):

    path='C:/Users/murasima/Desktop/'
    filename='rotation.txt'
    path_to_file=path+filename
    rotdata=[]
    with open(path_to_file,'r') as f:
        il=0
        for line in f:
            line2=line.split()
            il += 1
            if il == 1:
                print(line)
            elif il == 2:
                print(line)
            elif il > 2:
                rotdata.append(line2)
    print ("=== %s ===" % frame)
    iframe=int(frame)
    rot=np.zeros((3,3))
    rot[0,0]=float(rotdata[iframe][1])
    rot[0,1]=float(rotdata[iframe][2])
    rot[0,2]=float(rotdata[iframe][3])
    rot[1,0]=float(rotdata[iframe][4])
    rot[1,1]=float(rotdata[iframe][5])
    rot[1,2]=float(rotdata[iframe][6])
    rot[2,0]=float(rotdata[iframe][7])
    rot[2,1]=float(rotdata[iframe][8])
    rot[2,2]=float(rotdata[iframe][9])

    
    cell=data.cell
    ocell=odata.cell_

    pos=data.particles['Position']
    pbc=data.particles['Periodic Image']
    upos=odata.particles_['Position_']
    mol=odata.particles['Molecule Identifier']
    color=odata.particles_.create_property('Color')
    M=max(mol[:])
    NM=mol.size
    N=int(NM/M)
    print(N,M,NM)

    mol_pbc_max=np.zeros((M,3))
    mol_pbc_min=np.zeros((M,3))
    mol_pbc_diff=np.zeros((M,3))
    mol_pbc_flag=np.zeros(M)
    for i in range(M):
        for k in range(3):
            mol_pbc_max[i,k]=-10000
            mol_pbc_min[i,k]=10000
    image_pbc=np.zeros(3)
    for i in range(NM):
        imol=mol[i]-1
        for k in range(3):
            image_pbc[k]=pbc[i,k]
            if image_pbc[k] > mol_pbc_max[imol,k]:
                mol_pbc_max[imol,k]=image_pbc[k]
            if image_pbc[k] < mol_pbc_min[imol,k]:
                mol_pbc_min[imol,k]=image_pbc[k]
    for i in range(M):
        imol=mol[i]-1
        mol_pbc_diff[i]=mol_pbc_max[i]-mol_pbc_min[i]
        for k in range(3):
            if mol_pbc_diff[i,k]>10:
                mol_pbc_flag[i] += 1

    mol_center=np.zeros((M,3))
    mol_size=np.zeros(M)
    mol_center_wrapped=np.zeros((M,3))

    for i in range(NM):
        imol=mol[i]-1
        icolor=color_rainbow(imol,M-1)
        color[i]=icolor

    
    
    for k in range(3):
        print(cell[k,0],cell[k,1],cell[k,2])

    e0=cell[:,0]
    e1=cell[:,1]
    e2=cell[:,2]
    N0=0
    N1=0
    N2=0

    for k in range(3):
        print(e0[k],e1[k],e2[k])
        
    for i in range(NM):
        imol=mol[i]-1
        for k in range(3):
            image_pbc[k]=pbc[i,k]
        
        if mol_pbc_flag[imol] > 0:
            print(imol,mol_pbc_flag[imol],pbc[i])
            if pbc[i,0]<-508:
                image_pbc[0]+=1024
            if pbc[i,1]<-508:
                image_pbc[1]+=1024
            if pbc[i,2]<-508:
                image_pbc[2]+=1024
        for k in range(3):
            upos[i,k]=pos[i,k]+cell[k,0]*image_pbc[0]+cell[k,1]*image_pbc[1]+cell[k,2]*image_pbc[2]-cell[k,3]
    
    #for k in range(3):
    #    ocell[k,3]=ocell[k,3]-cell[k,3]

    for i in range(NM):
        imol=mol[i]-1
        mol_size[imol] += 1.0
        for k in range(3):
            mol_center[imol,k] += upos[i,k]

    
    for i in range(M):
        for k in range(3):
            mol_center[i,k] /= mol_size[i]
        get_wrapped_pos(mol_center[i],e0,e1,e2,N0,N1,N2,mol_center_wrapped[i])

    for i in range(NM):
        imol=mol[i]-1
        for k in range(3):
            upos[i,k]=upos[i,k]-mol_center[imol,k]+mol_center_wrapped[imol,k]


    org=[cell[0][3],cell[1][3],cell[2][3]]
    er0=np.zeros(3)
    er1=np.zeros(3)
    er2=np.zeros(3)
    # e0*rot
    vec_rot(rot,e0-org,er0)
    # e1*rot
    vec_rot(rot,e1,er1)
    # e2*rot
    vec_rot(rot,e2,er2)
    for k in range(3):
        ocell[k][0]=er0[k]
        ocell[k][1]=er1[k]
        ocell[k][2]=er2[k]
    ocell.vis.rendering_color = (0.0, 0.0, 0.0)


    r_tmp=np.zeros(3)
    r_rot=np.zeros(3)
    for i in range(NM):
        r_tmp=upos[i]
        vec_rot(rot,r_tmp,r_rot)
        upos[i]=r_rot

    if data.particles != None:
        print("There are %i particles with the following properties:" % data.particles.count)
        for property_name in data.particles.keys():
            print("  '%s'" % property_name)
